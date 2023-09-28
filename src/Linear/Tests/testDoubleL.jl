using Mortar2D, ForwardDiff, Ferrite, FerriteGmsh, FerriteMeshParser
using LinearSolve, SparseArrays, IterativeSolvers, IncompleteLU
using SparseDiffTools, Plots, Printf, JLD2, Statistics

include("..//mesh_reader.jl")
include("Contact//contact_help.jl")
include("assemLin.jl")
include("assemElemLin.jl")
include("..//material.jl")
include("..//fem.jl")
include("run_linear.jl")
include("sensitivitiesLin.jl")
include("..//mma.jl")


xl = 0.0
yl = 0.0
xr = -0.4999
yr = 1.4
Δx = 1.0
Δy = 1.0
th = 0.25
r1 = 0.05
r2 = 0.075

# # # # # # # # # #
# Finite element  #
# # # # # # # # # #
ip      = Lagrange{2,RefTetrahedron,1}()
qr      = QuadratureRule{2,RefTetrahedron}(1)
qr_face = QuadratureRule{1,RefTetrahedron}(1)
cv      = CellVectorValues(qr, ip)
fv      = FaceVectorValues(qr_face, ip)

# # # # # # # # #
# Create grids  #
# # # # # # # # #
grid1    = createLMesh("mesh_1", xl, yl, Δx, Δy, th, r1, r2, 0.05);
Γ_1      = getBoundarySet(grid1);
grid2    = createLMeshRev("mesh_2", xr, yr, Δx, Δy, th, r1, r2, 0.05);
Γ_2      = getBoundarySet(grid2);
grid_tot = merge_grids(grid1, grid2; tol=1e-6);
grid1    = nothing;
grid2    = nothing;
# Create dofhandler with displacement field u
global dh = DofHandler(grid_tot);
add!(dh, :u, 2);
close!(dh);

# Extract CALFEM-style matrices
global coord, enod = getTopology(dh);
global register = index_nod_to_grid(dh, coord);

# Exrtact full boundary
Γ_all = Ferrite.__collect_boundary_faces(dh.grid);
addfaceset!(dh.grid, "Γ_all", Γ_all);
Γ_all = getfaceset(dh.grid, "Γ_all");

n_all = getBoundarySet(dh.grid, Γ_all);
addnodeset!(dh.grid, "n_all", n_all);

Γ_all_dofs = Vector{Int64}()
# --------------
# Master
addfaceset!(dh.grid, "Γ_master", x -> x ∈ Γ_1);
Γm = getfaceset(dh.grid, "Γ_master");
Γm = intersect(Γm, Γ_all);

nₘ = getBoundarySet(dh.grid, Γm);
addnodeset!(dh.grid, "nₘ", nₘ);

# ---------------
# Slave
addfaceset!(dh.grid, "Γ_slave", x -> x ∈ Γ_2);
Γs = getfaceset(dh.grid, "Γ_slave");
Γs = intersect(Γs, Γ_all);

nₛ = getBoundarySet(dh.grid, Γs)
addnodeset!(dh.grid, "nₛ", nₛ)

# ---------------
# Displacement bc boundary u(x) = Δ ∀ x ∈ Γ_Δ
addfaceset!(dh.grid, "Γ_right", x -> x[1] ≈ xl + Δx)
Γ_right = getfaceset(dh.grid, "Γ_right")

addnodeset!(dh.grid, "n_right", x -> x[1] ≈ xl + Δx)
n_right = getnodeset(dh.grid, "n_right")

# --------------
# Displacement bc boundary u(x) = 0 ∀ x ∈ Γ_0
addfaceset!(dh.grid, "Γ_left", x -> x[1] ≈ xr)
Γ_left = getfaceset(dh.grid, "Γ_left")

addnodeset!(dh.grid, "n_left", x -> x[1] ≈ xr)
n_left = getnodeset(dh.grid, "n_left")

# --------------
# bottom
addfaceset!(dh.grid, "Γ_bot", x -> x[2] ≈ yl)
Γ_bot = getfaceset(dh.grid, "Γ_bot")

addnodeset!(dh.grid, "n_bot", x -> x[2] ≈ yl)
n_bot = getnodeset(dh.grid, "n_bot")

# --------------
# Top
addfaceset!(dh.grid, "Γ_top", x -> x[2] ≈ yr)
Γ_top = getfaceset(dh.grid, "Γ_top")

addnodeset!(dh.grid, "n_top", x -> x[2] ≈ yr)
n_top = getnodeset(dh.grid, "n_top")

# ---------------
# Design boundaries
Γ_robin = setdiff(Γ_all, union(Γ_left, Γ_right, Γ_bot, Γ_top))
addfaceset!(dh.grid, "Γ_robin", Γ_robin)

n_robin = getBoundarySet(dh.grid, Γ_robin)
addnodeset!(dh.grid, "n_robin", n_robin)


# # # # # # # # # # # # #
# Collect contact dofs  #
# # # # # # # # # # # # #
global contact_dofs = getContactDofs(nₛ, nₘ)
global contact_nods = getContactNods(nₛ, nₘ)
global order = Dict{Int64,Int64}()
for (i, nod) ∈ enumerate(contact_nods)
    push!(order, nod => i)
end
global freec_dofs = setdiff(1:dh.ndofs.x, contact_dofs)

# Final preparations for contact
global register = getNodeDofs(dh)
global X = getX(dh)
global coord = getCoordfromX(X)

# # # # # # # # #
# Init fictious #
# # # # # # # # #
global coord₀ = deepcopy(coord)

# # # # # # # # # # # #
# Collect design dofs #
# # # # # # # # # # # #
global free_d = []
for jnod in n_robin
    append!(free_d, register[jnod, 1])
    append!(free_d, register[jnod, 2])
end
#global locked_d = setdiff(1:dh.ndofs.x, free_d)
global locked_d = Vector{Int64}()
for n ∈ n_left
    push!(locked_d, register[n, 1])
    push!(locked_d, register[n, 2])
end
for n ∈ n_right
    push!(locked_d, register[n, 1])
    push!(locked_d, register[n, 2])
end

# Initialize tangents
global K = create_sparsity_pattern(dh)
global Kψ = create_sparsity_pattern(dh)
global a = zeros(dh.ndofs.x)
global d = zeros(dh.ndofs.x)
global Ψ = zeros(dh.ndofs.x)
global Fᵢₙₜ = zeros(dh.ndofs.x)
global rc = zeros(dh.ndofs.x)
global Fₑₓₜ = zeros(dh.ndofs.x)
global a = zeros(dh.ndofs.x)
global Δa = zeros(dh.ndofs.x)
global res = zeros(dh.ndofs.x)
global dr_dd = similar(K)
global ∂rψ_∂d = similar(K)
global ∂g_∂x = zeros(size(a)) # behövs inte om vi har lokal funktion?
global ∂g_∂u = zeros(size(d)) # behövs inte om vi har lokal funktion?
global λᵤ = similar(a)
global λψ = similar(a)

global Δ = 0.1
global nloadsteps = 10

include("initOptLinHook.jl")

# ------------------- #
# Boundary conditions #
# ------------------- #
bcdof_left, _  = setBCX(0.0, dh, n_left)
bcdof_right, _ = setBCX(0.0, dh, n_right)
bcdof_bot, _   = setBCY(0.0, dh, n_bot)
bcdof_top, _   = setBCY(0.0, dh, n_top)
# Kommentera om låsning i y-led
# bcdof_top      = Vector{Int64}();
# bcdof_bot      = Vector{Int64}();

bcdofs_opt        = [bcdof_left; bcdof_right; bcdof_bot; bcdof_top];
ϵᵢⱼₖ             = sortperm(bcdofs_opt)
global bcdofs_opt = bcdofs_opt[ϵᵢⱼₖ]
global bcval_opt  = bcdofs_opt .* 0.0


test         = [0.0, 0.0]
testvar      = register[13,1] #443 # 286 # 763 # 362
perturbation = -1e-6
global T     = zeros(size(a))
global T[bcdof_right[isodd.(bcdof_right)]] .=  1.0
#global T[bcdof_left[isodd.(bcdof_left)]]   .= -1.0

dh0 = deepcopy(dh)

# Test sensitivity
for pert in 1:2
    if pert == 1 #  perturbera d

        dh = deepcopy(dh0)
        global d[testvar] = d[testvar] + perturbation
    else # perturbera d och resetta dh
        dh = deepcopy(dh0)
        global d[testvar] = d[testvar] - perturbation
    end

    global coord₀ = getCoord(getX(dh0), dh0) # x₀
    # Check that grid is updated correctly
    Ψ, _, Kψ, _, λ = fictitious_solver_with_contact_hook(d, dh0, coord₀, 20)

    global dh = deepcopy(dh0)
    updateCoords!(dh, Ψ) # x₀ + Ψ = x
    global coord = getCoord(getX(dh), dh)

    #global X_ordered      = getXfromCoord(coord) # ta bort och byta då målfunk anropas
    #
    global ε = 1e4
    a, _, Fₑₓₜ, Fᵢₙₜ, K, traction = solver_C_hook(dh, coord, Δ, 10)
    #test[pert] = a' * Fₑₓₜ

    test[pert] = -T' * Fᵢₙₜ
    #test[pert] = contact_pnorm(X_ordered, a, ε, p)
end

∂rᵤ_∂x = similar(K)
∂rᵤ_∂x = drᵤ_dx_c(∂rᵤ_∂x, dh, mp, t, a, coord, enod, ε)
dr_dd = drψ(dr_dd, dh0, Ψ, λ, d, Γ_robin, coord₀)

∂g_∂u = zeros(size(d))
∂g_∂x = zeros(size(d))

#assemGlobal!(K, Fᵢₙₜ, dh, mp, t, a, coord, enod)


∂g_∂x = -T' * ∂rᵤ_∂x#drᵤ_dx(dr, dh, mp, t, a, coord, enod) ## stämmer? innehåller kontaktkänslighet men dessa träffar bara kontaktdofs som inte är kopplade till bc.
∂g_∂u = -T' * K

solveq!(λᵤ, K', ∂g_∂u, bcdofs_opt, bcval_opt)
solveq!(λψ, Kψ', ∂g_∂x' - ∂rᵤ_∂x' * λᵤ, bcdofs_opt, bcval_opt)

∂g_∂d = -transpose(λψ) * dr_dd
asens = ∂g_∂d[testvar]

numsens = (test[1] - test[2]) / perturbation

println("numsens: $numsens")
println("asens:   $asens")

#=
sorted = similar(∂g_∂d)

for (idx,idy) ∈ enumerate(register[:,1])
   @show idx #idy
   sorted[register[idx,1]] = ∂g_∂d[idx*2-1]
   sorted[register[idx,2]] = ∂g_∂d[idx*2]
end

postprocess_opt(sorted', dh0, "test")
=#
