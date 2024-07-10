using Pkg
Pkg.activate()
# kolla Pkg.status() vid problem / jämför med att bara starta julia i en terminal
using ForwardDiff, Ferrite, FerriteGmsh, FerriteMeshParser
using LinearSolve, SparseArrays, IterativeSolvers, IncompleteLU
using SparseDiffTools, Plots, Printf, JLD2, Statistics, AlgebraicMultigrid
# kan behöva köra export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib
plotlyjs()
#
include("Contact//Mortar2D//Mortar2D.jl")
include("..//mesh_reader.jl")
include("Contact//contact_help.jl")
include("assemLin.jl")
include("assemElemLin.jl")
include("..//material.jl")
include("..//fem.jl")
include("run_linear.jl")
include("sensitivitiesLin.jl")
include("..//mma.jl")
# ------------------- #
# Geometry parameters #
# ------------------- #
# - Block - #
th = 0.1
x₁ = 0.0
y₁ = 0.25001
Δx = 0.5
Δy = 0.1
# - Seal - #
x₀ = 0.0
y₀ = 0.0
B  = 0.25
b  = 0.1
Δl = (Δx - B)  #0.05
H  = 0.15
r = 0.025
# grid size
h = 0.075
# # # # # # # # # #
# Finite element  #
# # # # # # # # # #
ip = Lagrange{2,RefTetrahedron,1}()
qr = QuadratureRule{2,RefTetrahedron}(3)
qr_face = QuadratureRule{1,RefTetrahedron}(2)
cv = CellVectorValues(qr, ip)
fv = FaceVectorValues(qr_face, ip)
# # # # # # # # #
# Create grids  #
# # # # # # # # #
grid1 = createQuarterLabyrinthMeshRounded("mesh_1", x₀, y₀, th, B, b, Δl, H, r, h);
Γ_1 = getBoundarySet(grid1);
grid2 = createBoxMeshRev("mesh_2", x₁, y₁, Δx, Δy, h/2);
Γ_2 = getBoundarySet(grid2);
grid_tot = merge_grids(grid1, grid2; tol=1e-8);
grid1 = nothing;
grid2 = nothing;
# ------------------------------------------- #
# Create dofhandler with displacement field u #
# ------------------------------------------- #
global dh = DofHandler(grid_tot);
add!(dh, :u, 2);
close!(dh);
# Extract CALFEM-style matrices
global coord, enod = getTopology(dh);
global register = index_nod_to_grid(dh, coord);
# - - - - - - - - #
# Create cellsets #
# - - - - - - - - #
top_mesh = addcellset!(dh.grid, "top mesh", x -> x[2] > th + H)
bot_mesh = addcellset!(dh.grid, "bot mesh", x -> x[2] ≤ th + H)
# Exrtact full boundary
Γ_all = Ferrite.__collect_boundary_faces(dh.grid);
addfaceset!(dh.grid, "Γ_all", Γ_all);
Γ_all = getfaceset(dh.grid, "Γ_all");
#
n_all = getBoundarySet(dh.grid, Γ_all);
addnodeset!(dh.grid, "n_all", n_all);
#
Γ_all_dofs = Vector{Int64}()
# ------ #
# Master #
# ------ #
addfaceset!(dh.grid, "Γ_master", x -> x[2] ≈ y₁);
Γm = getfaceset(dh.grid, "Γ_master");
Γm = intersect(Γm, Γ_all);
#
nₘ = getBoundarySet(dh.grid, Γm);
addnodeset!(dh.grid, "nₘ", nₘ);
#
# ----- #
# Slave #
# ----- #
addfaceset!(dh.grid, "Γ_slave", x ->  x ∈ Γ_1 );
Γs = getfaceset(dh.grid, "Γ_slave");
Γs = intersect(Γs, Γ_all);
#
global nₛ = getBoundarySet(dh.grid, Γs)
addnodeset!(dh.grid, "nₛ", nₛ)

# ------ #
# bottom #
# ------ #
addfaceset!(dh.grid, "Γ_bot", x -> x[2] ≈ y₀)
Γ_bot = getfaceset(dh.grid, "Γ_bot")

addnodeset!(dh.grid, "n_bot", x -> x[2] ≈ y₀)
n_bot = getnodeset(dh.grid, "n_bot")
# --- #
# Top #
# --- #
addfaceset!(dh.grid, "Γ_top", x -> x[2] ≈ y₁ + Δy)
Γ_top = getfaceset(dh.grid, "Γ_top")

addnodeset!(dh.grid, "n_top", x -> x[2] ≈ y₁ + Δy)
n_top = getnodeset(dh.grid, "n_top")

# ----------------------------- #
# left and right sides of block #
# ----------------------------- #
#addnodeset!(dh.grid, "n_lr", x -> ( x[2]≥ y₁ && (x[1] ≈ x₁ || x[1]≈ x₁ + Δx) ))
addnodeset!(dh.grid, "n_lr", x -> ( x[2]≥ y₁ && x[1] ≈ x₁  ))
nₗᵣ = getnodeset(dh.grid, "n_lr")

addfaceset!(dh.grid, "Γ_lr", x -> ( x[2]≥ y₁ && x[1] ≈ x₁  ))
Γ_lr = getfaceset(dh.grid, "Γ_lr")

# ------------------------------ #
# Middle nodes on top and bottom #
# ------------------------------ #
addnodeset!(dh.grid,"n_sym", x->x[1] ≈ 0.5)
n_sym = getnodeset(dh.grid, "n_sym")

addfaceset!(dh.grid, "Γ_sym", x->x[1] ≈ 0.5)
Γ_sym = getfaceset(dh.grid, "Γ_sym")

# ----------------- #
# Design boundaries #
# ----------------- #
Γ_robin = setdiff(Γ_all, union(Γ_top, Γ_bot, Γm, Γ_sym, Γ_lr))
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
global ∂g_∂x = zeros(size(a))
global ∂g_∂u = zeros(size(d))
global ∂g₂_∂x = zeros(size(a))
global ∂g₂_∂u = zeros(size(d))
global λᵤ = similar(a)
global λψ = similar(a)
global Δ = -0.025
global nloadsteps = 10
global g = 0.0

# # # # # # # # # # # # # # # #
# Init optimization variables #
# # # # # # # # # # # # # # # #
include("initLab.jl")
# ------------------- #
# Boundary conditions #
# ------------------- #
bcdof_bot, _ = setBCY(0.0, dh, n_bot)
bcdof_top, _ = setBCY(Δ, dh, n_top)

bcdof_right, _ = setBCX(0.0, dh, n_sym)


# - - - - - - - - #
# Lås master dofs #
# - - - - - - - - #
bcdof_contact, _ = setBCXY_both(0.0, dh, nₘ) # union(n,n,n) om flera set skall slås samman
bcdofs_opt = [bcdof_bot; bcdof_top; bcdof_contact; bcdof_right];
ϵᵢⱼₖ      = sortperm(bcdofs_opt)
global bcdofs_opt  = bcdofs_opt[ϵᵢⱼₖ]
global bcval_opt   = bcdofs_opt .* 0.0
global asy_counter = zeros(dh.ndofs.x, 300)

global low_hist = zeros(length(d), 300)
global upp_hist = zeros(length(d), 300)
global d_hist2  = zeros(length(d), 300)

global low_hist = zeros(length(d), 300)
global upp_hist = zeros(length(d), 300)
global d_hist2  = zeros(length(d), 300)

test         = [0.0, 0.0]
testvar      = 17#register[19,1] # register[21,1] #443 # 286 # 763 # 362
perturbation = 1e-5

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

    global nloadsteps = 10 #10
    global μ = 1e4 #1e3
    global coord₀ = getCoord(getX(dh0), dh0) # x₀
    Ψ, _, Kψ, _, λ = fictitious_solver_with_contact_lab(d, dh0, coord₀, nloadsteps)

    global dh = deepcopy(dh0)
    updateCoords!(dh, Ψ) # x₀ + Ψ = x
    global coord = getCoord(getX(dh), dh)

    #global X_ordered      = getXfromCoord(coord) # ta bort och byta då målfunk anropas
    #
    global nloadsteps = 10
    global ε = 1e5
    a, _, Fₑₓₜ, Fᵢₙₜ, K = solver_Lab(dh, coord, Δ, nloadsteps)
    #test[pert] = a' * Fₑₓₜ

    # test[pert] = -T' * Fᵢₙₜ
    p = 1
    X_ordered = getXfromCoord(coord)
    test[pert] = 1-contact_area(X_ordered, a, 1.0)./5
end

∂rᵤ_∂x = similar(K)
∂rᵤ_∂x = drᵤ_dx_c(∂rᵤ_∂x, dh, t, a, coord, enod, ε, mp₁, mp₂)
dr_dd = drψ(dr_dd, dh0, Ψ, λ, d, Γ_robin, coord₀)

∂g_∂u = zeros(size(d))
∂g_∂x = zeros(size(d))

#assemGlobal!(K, Fᵢₙₜ, dh, mp, t, a, coord, enod)


# ∂g_∂x = -T' * ∂rᵤ_∂x#drᵤ_dx(dr, dh, mp, t, a, coord, enod) ## stämmer? innehåller kontaktkänslighet men dessa träffar bara kontaktdofs som inte är kopplade till bc.
# ∂g_∂u = -T' * K
X_ordered = getXfromCoord(coord)
p = 1
∂g_∂x = ForwardDiff.gradient(x -> contact_area_ordered(x, a, 1.0, ), getXinDofOrder(dh, X_ordered, coord))
∂g_∂u = ForwardDiff.gradient(u -> contact_area(X_ordered, u, 1.0, ), a)

solveq!(λᵤ, K', -∂g_∂u./5, bcdofs_opt, bcval_opt.*0)
solveq!(λψ, Kψ', -∂g_∂x./5 - ∂rᵤ_∂x' * λᵤ, bcdofs_opt, bcval_opt.*0)

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
