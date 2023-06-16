using Mortar2D, ForwardDiff
using Ferrite, FerriteGmsh, FerriteMeshParser
using LinearSolve, SparseArrays # LinearSolvePardiso
using IterativeSolvers, IncompleteLU    # AlgebraicMultigrid

include("..//mesh_reader.jl")
include("initLin.jl")
include("Contact//contact_help.jl")
include("assemLin.jl")
include("assemElemLin.jl")
include("..//material.jl")
include("..//fem.jl")
include("run_linear.jl")
include("sensitivitiesLin.jl")

gridC = createCircleMesh("circle",0.5, 1.5, 0.5, 0.1)

gridB = createBoxMeshRev("box", 0.0, 0.5, 1.0, 0.501, 0.01)

grid_tot = merge_grids(gridB, gridC; tol=1e-6)

dh = DofHandler(grid_tot)

add!(dh, :u, 2)
close!(dh)

# Extract CALFEM-style matrices
coord, enod = getTopology(dh)
register = index_nod_to_grid(dh, coord)


# ------------------ #
# Create master sets #
# ------------------ #
addfaceset!(dh.grid, "Γ_master", x -> ((x[1]-0.5)^2 + (x[2]-1.5)^2 ) ≈ 0.5^2 && x[2] < 1.5 )
Γm = getfaceset(dh.grid, "Γ_master")

addnodeset!(dh.grid, "nₘ", x -> ((x[1] - 0.5)^2 + (x[2] - 1.5)^2) ≈ 0.5^2 && x[2] < 1.5)
nₘ = getnodeset(dh.grid, "nₘ")

# ----------------- #
# Create slave sets #
# ----------------- #
addfaceset!(dh.grid, "Γ_slave", x -> x[2] ≈ 1.001)
Γs = getfaceset(dh.grid, "Γ_slave")

addnodeset!(dh.grid, "nₛ", x -> x[2] ≈ 1.001)
nₛ = getnodeset(dh.grid, "nₛ")

# Extract all nbr nodes and dofs
contact_dofs = getContactDofs(nₛ, nₘ)
contact_nods = getContactNods(nₛ, nₘ)

order = Dict{Int64,Int64}()
for (i, nod) ∈ enumerate(contact_nods)
    push!(order, nod => i)
end
freec_dofs = setdiff(1:dh.ndofs.x, contact_dofs)

# Define top nodeset for displacement controlled loading
addnodeset!(dh.grid, "Γ_top", x -> x[2] ≈ 1.5)
Γ_top = getnodeset(dh.grid, "Γ_top")


# Define bottom nodeset subject to  u(X) = 0 ∀ X ∈ Γ_bot
addnodeset!(dh.grid, "Γ_bot", x -> x[2] ≈ 0.5)
Γ_bot = getnodeset(dh.grid, "Γ_bot")


# Final preparations for contact
register = getNodeDofs(dh)
X = getX(dh)
coord = getCoordfromX(X)

# Init fictious

coord₀ = deepcopy(coord)
Γ_robin = union(
    getfaceset(dh.grid, "Γ_slave"),
    getfaceset(dh.grid, "Γ_master")
)

# ---------- #
# Set BCS    #
# ---------- #
# Set bcs - should be moved outside this function
bcdof_top, bcval_top = setBCXY(-0.01, dh, Γ_top)
bcdof_bot, bcval_bot = setBCXY(0.0, dh, Γ_bot)
bcdof = [bcdof_top; bcdof_bot]
bcval = [bcval_top; bcval_bot]

ϵᵢⱼₖ = sortperm(bcdof)
bcdof = bcdof[ϵᵢⱼₖ]
bcval = bcval[ϵᵢⱼₖ]

# - For Linear solver..
pdofs = bcdof
fdofs = setdiff(1:dh.ndofs.x, pdofs)

# Initialize tangents
global K = create_sparsity_pattern(dh)
global Kψ = create_sparsity_pattern(dh)
global a = zeros(dh.ndofs.x)
global Ψ = zeros(dh.ndofs.x)
global Fᵢₙₜ = zeros(dh.ndofs.x)
global rc = zeros(dh.ndofs.x)
global Fₑₓₜ = zeros(dh.ndofs.x)
global a = zeros(dh.ndofs.x)
global Δa = zeros(dh.ndofs.x)
global res = zeros(dh.ndofs.x)


for node in nₘ
    dofs = register[node, :]
    a[dofs] .= 1.0
end

for node in nₛ
    dofs = register[node, :]
    a[dofs] .= 1.0
end


# Penalty parameter
ε = 2000

dh0 = deepcopy(dh)
d = zeros(size(a))
d .= 0.0
testvar = 290
perturbation = 1e-6
mp = [210 0.3] # [E ν]
test = zeros(2)
dFₑₓₜ_dx = similar(K)
C = zeros(2)
∂rᵤ_∂x = similar(K)
dr_dd = similar(K)
∂rψ_∂d = similar(K)
λᵤ = similar(a)
λψ = similar(a)
sens_test = 1

#a, _, Fₑₓₜ, Fᵢₙₜ, K, traction = solver_C2(dh, coord);

if sens_test == 1

    # Test sensitivity
    for pert in 1:2
        if pert == 1
            # perturbera d
            dh = deepcopy(dh0)
            #dh.grid.nodes = deepcopy(dh0.grid.nodes)
            d[testvar] = d[testvar] + perturbation
        else
            # perturbera d och resetta dh
            dh = deepcopy(dh0)
            #dh.grid.nodes = deepcopy(dh0.grid.nodes)
            d[testvar] = d[testvar] - perturbation
        end

        # Check that grid is updated correctly
        Ψ, _, Kψ, _, λ = fictitious_solver_C(d, dh0, coord₀);
        updateCoords!(dh, Ψ)

        coord = getCoord(getX(dh), dh)
        register = getNodeDofs(dh)
        X = getXfromCoord(coord)
        #coord          = getCoordfromX(X)
        vtk_grid("contact fictious", dh) do vtkfile
            vtk_point_data(vtkfile, dh0, Ψ) # displacement field
        end
        #
        a, _, Fₑₓₜ, Fᵢₙₜ, K, traction = solver_C2(dh, coord);
        #test[pert] = a' * Fₑₓₜ

        test[pert] = -a[pdofs]' * Fᵢₙₜ[pdofs]
    end

    ∂g_∂u = zeros(size(d))
    ∂g_∂u[fdofs] = -a[pdofs]' * K[pdofs, fdofs]

    #dFₑₓₜ_dx = dFext_dx(dFₑₓₜ_dx, dh, mp, t, a, coord, enod, τ, Γt)
    ∂rᵤ_∂x = drᵤ_dx_c(∂rᵤ_∂x, dh, mp, t, a, coord, enod, ε)
    dr_dd = drψ(dr_dd, dh0, Ψ, λ, d, Γ_robin, coord₀)

    ∂g_∂x = zeros(size(d))

    ∂g_∂x[fdofs] = -a[pdofs]' * ∂rᵤ_∂x[pdofs, fdofs]

    solveq!(λᵤ, K', ∂g_∂u, bcdof, bcval)  # var Fₑₓₜ;
    solveq!(λψ, Kψ', ∂g_∂x - ∂rᵤ_∂x' * λᵤ, bcdof, bcval)

    ∂g_∂d = -transpose(λψ) * dr_dd
    asens = ∂g_∂d[testvar]

    numsens = (test[1] - test[2]) / perturbation
    numsens / asens

    println("numsens: $numsens")
    println("asens: $asens")

    vtk_grid("contact fictious", dh) do vtkfile
        vtk_point_data(vtkfile, dh, Ψ) # displacement field
    end

end

a, _, Fₑₓₜ, Fᵢₙₜ, K, traction = solver_C2(dh, coord);

if 1 == 1

    using Plots
    X_c = []
    tract = []
    for (key, val) ∈ traction
        append!(X_c, coord[key, 1])
        append!(tract, val)
    end

    ϵᵢⱼₖ = sortperm(X_c)
    tract = tract[ϵᵢⱼₖ]
    X_c = X_c[ϵᵢⱼₖ]

    plot(X_c, tract, legend=false, marker=4, lc=:tomato, mc=:tomato)

end
