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

gridC = createCircleMesh("circle",0.5, 1.5, 0.5, 0.4)

gridB = createBoxMeshRev("box", 0.0, 0.5, 1.0, 0.501, 0.1)

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
ε = 200

a, _, Fₑₓₜ, Fᵢₙₜ, K, traction = solver_C2(dh, coord);

if 1 == 1

    using Plots

    contact_dofs = findall(t -> t != 0, traction)

    X_c = coord[contact_dofs, 1]
    tract = traction[contact_dofs]
    ϵᵢⱼₖ = sortperm(X_c)
    X_c = X_c[ϵᵢⱼₖ]
    tract = tract[ϵᵢⱼₖ]

    plot(X_c, tract, legend=false, marker=4, lc=:tomato, mc=:tomato)

end
