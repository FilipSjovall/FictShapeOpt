using Mortar2D, ForwardDiff
using Ferrite, FerriteGmsh, FerriteMeshParser
using LinearSolve, SparseArrays # LinearSolvePardiso
using IterativeSolvers, IncompleteLU    # AlgebraicMultigrid
#using SparseDiffTools
using Plots, Printf

include("..//mesh_reader.jl")
include("initLin.jl")
include("Contact//contact_help.jl")
include("assemLin.jl")
include("assemElemLin.jl")
include("..//material.jl")
include("..//fem.jl")
include("run_linear.jl")
include("sensitivitiesLin.jl")
include("..//mma.jl")

r  = 0.49
y₀ = 0.999
h  = 0.1

grid1   = createBoxMeshRounded("rounded", r, h)
grid2   = createBoxMeshRounded_Flipped("rounded_and_flipped", r, 1.0, 1.49, h)

dh2 = DofHandler(grid2)
add!(dh2, :u, 2)
close!(dh2)

dh1 = DofHandler(grid1)
add!(dh1, :u, 2)
close!(dh1)

slaves  = getSlaveCoord_remeshed(dh1)
masters = getMasterCoord_remeshed(dh2)
# Merge into one grid
grid_tot = merge_grids(grid1, grid2; tol=1e-6)
grid1    = nothing
grid2    = nothing
# Create dofhandler with displacement field u
global dh = DofHandler(grid_tot)
add!(dh, :u, 2)
close!(dh)
# Extract CALFEM-style matrices
global coord, enod = getTopology(dh)
global register = index_nod_to_grid(dh, coord)
# Get master and slave sets for combined grid
addfaceset!(dh.grid, "Γ_slave", x -> x ∈ eachrow(slaves))
global Γs = getfaceset(dh.grid, "Γ_slave")
addnodeset!(dh.grid, "nₛ", x -> x ∈ eachrow(slaves))
global nₛ = getnodeset(dh.grid, "nₛ")

addfaceset!(dh.grid, "Γ_master", x -> x ∈ eachrow(masters))
global Γm = getfaceset(dh.grid, "Γ_master")
addnodeset!(dh.grid, "nₘ", x -> x ∈ eachrow(masters))
global nₘ = getnodeset(dh.grid, "nₘ")

# Extract all nbr nodes and dofs
global contact_dofs = getContactDofs(nₛ, nₘ)
global contact_nods = getContactNods(nₛ, nₘ)
global order = Dict{Int64,Int64}()
for (i, nod) ∈ enumerate(contact_nods)
    push!(order, nod => i)
end
global freec_dofs = setdiff(1:dh.ndofs.x, contact_dofs)
# Define top nodeset for displacement controlled loading
addnodeset!(dh.grid, "Γ_top", x -> x[2] ≈ 2.0)
global Γ_top = getnodeset(dh.grid, "Γ_top")

addnodeset!(dh.grid, "n_top", x -> x[2] ≈ 2.0)
global n_top = getnodeset(dh.grid, "n_top")

# Define bottom nodeset subject to  u(X) = 0 ∀ X ∈ Γ_bot
addnodeset!(dh.grid, "Γ_bot", x -> x[2] ≈ 0.0)
global Γ_bot = getnodeset(dh.grid, "Γ_bot")

addnodeset!(dh.grid, "n_bot", x -> x[2] ≈ 0.0)
global n_bot = getnodeset(dh.grid, "n_bot")

# Final preparations for contact
global register = getNodeDofs(dh)
global X = getX(dh)
global coord = getCoordfromX(X)

# Init fictious

global coord₀ = deepcopy(coord)
global Γ_robin = union(
    getfaceset(dh.grid, "Γ_slave"),
    getfaceset(dh.grid, "Γ_master")
)
global n_robin = union(
    getnodeset(dh.grid, "nₛ"),
    getnodeset(dh.grid, "nₘ")
)

global K = create_sparsity_pattern(dh)
global a = zeros(dh.ndofs.x)
global Fᵢₙₜ = zeros(dh.ndofs.x)
global a = zeros(dh.ndofs.x)
global Δa = zeros(dh.ndofs.x)

# boundary conditions for contact analysis
bcdof_top_o, _ = setBCXY_both(-0.01, dh, Γ_top)
bcdof_bot_o, _ = setBCXY_both(0.0, dh, Γ_bot)
bcdof_o = [bcdof_top_o; bcdof_bot_o]
ϵᵢⱼₖ = sortperm(bcdof_o)
global bcdof_o = bcdof_o[ϵᵢⱼₖ]
global bcval_o = bcdof_o .* 0.0

global Δ = -0.1
global nloadsteps = 10
global ε = 1e6

a, _, Fₑₓₜ, Fᵢₙₜ, K, traction = solver_C(dh, coord, Δ, nloadsteps)
