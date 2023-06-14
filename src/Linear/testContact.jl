
using Mortar2D, ForwardDiff
using Ferrite, FerriteGmsh, FerriteMeshParser
using LinearSolve, SparseArrays # LinearSolvePardiso
using IterativeSolvers, IncompleteLU    # AlgebraicMultigrid
using SparseDiffTools, Symbolics

include("..//mesh_reader.jl")
include("initLin.jl") # initieras massa skit
include("Contact//contact_help.jl")
include("assemLin.jl")
include("assemElemLin.jl")
include("..//material.jl")
include("..//fem.jl")
#include("runLinContact.jl")
include("run_linear.jl")
include("sensitivitiesLin.jl")

# Create two grids
grid1 = createBoxMesh("box_1", 0.0, 0.0, 1.0, 1.0, 0.05)
grid2 = createBoxMesh("box_2", 0.0, 0.99, 1.0, 1.0, 0.03)

# Merge into one grid
grid_tot = merge_grids(grid1, grid2; tol=1e-6)

# Create dofhandler with displacement field u
dh = DofHandler(grid_tot)



add!(dh, :u, 2)
close!(dh)


# Extract CALFEM-style matrices
coord, enod = getTopology(dh)

register = index_nod_to_grid(dh, coord)

# --------------------------------------------------------------------------- #
# These are useful if Shape functions/gradients defined by Ferrite are needed #
# --------------------------------------------------------------------------- #
ip = Lagrange{2,RefTetrahedron,1}();
qr = QuadratureRule{2,RefTetrahedron}(1);
qr_face = QuadratureRule{1,RefTetrahedron}(1);
cv = CellVectorValues(qr, ip);
fv = FaceVectorValues(qr_face, ip);

# ------------------ #
# Create master sets #
# ------------------ #
addfaceset!(dh.grid, "Γ_master", x -> x[2] ≈ 1.0)
Γm = getfaceset(dh.grid, "Γ_master")

addnodeset!(dh.grid, "nₘ", x -> x[2] ≈ 1.0)
nₘ = getnodeset(dh.grid, "nₘ")

# ----------------- #
# Create slave sets #
# ----------------- #
addfaceset!(dh.grid, "Γ_slave", x -> x[2] ≈ 0.99)
Γs = getfaceset(dh.grid, "Γ_slave")

addnodeset!(dh.grid, "nₛ", x -> x[2] ≈ 0.99)
nₛ = getnodeset(dh.grid, "nₛ")

# Extract all nbr nodes and dofs
contact_dofs = getContactDofs(nₛ, nₘ)
contact_nods = getContactNods(nₛ, nₘ)
order = Dict{Int64,Int64}()
for (i, nod) ∈ enumerate(contact_nods)
    push!(order, nod => i)
end
free_dofs    = setdiff(1:dh.ndofs.x,contact_dofs)

# Define top nodeset for displacement controlled loading
addnodeset!(dh.grid, "Γ_top", x -> x[2] ≈ 1.49)
Γ_top = getnodeset(dh.grid, "Γ_top")


# Define bottom nodeset subject to  u(X) = 0 ∀ X ∈ Γ_bot
addnodeset!(dh.grid, "Γ_bot", x -> x[2] ≈ 0.0)
Γ_bot = getnodeset(dh.grid, "Γ_bot")


# Final preparations for contact
register = getNodeDofs(dh)
X = getX(dh)
coord = getCoordfromX(X)

# Init fictious

coord₀  = deepcopy(coord)
Γ_robin = union(
    getfaceset(dh.grid, "Γ_slave"),
    getfaceset(dh.grid, "Γ_master")
)

# Penalty parameter
ε = 200.0

# boundary conditions for contact analysis
bcdof_top_o, _ = setBCXY(-0.01, dh, Γ_top)
bcdof_bot_o, _ = setBCXY(0.0, dh, Γ_bot)
bcdof_o        = [bcdof_top_o; bcdof_bot_o]


ϵᵢⱼₖ    = sortperm(bcdof_o)
bcdof_o = bcdof_o[ϵᵢⱼₖ]
bcval_o = bcdof_o .* 0.0

# - For Linear solver..
pdofs = bcdof_o
fdofs = setdiff(1:dh.ndofs.x, pdofs)

# Initialize tangents
global K        = create_sparsity_pattern(dh)
global Kψ       = create_sparsity_pattern(dh)
global a        = zeros(dh.ndofs.x)
global Ψ        = zeros(dh.ndofs.x)
global Fᵢₙₜ     = zeros(dh.ndofs.x)
global rc       = zeros(dh.ndofs.x)
global Fₑₓₜ     = zeros(dh.ndofs.x)
global a        = zeros(dh.ndofs.x)
global Δa       = zeros(dh.ndofs.x)
global res      = zeros(dh.ndofs.x)
dh0      = deepcopy(dh)
d        = zeros(size(a))
d       .= 0.0
testvar  = 114
perturbation        = 1e-6
mp       = [210 0.3] # [E ν]
test     = zeros(2)
dFₑₓₜ_dx = similar(K)
C        = zeros(2)
∂rᵤ_∂x   = similar(K)
dr_dd    = similar(K)
∂rψ_∂d   = similar(K)
λᵤ       = similar(a)
λψ       = similar(a)

X_ordered = getXfromCoord(coord)

### "Vanlig"
@time Kc2                       = ForwardDiff.jacobian( u -> contact_residual(X_ordered,u,ε), a);


println("nu kommer det roliga")

###
rc  = contact_residual(X_ordered, a, ε)
rcc = contact_residual_reduced(X_ordered, a[contact_dofs], a[free_dofs], ε)

@time Kc = ForwardDiff.jacobian(u -> contact_residual_reduced(X_ordered, u, a[free_dofs], ε), a[contact_dofs]);

a, _, Fₑₓₜ, Fᵢₙₜ, K, traction = solver_C(dh, coord) # behövs "local" här?





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
    Ψ, _, Kψ, _, λ = fictitious_solver_C(d, dh0, coord₀)
    updateCoords!(dh, Ψ)

    coord          = getCoord(getX(dh), dh)
    register       = getNodeDofs(dh)
    X              = getXfromCoord(coord)
    #coord          = getCoordfromX(X)
    vtk_grid("contact fictious", dh) do vtkfile
        vtk_point_data(vtkfile, dh0, Ψ) # displacement field
    end
    #
    a, _, Fₑₓₜ, Fᵢₙₜ,  K, traction = solver_C(dh, coord) # behövs "local" här?
    #test[pert] = a' * Fₑₓₜ

    test[pert]          = -a[pdofs]' * Fᵢₙₜ[pdofs]
end

∂g_∂u = zeros(size(d))
∂g_∂u[fdofs] = -a[pdofs]' * K[pdofs, fdofs]

#dFₑₓₜ_dx = dFext_dx(dFₑₓₜ_dx, dh, mp, t, a, coord, enod, τ, Γt)
∂rᵤ_∂x = drᵤ_dx_c(∂rᵤ_∂x, dh, mp, t, a, coord, enod, ε)
dr_dd  = drψ(dr_dd, dh0, Ψ, λ, d, Γ_robin, coord₀)

∂g_∂x = zeros(size(d))

∂g_∂x[fdofs] = -a[pdofs]' * ∂rᵤ_∂x[pdofs, fdofs]

solveq!(λᵤ, K', ∂g_∂u, bcdof_o, bcval_o)  # var Fₑₓₜ;
solveq!(λψ, Kψ', ∂g_∂x - ∂rᵤ_∂x' * λᵤ, bcdof_o, bcval_o)

∂g_∂d = -transpose(λψ) * dr_dd
asens = ∂g_∂d[testvar]

numsens = (test[1] - test[2]) / perturbation
numsens / asens

println("numsens: $numsens")
println("asens: $asens")

vtk_grid("contact fictious", dh) do vtkfile
    vtk_point_data(vtkfile, dh, Ψ) # displacement field
end



if 1 == 1

    using Plots

    contact_dofs = findall(t -> t != 0, traction)

    X_c   = coord[contact_dofs, 1]
    tract = traction[contact_dofs]
    ϵᵢⱼₖ  = sortperm(X_c)
    X_c   = X_c[ϵᵢⱼₖ]
    tract = tract[ϵᵢⱼₖ]

    plot(X_c, tract,legend= false, marker= 4, lc= :tomato, mc=:tomato )

end
