
using Mortar2D, ForwardDiff
using Ferrite, FerriteGmsh, FerriteMeshParser
using LinearSolve, SparseArrays # LinearSolvePardiso
using IterativeSolvers, IncompleteLU    # AlgebraicMultigrid
using SparseDiffTools, Symbolics

include("..//mesh_reader.jl")
#include("initLin.jl") # initieras massa skit
include("Contact//contact_help.jl")
include("assemLin.jl")
include("assemElemLin.jl")
include("..//material.jl")
include("..//fem.jl")
#include("runLinContact.jl")
include("run_linear.jl")
include("sensitivitiesLin.jl")


# Create two grids
grid1 = createCircleMesh("circle", 0.5, 1.5, 0.5, 0.1)
#grid2 = createBoxMeshRev("box_1", 0.0, 0.0, 1.0, 1.001, 0.05)
grid2 = createBoxMeshRev("box_1", -0.5, 0.0, 2.0, 1.001, 0.05)

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
addfaceset!(dh.grid, "Γ_slave", x -> x[2] ≈ 1.001)
Γs = getfaceset(dh.grid, "Γ_slave")

addnodeset!(dh.grid, "nₛ", x -> x[2] ≈ 1.001)
nₛ = getnodeset(dh.grid, "nₛ")

# ----------------- #
# Create slave sets #
# ----------------- #
addfaceset!(dh.grid, "Γ_master", x -> ((x[1] - 0.5)^2 + (x[2] - 1.5)^2) ≈ 0.5^2 && x[2] < 1.5)
Γm = getfaceset(dh.grid, "Γ_master")

addnodeset!(dh.grid, "nₘ", x -> ((x[1] - 0.5)^2 + (x[2] - 1.5)^2) ≈ 0.5^2 && x[2] < 1.5)
nₘ = getnodeset(dh.grid, "nₘ")

# Extract all nbr nodes and dofs
contact_dofs = getContactDofs(nₛ, nₘ)
contact_nods = getContactNods(nₛ, nₘ)
order = Dict{Int64,Int64}()
for (i, nod) ∈ enumerate(contact_nods)
    push!(order, nod => i)
end
freec_dofs    = setdiff(1:dh.ndofs.x,contact_dofs)

# Define top nodeset for displacement controlled loading
addnodeset!(dh.grid, "Γ_top", x -> x[2] ≈ 1.5)
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
mp       = [175 80.769230769230759] # [E ν]
mp₀      = [1.0 1.0]
t        = 1.0
test     = zeros(2)
dFₑₓₜ_dx = similar(K)
C        = zeros(2)
∂rᵤ_∂x   = similar(K)
dr_dd    = similar(K)
∂rψ_∂d   = similar(K)
∂rΨ_∂x   = similar(K)
λᵤ       = similar(a)
λψ       = similar(a)
#traction = similar(contact_dofs)
X_ordered = getXfromCoord(coord)
include("initOptLin.jl")
ε = 5e4
μ = 5e4
p = 2
testvar  = 414

sens_test = 1
perturbation        = 1e-6
X = getXfromCoord(coord)
if sens_test==1

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
        #Ψ, _, Kψ, _, λ = fictitious_solver_C(d, dh0, coord₀)
        Ψ, _, Kψ, _, λ = fictitious_solver_with_contact(d, dh0, coord₀)
        updateCoords!(dh, Ψ)

        coord          = getCoord(getX(dh), dh)
        #register       = getNodeDofs(dh)
        X              = getXfromCoord(coord)
        X_ordered      = getXfromCoord(coord) # ta bort och byta då målfunk anropas
        #
        a, _, Fₑₓₜ, Fᵢₙₜ,  K, traction = solver_C(dh, coord)
        #test[pert] = a' * Fₑₓₜ

        test[pert]          = -a[pdofs]' * Fᵢₙₜ[pdofs]
        #test[pert] = contact_pnorm(X_ordered, a, ε, p)
    end
    ∂rᵤ_∂x = similar(K)
    #∂rΨ_∂x = similar(K)
    ∂rᵤ_∂x = drᵤ_dx_c(∂rᵤ_∂x, dh, mp, t, a, coord, enod, ε)
    #∂rΨ_∂x = drΨ_dx_c(∂rΨ_∂x, dh0, mp₀, t, Ψ, coord₀, enod, λ, d, Γ_robin, μ)
    dr_dd  = drψ(dr_dd, dh0, Ψ, λ, d, Γ_robin, coord₀)

    ∂g_∂u  = zeros(size(d))
    ∂g_∂x  = zeros(size(d))

    #∂g_∂x = ForwardDiff.gradient(x -> contact_pnorm_ordered(x, a, ε, p), getXinDofOrder(dh, X_ordered, coord))
    #∂g_∂u = ForwardDiff.gradient(u -> contact_pnorm(X_ordered, u, ε, p), a)

    ∂g_∂u[fdofs] = -a[pdofs]' * K[pdofs, fdofs]
    ∂g_∂x[fdofs] = -a[pdofs]' * ∂rᵤ_∂x[pdofs, fdofs]

    solveq!(λᵤ, K',  ∂g_∂u, bcdof_o, bcval_o)
    solveq!(λψ, Kψ', ∂g_∂x - ∂rᵤ_∂x' * λᵤ, bcdof_o, bcval_o)
    #solveq!(λψ, Kψ', ∂g_∂x - ∂rᵤ_∂x' * λᵤ, bcdof_o, bcval_o)

    ∂g_∂d   = -transpose(λψ) * dr_dd
    asens   =  ∂g_∂d[testvar]

    numsens = (test[1] - test[2]) / perturbation

    println("numsens: $numsens")
    println("asens: $asens")
end



g_test  = 0
Fψ = similar(Fᵢₙₜ)
#perturbation = 1e-5
#testvar = 30
if g_test == 1
    X_ordered = getXfromCoord(coord₀)
    τ = [1.0 1.0]
    testf = zeros(dh.ndofs.x,2)
    ∂rᵤ_∂x = similar(K)
    for pert in 1:2
        if pert == 1
            #X_ordered[testvar] += perturbation
            #a[testvar] += perturbation
            Ψ[testvar] += perturbation
        else
            #X_ordered[testvar] -= perturbation
            #a[testvar] -= perturbation
            Ψ[testvar] -= perturbation
        end
        #test[pert] = contact_pnorm(X_ordered, a, ε, p)
        #assemGlobal!(K, Fᵢₙₜ, rc, dh, mp, t, a, coord, enod, ε, Γ_top, τ)
        #test[pert] = -a[pdofs]' * Fᵢₙₜ[pdofs]
        assemGlobal!(Kψ, Fψ, dh0, mp₀, t, Ψ, coord₀, enod, λ, d, Γ_robin, μ)
        testf[:, pert] = Fψ
    end
    numsens = (testf[:,1] - testf[:,2]) / perturbation
    #∂rᵤ_∂x  = drᵤ_dx_c(∂rᵤ_∂x, dh, mp, t, a, coord, enod, ε)
    #asens   = -a[pdofs]' * ∂rᵤ_∂x[pdofs, fdofs]
    #@show asens[testvar]
    #asens   = ForwardDiff.gradient(u -> contact_pnorm(X_ordered, u, ε, p), a)
    #asens   = ForwardDiff.gradient(x -> contact_pnorm(x, a, ε, p), X_ordered)
    #∂rΨ_∂x = drΨ_dx_c(∂rΨ_∂x, dh0, mp₀, t, Ψ, coord₀, enod, λ, d, Γ_robin, μ)
    #asens  = ∂rΨ_∂x
    asens = Kψ[:,testvar]
    @show hello = asens[findall(x -> x != 0, numsens)] ./ numsens[findall(x -> x != 0, numsens)]
end

if 1 == 1
    a, _, Fₑₓₜ, Fᵢₙₜ, K, traction = solver_C(dh, coord);
    using Plots
    X_c   = []
    tract = []
    for (key,val) ∈ traction
        append!(X_c,coord[key,1])
        append!(tract,val)
    end
    ϵᵢⱼₖ  = sortperm(X_c)
    tract = tract[ϵᵢⱼₖ]
    X_c   = X_c[ϵᵢⱼₖ]
    plot(X_c, tract,legend= false, marker= 4, lc= :tomato, mc=:tomato )
end
