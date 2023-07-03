
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


r₀ = 0.5
# Create two grids
grid1 = createCircleMesh("circle", 0.5, 1.5, r₀, 0.15)
#grid1 = createBoxMeshRev("box_2", 0.0, 1.0, 1.0, 0.5, 0.08)
grid2 = createBoxMeshRev("box_1", 0.0, 0.0, 1.0, 1.001, 0.05)

# Merge into one grid
grid_tot = merge_grids(grid1, grid2; tol=1e-6)

grid1 = nothing
grid2 = nothing

# Create dofhandler with displacement field u
global dh = DofHandler(grid_tot)

add!(dh, :u, 2)
close!(dh)


# Extract CALFEM-style matrices
global coord, enod = getTopology(dh)
global register = index_nod_to_grid(dh, coord)


# ------------------ #
# Create master sets #
# ------------------ #
addfaceset!(dh.grid, "Γ_slave", x -> x[2] ≈ 1.001)
global Γs = getfaceset(dh.grid, "Γ_slave")

addnodeset!(dh.grid, "nₛ", x -> x[2] ≈ 1.001)
global nₛ = getnodeset(dh.grid, "nₛ")

# ----------------- #
# Create slave sets #
# ----------------- #

addfaceset!(dh.grid, "Γ_master", x -> ((x[1] - r₀)^2 + (x[2] - 1.5)^2) ≈ r₀^2 && x[2] < 1.5)
global Γm = getfaceset(dh.grid, "Γ_master")

addnodeset!(dh.grid, "nₘ", x -> ((x[1] - r₀)^2 + (x[2] - 1.5)^2) ≈ r₀^2 && x[2] < 1.5)
global nₘ = getnodeset(dh.grid, "nₘ")

#=
addfaceset!(dh.grid, "Γ_master", x -> x[2]≈1.0)
global Γm = getfaceset(dh.grid, "Γ_master")

addnodeset!(dh.grid, "nₘ", x -> x[2] ≈ 1.0)
global nₘ = getnodeset(dh.grid, "nₘ")
=#

# Extract all nbr nodes and dofs
global contact_dofs = getContactDofs(nₛ, nₘ)
global contact_nods = getContactNods(nₛ, nₘ)
global order = Dict{Int64,Int64}()
for (i, nod) ∈ enumerate(contact_nods)
    push!(order, nod => i)
end
global freec_dofs = setdiff(1:dh.ndofs.x, contact_dofs)

# Define top nodeset for displacement controlled loading
addnodeset!(dh.grid, "Γ_top", x -> x[2] ≈ 1.5)
global Γ_top = getnodeset(dh.grid, "Γ_top")

addnodeset!(dh.grid, "n_top", x -> x[2] ≈ 1.5)
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

#for inod in nodx
#   append!(free_d,register[inod,2]*2-1)
#end
global free_d = []
for jnod in n_robin
    append!(free_d, register[jnod, 2])
end

global locked_d = setdiff(1:dh.ndofs.x, free_d)




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

# boundary conditions for contact analysis
bcdof_top_o, _ = setBCXY(-0.01, dh, Γ_top)
bcdof_bot_o, _ = setBCXY(0.0, dh, Γ_bot)
bcdof_o = [bcdof_top_o; bcdof_bot_o]
ϵᵢⱼₖ = sortperm(bcdof_o)
global bcdof_o = bcdof_o[ϵᵢⱼₖ]
global bcval_o = bcdof_o .* 0.0

#bcdof_top_o2, _ = setBCXY_both(0.0, dh, Γ_top)
#bcdof_bot_o2, _ = setBCXY_both(0.0, dh, Γ_bot)
bcdof_top_o2, _ = setBCXY(0.0, dh, Γ_top)
bcdof_bot_o2, _ = setBCXY(0.0, dh, Γ_bot)
bcdof_o2 = [bcdof_top_o2; bcdof_bot_o2]
ϵᵢⱼₖ = sortperm(bcdof_o)
global bcdof_o2 = bcdof_o2[ϵᵢⱼₖ]
global bcval_o2 = bcdof_o2 .* 0.0


# - For Linear solver..
global pdofs = bcdof_o
global fdofs = setdiff(1:dh.ndofs.x, pdofs)


global dh0 = deepcopy(dh)
global d = zeros(size(a))
global d .= 0.0
global ∂rᵤ_∂x = similar(K)
global dr_dd = similar(K)
global ∂rψ_∂d = similar(K)
global λᵤ = similar(a)
global λψ = similar(a)
global λ = 0


global Δ = -0.1
global nloadsteps = 10

include("initOptLin.jl")
ε = 1e5
μ = 1e3
p = 2
testvar  = 191

sens_test = 1
perturbation        = 1e-6
X = getXfromCoord(coord)

remesh = 2

test = [0.0,0.0]

if remesh == 1

    reMeshGrids!(0.05, dh, coord, enod, register, Γs, nₛ, Γm, nₘ, contact_dofs, contact_nods, order, freec_dofs, free_d, locked_d, bcdof_o, bcval_o, d, dh0, coord₀)
    # Initialize tangents
    global K = create_sparsity_pattern(dh) # behövs
    global Kψ = create_sparsity_pattern(dh) # behövs
    global a = zeros(dh.ndofs.x) # behövs
    global Ψ = zeros(dh.ndofs.x) # behövs
    global Fᵢₙₜ = zeros(dh.ndofs.x) # behövs?
    global rc = zeros(dh.ndofs.x) # behövs?
    global Fₑₓₜ = zeros(dh.ndofs.x) # behövs ?
    global a = zeros(dh.ndofs.x) # behövs ?
    global d = zeros(dh.ndofs.x)
    global Δa = zeros(dh.ndofs.x) # behövs inte
    global res = zeros(dh.ndofs.x) # behövs inte
    global ∂rᵤ_∂x = similar(K) # behövs inte om vi har lokal funktion?
    global dr_dd = similar(K) # behövs inte om vi har lokal funktion?
    global ∂rψ_∂d = similar(K) # behövs inte om vi har lokal funktion?
    global ∂g_∂x = zeros(size(a)) # behövs inte om vi har lokal funktion?
    global ∂g_∂u = zeros(size(d)) # behövs inte om vi har lokal funktion?
    global ∂rᵤ_∂x = similar(K) # behövs inte om vi har lokal funktion?
    global λᵤ = similar(a) # behövs inte om vi har lokal funktion?
    global λψ = similar(a) # behövs inte om vi har lokal funktion?
    global λᵥₒₗ = similar(a) # behövs inte om vi har lokal funktion?
    global m = 1 # behöver inte skrivas över
    global n_mma = length(d) # behöver skrivas över
    global epsimin = 0.0000001 # behöver inte skrivas över
    global xval = d[:] # behöver skrivas över
    global xold1 = xval # behöver skrivas över
    global xold2 = xval # behöver skrivas över
    global xmin = -ones(n_mma) / 20 # behöver skrivas över
    global xmax = ones(n_mma) / 20 # behöver skrivas över
    global C = 1000 * ones(m) # behöver inte skrivas över
    global d2 = zeros(m) # behöver inte skrivas över
    global a0 = 1 # behöver inte skrivas över
    global am = zeros(m) # behöver inte skrivas över
    global outeriter = 0 # behöver inte skrivas över
    global kkttol = 0.001 # behöver inte skrivas över
    global changetol = 0.001 # behöver inte skrivas över
    global kktnorm = kkttol + 10 # behöver inte skrivas över
    global outit = 0 # behöver inte skrivas över
    global change = 1 # behöver inte skrivas över
    global xmin[contact_dofs] .= -0.2 # behöver skrivas över
    global xmax[contact_dofs] .= 0.2 # behöver skrivas över
    global xmin[contact_dofs[findall(x -> x % 2 == 0, contact_dofs)]] .= -0.2 # behöver skrivas över
    global xmax[contact_dofs[findall(x -> x % 2 == 0, contact_dofs)]] .= 0.2 # behöver skrivas över
    global low = xmin # behöver skrivas över
    global upp = xmax # behöver skrivas över
    global d .= 0
    #global d[free_d] .= 0.05
    global bcdof_top_o, _ = setBCXY(-0.01, dh, Γ_top)
    global bcdof_bot_o, _ = setBCXY(0.0, dh, Γ_bot)
    global bcdof_o = [bcdof_top_o; bcdof_bot_o]
    ϵᵢⱼₖ = sortperm(bcdof_o)
    global bcdof_o = bcdof_o[ϵᵢⱼₖ]
    global bcval_o = bcdof_o .* 0.0

    #bcdof_top_o2, _ = setBCXY_both(0.0, dh, Γ_top)
    #bcdof_bot_o2, _ = setBCXY_both(0.0, dh, Γ_bot)
    bcdof_top_o2, _ = setBCXY(0.0, dh, Γ_top)
    bcdof_bot_o2, _ = setBCXY(0.0, dh, Γ_bot)

    bcdof_o2 = [bcdof_top_o2; bcdof_bot_o2]
    ϵᵢⱼₖ = sortperm(bcdof_o)
    global bcdof_o2 = bcdof_o2[ϵᵢⱼₖ]
    global bcval_o2 = bcdof_o2 .* 0.0

    global T = zeros(size(a))
    global T[bcdof_bot_o] .= 1.0
    global nloadsteps =  10
end


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

        global coord₀ = getCoord(getX(dh0), dh0) # x₀
        # Check that grid is updated correctly
        #Ψ, _, Kψ, _, λ = fictitious_solver_C(d, dh0, coord₀)
        Ψ, _, Kψ, _, λ = fictitious_solver_with_contact(d, dh0, coord₀, 10)

        global dh = deepcopy(dh0)
        updateCoords!(dh, Ψ) # x₀ + Ψ = x
        global coord = getCoord(getX(dh), dh)

        p = 2
        #global X_ordered      = getXfromCoord(coord) # ta bort och byta då målfunk anropas
        #
        a, _, Fₑₓₜ, Fᵢₙₜ, K, traction = solver_C(dh, coord, Δ/10, 1)
        #test[pert] = a' * Fₑₓₜ

        test[pert]          = -a[pdofs]' * Fᵢₙₜ[pdofs]
        #test[pert] = contact_pnorm(X_ordered, a, ε, p)
    end
    ∂rᵤ_∂x = similar(K)
    ∂rᵤ_∂x = drᵤ_dx_c(∂rᵤ_∂x, dh, mp, t, a, coord, enod, ε)
    dr_dd  = drψ(dr_dd, dh0, Ψ, λ, d, Γ_robin, coord₀)

    ∂g_∂u  = zeros(size(d))
    ∂g_∂x  = zeros(size(d))

    #∂g_∂x = ForwardDiff.gradient(x -> contact_pnorm_ordered(x, a, ε, p), getXinDofOrder(dh, X_ordered, coord))
    #∂g_∂u = ForwardDiff.gradient(u -> contact_pnorm(X_ordered, u, ε, p), a)

    ∂g_∂x[fdofs] = -a[pdofs]' * ∂rᵤ_∂x[pdofs, fdofs]
    ∂g_∂u[fdofs] = -a[pdofs]' * K[pdofs, fdofs]

    solveq!(λᵤ, K', ∂g_∂u, bcdof_o, bcval_o)
    solveq!(λψ, Kψ', ∂g_∂x - ∂rᵤ_∂x' * λᵤ, bcdof_o2, bcval_o2)
    #solveq!(λψ, Kψ', ∂g_∂x - ∂rᵤ_∂x' * λᵤ, bcdof_o, bcval_o)

    ∂g_∂d   = -transpose(λψ) * dr_dd
    asens   =  ∂g_∂d[testvar]

    numsens = (test[1] - test[2]) / perturbation

    println("numsens: $numsens")
    println("asens: $asens")
end



g_test  = 0
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

if 1 == 2
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
