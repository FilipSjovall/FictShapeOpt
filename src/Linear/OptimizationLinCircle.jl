using Mortar2D, ForwardDiff
using Ferrite, FerriteGmsh, FerriteMeshParser
using LinearSolve, SparseArrays # LinearSolvePardiso
using IterativeSolvers, IncompleteLU    # AlgebraicMultigrid
using SparseDiffTools
using Plots

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

include("..//mma.jl")

#include("initLin.jl")


# Create two grids
grid1 = createCircleMesh("circle", 0.5, 1.5, 0.5, 0.15)
grid2 = createBoxMeshRev("box_1", 0.0, 0.0, 1.0, 1.001, 0.05)

# Merge into one grid
grid_tot = merge_grids(grid1, grid2; tol=1e-6)

# Create dofhandler with displacement field u
dh = DofHandler(grid_tot)



add!(dh, :u, 2)
close!(dh)


# Extract CALFEM-style matrices
coord, enod = getTopology(dh)
register = index_nod_to_grid(dh, coord)


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

coord₀ = deepcopy(coord)
Γ_robin = union(
    getfaceset(dh.grid, "Γ_slave"),
    getfaceset(dh.grid, "Γ_master")
)
n_robin = union(
    getnodeset(dh.grid, "nₛ"),
    getnodeset(dh.grid, "nₘ")
)

#for inod in nodx
#   append!(free_d,register[inod,2]*2-1)
#end
free_d = []
for jnod in n_robin
    append!(free_d, register[jnod, 2] )
end

locked_d = setdiff(1:dh.ndofs.x,free_d)




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
bcdof_o = bcdof_o[ϵᵢⱼₖ]
bcval_o = bcdof_o .* 0.0

# - For Linear solver..
pdofs = bcdof_o
fdofs = setdiff(1:dh.ndofs.x, pdofs)


dh0 = deepcopy(dh)
d   = zeros(size(a))
d .= 0.0
global ∂rᵤ_∂x = similar(K)
dr_dd = similar(K)
∂rψ_∂d = similar(K)
λᵤ = similar(a)
λψ = similar(a)


include("initOptLin.jl")


function Optimize(dh)

    # Flytta allt nedan till init_opt?
        dh0 = deepcopy(dh)
        λψ   = similar(a)
        λᵤ   = similar(a)
        λᵥₒₗ = similar(a)
        Vₘₐₓ = 1.2 * volume(dh, coord, enod)
        global ε    = 1e5
        global μ    = 1e4
        l    = similar(a)
        l   .= 0.5
        tol  = 1e-6
        OptIter = 0
        coord₀  = coord
        global g_hist = zeros(200)
        global v_hist = zeros(200)
        global T = zeros(size(a))
        T[bcdof_bot_o] .= 1.0
    #
    while kktnorm > tol || OptIter < 3 #&& OptIter < 50

        # # # # # # # # # # # # # #
        #       Definitions       #
        # # # # # # # # # # # # # #
            global d
            global Ψ
            global a
            global Fₑₓₜ
            global K
            global Kψ
            global ∂rᵤ_∂x
            global dr_dd
            global ∂rψ_∂d
            global mp
            global mp₀
            global t
            global m
            global n
            global epsimin
            global xval
            global xold1
            global xold2
            global xmin
            global xmax
            global low#A() = A(Float64[],[]).
            global C
            global d2
            global a0
            global outeriter
            global am
            global kkttol
            global changetol
            global kktnorm
            global outit
            global change
            global λ
            global g_ini
            global pdofs = bcdof_o
            global fdofs = setdiff(1:length(a), pdofs)
            global locked_d = setdiff(1:length(a),free_d)
            global low
            global upp
            global traction
        # # # # # # # # # # # # # #
        OptIter += 1

        # # # # # # # # # # # # # #
        # Fictitious equillibrium #
        # # # # # # # # # # # # # #
        coord₀ = getCoord(getX(dh0), dh0) # x₀
        #Ψ, _, Kψ, _, λ = fictitious_solver_C(d, dh0, coord₀)
        Ψ, _, Kψ, _, λ = fictitious_solver_with_contact(d, dh0, coord₀)

        # # # # # #
        # Filter  #
        # # # # # #
        dh    = deepcopy(dh0)
        updateCoords!(dh, Ψ) # x₀ + Ψ = x
        coord = getCoord(getX(dh), dh)

        # # # # # # # # #
        # Equillibrium  #
        # # # # # # # # #
        a, _, Fₑₓₜ, Fᵢₙₜ, K, traction = solver_C(dh, coord)

        # # # # # # #
        # Objective #
        # # # # # # #
        # Max reaction force
        g            = -T' * Fᵢₙₜ
        ∂g_∂x = -(T' * ∂rᵤ_∂x)'
        ∂g_∂u = -(T' * K)'
        # Compliance
        #g = -a[pdofs]' * Fᵢₙₜ[pdofs]
        #∂g_∂x[fdofs] = -a[pdofs]' * ∂rᵤ_∂x[pdofs, fdofs]
        #∂g_∂u[fdofs] = -a[pdofs]' * K[pdofs, fdofs]
        # Max/Min λ
        #p = 2
        #X_ordered = getXfromCoord(coord)
        #g         = contact_pnorm(X_ordered, a, ε, p)
        #∂g_∂x     = ForwardDiff.gradient(x -> contact_pnorm_ordered(x, a, ε, p), getXinDofOrder(dh, X_ordered, coord))
        #∂g_∂u     = ForwardDiff.gradient(u -> contact_pnorm(X_ordered, u, ε, p), a)

        # # # # # # # # #
        # Sensitivities #
        # # # # # # # # #
        ∂rᵤ_∂x = similar(K)
        ∂rᵤ_∂x       = drᵤ_dx_c(∂rᵤ_∂x, dh, mp, t, a, coord, enod, ε)
        dr_dd        = drψ(dr_dd, dh0, Ψ, λ, d, Γ_robin, coord₀)

        # # # # # # #
        # Adjoints  #
        # # # # # # #
        solveq!(λᵤ, K', ∂g_∂u, bcdof_o, bcval_o)
        solveq!(λψ, Kψ', ∂g_∂x - ∂rᵤ_∂x' * λᵤ, bcdof_o, bcval_o)

        # # # # # # # # # # #
        # Full sensitivity  #
        # # # # # # # # # # #
        ∂g_∂d = (-transpose(λψ) * dr_dd)'
        ∂g_∂d[locked_d] .= 0.0 # fulfix?

        # # # # # # # # # # #
        # Volume constraint #
        # # # # # # # # # # #
        g₁    = volume(dh,coord,enod) / Vₘₐₓ - 1
        ∂Ω_∂x = volume_sens(dh,coord)
        solveq!(λᵥₒₗ, Kψ, ∂Ω_∂x, bcdof_o, bcval_o.*0);
        ∂Ω∂d  = -transpose(λᵥₒₗ)*dr_dd ./ Vₘₐₓ;
        ∂Ω∂d[locked_d] .= 0.0

        # # # # #
        # M M A #
        # # # # #
        X, ymma, zmma, lam, xsi, eta, mu, zet, S, low, upp = mmasub(m, n_mma, OptIter, d, xmin, xmax, xold1, xold2, -10 * g, -10 * ∂g_∂d, g₁, ∂Ω∂d', low, upp, a0, am, C, d2)
        xold2  = xold1
        xold1  = d
        d      = X
        change = norm(d .- xold1)

        # # # # # # # # # #
        # Postprocessing  #
        # # # # # # # # # #
        global v_hist[OptIter] = g₁
        global g_hist[OptIter] = g

        #The residual vector of the KKT conditions is calculated:
        #residu,kktnorm,residumax = kktcheck(m,n,X,ymma,zmma,lam,xsi,eta,mu,zet,S, xmin,xmax,∂g_∂d,[0.0],zeros(size(d)),a0,a,C,d2);
        kktnorm = change
        println("Iter: ", OptIter, " Norm of change: ", kktnorm, " Objective: ", g)
        if mod(OptIter,1) == 0
            coord = getCoord(getX(dh0), dh0)
            postprocess_opt(Ψ, dh0, "results/Shape_with_contact" * string(OptIter))
            coord = getCoord(getX(dh), dh)
            postprocess_opt(a, dh, "results/DeformationC" * string(OptIter))
        end
        println("Objective: ", g_hist, " Constraint: ", g₁)
        if OptIter == 100
            break
        end

    end

    return g_hist, v_hist, OptIter, traction
end



Optimize(dh)

X_c = []
tract = []
for (key, val) ∈ traction
    append!(X_c, coord[key, 1])
    append!(tract, val)
end
ϵᵢⱼₖ = sortperm(X_c)
tract = tract[ϵᵢⱼₖ]
X_c = X_c[ϵᵢⱼₖ]
Plots.plot(X_c, tract, legend=false, marker=4, lc=:tomato, mc=:tomato)
OptIter = 2
Plots.plot(collect(1:OptIter), g_hist[1:OptIter],lc =:red, label="Objective",linewidth=3)
Plots.plot!(collect(1:OptIter), v_hist[1:OptIter], lc = :blue, label="Constraint",linewidth=3)

#using Makie
#OptIter=37
#fig = Figure()
#ax1, l1 = lines(fig[1, 1], 1..OptIter, g_hist[1:OptIter], color = :red)
#ax2, l2 = lines(fig[2, 1], 1..OptIter, v_hist[1:OptIter], color = :blue)
#Legend(fig[1:2, 2], [l1, l2], ["Objective", "Constraint"])
#fig
