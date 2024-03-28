using Mortar2D, ForwardDiff
using Ferrite, FerriteGmsh, FerriteMeshParser
using LinearSolve, SparseArrays # LinearSolvePardiso
using IterativeSolvers, IncompleteLU    # AlgebraicMultigrid
using SparseDiffTools
using Plots
using Printf
using JLD2
using Statistics # för var(λ)<


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

include("..//mma.jl")

#include("initLin.jl")

r₀ = 0.5
# Create two grids
y₁ = 0.999# 0.98
Δy = 1.0
Δx = 1.0
case = "box"
grid1 = createBoxMeshRev("box_1",  0.0, y₁, Δx, Δy, 1/10)
#grid2 = createBoxMeshRev("box_2",  -0.1, 0.0, 1.2, 1.0, 0.081)
grid2 = createBoxMeshRev("box_2",  0.0, 0.0, Δx, Δy, 1/11)
#_bothgrid1 = createBoxMeshRev("box_2", 0.0, 1.0, 1.0, 0.5, 0.08)

# # # # # # # # # #
# Finite element  #
# # # # # # # # # #
ip      = Lagrange{2,RefTetrahedron,1}()
qr      = QuadratureRule{2,RefTetrahedron}(3)
qr_face = QuadratureRule{1,RefTetrahedron}(2)
cv      = CellVectorValues(qr, ip)
fv      = FaceVectorValues(qr_face, ip)

# case  = "circle"
# grid1 = createCircleMesh("circle", 0.5, 1.5, r₀, 0.01)
# grid2 = createCircleMeshUp("circle2",0.5, 0.5001, r₀, 0.01) # inte rätt

# Merge into one grid
grid_tot = merge_grids(grid1, grid2; tol=1e-8)

grid1 = nothing
grid2 = nothing

# Create dofhandler with displacement field u
global dh = DofHandler(grid_tot)
add!(dh, :u, 2)
close!(dh)

# Extract CALFEM-style matrices
global coord, enod = getTopology(dh)
global register = index_nod_to_grid(dh, coord)

    if case == "box"
        # ------------------ #
        # Create master sets #
        # ------------------ #
        addfaceset!(dh.grid, "Γ_master", x -> x[2] ≈ 1.0)
        global Γs = getfaceset(dh.grid, "Γ_master")

        addnodeset!(dh.grid, "nₘ", x -> x[2] ≈ 1.0)
        global nₛ = getnodeset(dh.grid, "nₘ")

        # ----------------- #
        # Create slave sets #
        # ----------------- #
        addfaceset!(dh.grid, "Γ_slave", x -> ( x[2] ≈ y₁ || x[2] > 1.0 && ( x[1] ≈ 0.25 || x[1] ≈ 0.75 ) ) )
        global Γm = getfaceset(dh.grid, "Γ_slave")

        addnodeset!(dh.grid, "nₛ", x -> ( x[2] ≈ y₁ || x[2] > 1.0 && ( x[1] ≈ 0.25 || x[1] ≈ 0.75 ) ) )
        global nₘ = getnodeset(dh.grid, "nₛ")
    else
                # ------------------ #
        # Create master sets #
        # ------------------ #
        addfaceset!(dh.grid, "Γ_master",  x -> ((x[1]  - 0.5)^2 + (x[2] - 1.5  )^2) ≈ r₀^2)
        global Γs = getfaceset(dh.grid, "Γ_master")

        addnodeset!(dh.grid, "nₘ",  x -> ((x[1]  - 0.5)^2 + (x[2] - 1.5  )^2) ≈ r₀^2)
        global nₛ = getnodeset(dh.grid, "nₘ")

        # ----------------- #
        # Create slave sets #
        # ----------------- #
        #0.5, 0.5001
        addfaceset!(dh.grid, "Γ_slave", x -> ((x[1]  - 0.5)^2 + (x[2] - 0.5001  )^2) ≈ r₀^2)
        global Γm = getfaceset(dh.grid, "Γ_slave")

        addnodeset!(dh.grid, "nₛ", x -> ((x[1]  - 0.5)^2 + (x[2] - 0.5001 )^2) ≈ r₀^2)
        global nₘ = getnodeset(dh.grid, "nₛ")
    end
    #
    if case == "box"
        #=
        # ------------------ #
        # Create left | sets #
        # ------------------ #
        addfaceset!(dh.grid, "Γ_left", x -> x[2] < 0.99 && x[1] ≈ 0.0)
        global Γ_left = getfaceset(dh.grid, "Γ_left")

        addnodeset!(dh.grid, "nₗ", x -> x[2] < 0.99 && x[1] ≈ 0.0)
        global n_left = getnodeset(dh.grid, "nₗ")

        # ------------------ #
        # Create right  sets #
        # ------------------ #
        addfaceset!(dh.grid, "Γ_right", x -> x[2] < 0.99 && x[1] ≈ 1.0)
        global Γ_right = getfaceset(dh.grid, "Γ_right")

        addnodeset!(dh.grid, "nᵣ", x -> x[2] < 0.99 && x[1] ≈ 1.0)
        global n_right = getnodeset(dh.grid, "nᵣ")
        =#
        # % % % % % % % % % % % %
        # ------------------ #
        # Create left | sets #
        # ------------------ #
        addfaceset!(dh.grid, "Γ_left", x ->  x[1] ≈ 0.0)
        global Γ_left = getfaceset(dh.grid, "Γ_left")

        addnodeset!(dh.grid, "nₗ", x ->  x[1] ≈ 0.0)
        global n_left = getnodeset(dh.grid, "nₗ")

        # ------------------ #
        # Create right  sets #
        # ------------------ #
        addfaceset!(dh.grid, "Γ_right", x ->  x[1] ≈ 1.0)
        global Γ_right = getfaceset(dh.grid, "Γ_right")

        addnodeset!(dh.grid, "nᵣ", x ->  x[1] ≈ 1.0)
        global n_right = getnodeset(dh.grid, "nᵣ")
    end



# Extract all nbr nodes and dofs
global contact_dofs = getContactDofs(nₛ, nₘ)
global contact_nods = getContactNods(nₛ, nₘ)
global order = Dict{Int64,Int64}()
for (i, nod) ∈ enumerate(contact_nods)
    push!(order, nod => i)
end
global freec_dofs    = setdiff(1:dh.ndofs.x,contact_dofs)

# Define top nodeset for displacement controlled loading
if case == "box"
    addnodeset!(dh.grid, "Γ_top", x -> x[2] ≈ y₁ + Δy)
    global Γ_top = getnodeset(dh.grid, "Γ_top")

    addnodeset!(dh.grid, "n_top", x -> x[2] ≈ y₁ + Δy)
    global n_top = getnodeset(dh.grid, "n_top")

    # Define bottom nodeset subject to  u(X) = 0 ∀ X ∈ Γ_bot
    addnodeset!(dh.grid, "Γ_bot", x -> x[2] ≈ 0.0)
    global Γ_bot = getnodeset(dh.grid, "Γ_bot")

    addnodeset!(dh.grid, "n_bot", x -> x[2] ≈ 0.0)
    global n_bot = getnodeset(dh.grid, "n_bot")
    #
    global Γ_robin = union(
        getfaceset(dh.grid, "Γ_slave"),
        getfaceset(dh.grid, "Γ_left"),
        getfaceset(dh.grid, "Γ_right"),
        getfaceset(dh.grid, "Γ_master")
    )
    global n_robin = union(
        getnodeset(dh.grid, "nₛ"),
        getnodeset(dh.grid, "nₗ"),
        getnodeset(dh.grid, "nᵣ"),
        getnodeset(dh.grid, "nₘ")
    )
else
    addnodeset!(dh.grid, "Γ_top", x -> x[2] ≈ 1.5)
    global Γ_top = getnodeset(dh.grid, "Γ_top")

    addnodeset!(dh.grid, "n_top", x -> x[2] ≈ 1.5)
    global n_top = getnodeset(dh.grid, "n_top")

    # Define bottom nodeset subject to  u(X) = 0 ∀ X ∈ Γ_bot
    addnodeset!(dh.grid, "Γ_bot", x -> x[2] ≈ 0.5001)
    global Γ_bot = getnodeset(dh.grid, "Γ_bot")

    addnodeset!(dh.grid, "n_bot", x -> x[2] ≈ 0.5001)
    global n_bot = getnodeset(dh.grid, "n_bot")
    #
    #
    global Γ_robin = union(
        getfaceset(dh.grid, "Γ_slave"),
        #getfaceset(dh.grid, "Γ_left"),
        #getfaceset(dh.grid, "Γ_right"),
        getfaceset(dh.grid, "Γ_master")
    )
    global n_robin = union(
        getnodeset(dh.grid, "nₛ"),
        #getnodeset(dh.grid, "nₗ"),
        #getnodeset(dh.grid, "nᵣ"),
        getnodeset(dh.grid, "nₘ")
    )
end
# Final preparations for contact
global register = getNodeDofs(dh)
global X = getX(dh)
global coord = getCoordfromX(X)

# # # # # # # # #
# Init fictious #
# # # # # # # # #
global coord₀ = deepcopy(coord)



global free_d = []
for jnod in n_robin
    if in(jnod,n_left) || in(jnod,n_right)
        #append!(free_d, register[jnod, 1] )
    else
        #append!(free_d, register[jnod, 1] )
        append!(free_d, register[jnod, 2] )
    end
end
global locked_d = setdiff(1:dh.ndofs.x,free_d)

# Initialize tangents
global K    = create_sparsity_pattern(dh)
global Kψ   = create_sparsity_pattern(dh)
global a    = zeros(dh.ndofs.x)
global d    = zeros(dh.ndofs.x)
global Ψ    = zeros(dh.ndofs.x)
global Fᵢₙₜ = zeros(dh.ndofs.x)
global rc   = zeros(dh.ndofs.x)
global Fₑₓₜ = zeros(dh.ndofs.x)
global a    = zeros(dh.ndofs.x)
global Δa   = zeros(dh.ndofs.x)
global res  = zeros(dh.ndofs.x)

# boundary conditions for contact analysis
bcdof_top_o, _ = setBCXY_both(-0.01, dh, Γ_top)
bcdof_bot_o, _ = setBCXY_both(0.0, dh, Γ_bot)
#bcdof_top_o, _ = setBCXY(-0.01, dh, Γ_top)
#bcdof_bot_o, _ = setBCXY(0.0, dh, Γ_bot)
bcdof_o = [bcdof_top_o; bcdof_bot_o]
ϵᵢⱼₖ = sortperm(bcdof_o)
global bcdof_o = bcdof_o[ϵᵢⱼₖ]
global bcval_o = bcdof_o .* 0.0

bcdof_top_o2, _ = setBCXY_both(0.0, dh, Γ_top)
bcdof_bot_o2, _ = setBCXY_both(0.0, dh, Γ_bot)
#bcdof_top_o2, _ = setBCXY(0.0, dh, Γ_top)
#bcdof_bot_o2, _ = setBCXY(0.0, dh, Γ_bot)
bcdof_o2 = [bcdof_top_o2; bcdof_bot_o2]
ϵᵢⱼₖ = sortperm(bcdof_o)
global bcdof_o2 = bcdof_o2[ϵᵢⱼₖ]
global bcval_o2 = bcdof_o2 .* 0.0

# - For Linear solver..gmsh.model.add_physical_group(1, Lines[2:end-2], -1, "Γ_m")
global dr_dd      = similar(K)
global ∂rψ_∂d     = similar(K)
global ∂g_∂x      = zeros(size(a)) # behövs inte om vi har lokal funktion?
global ∂g_∂u      = zeros(size(d)) # behövs inte om vi har lokal funktion?
global ∂g₂_∂x     = zeros(size(a)) # behövs inte om vi har lokal funktion?
global ∂g₂_∂u     = zeros(size(d)) # behövs inte om vi har lokal funktion?
global λᵤ         = similar(a)
global λψ         = similar(a)
global Δ          = -0.05
global nloadsteps = 10
include("initOptLin.jl")
global asy_counter = zeros(dh.ndofs.x, 400)

function Optimize(dh)
    # Flytta allt nedan till init_opt?
        global dh0   = deepcopy(dh)
        global λψ    = similar(a)
        global λᵤ    = similar(a)
        global λᵥₒₗ  = similar(a)
        Vₘₐₓ         = 1.78 #1.1 * volume(dh, coord, enod)
       # global ε     = 1e6
       # global μ     = 1e3
        #l    = similar(a)
        #l   .= 0.5
        tol     = 1e-6
        OptIter = 0
        global true_iteration = 0
        global coord₀
        v_hist         = zeros(200)
        p_hist         = zeros(200)
        g_hist         = zeros(200)
        historia = zeros(200,4)
        global T = zeros(size(a))
        global T[bcdof_bot_o[bcdof_bot_o .% 2 .==0]] .= 1.0
        g₁ = 0.0
        g₂ = 0.0
    #
    while kktnorm > tol || OptIter < 200

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
            global ∂g_∂d
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
            global low #A() = A(Float64[],[]).
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

            global g_ini
            global pdofs    = bcdof_o
            global fdofs    = setdiff(1:length(a), pdofs)
            global locked_d = setdiff(1:length(a),free_d)
            global low
            global upp
            global traction
        # # # # # # # # # # # # # #
        OptIter += 1
        global true_iteration +=1

        if OptIter % 5 == 0 && g₂ < 0
            dh0 = deepcopy(dh)
            global d          = zeros(dh.ndofs.x)
            global xold1      = d[:]
            global xold2      = d[:]
            global low        = xmin
            global upp        = xmax
            OptIter           = 1
        end

       # if true_iteration % 5 != 0

            # # # # #
            # test  #
            # # # # #
            global nloadsteps = 10
            global μ = 0.0 # = 0.0 för bild # 1e3 # var μ = 1e4

            # # # # # # # # # # # # # #
            # Fictitious equillibrium #
            # # # # # # # # # # # # # #
            global coord₀ = getCoord(getX(dh0), dh0) # x₀
            Ψ, _, Kψ, _, λ = fictitious_solver_with_contact(d, dh0, coord₀, nloadsteps)

            # # # # # #
            # Filter  #
            # # # # # #
            global dh    = deepcopy(dh0)
            updateCoords!(dh, Ψ) # x₀ + Ψ = x
            global coord = getCoord(getX(dh), dh)
       # else
       #     λ = 1.0
       # end

        # # # # #
        # test  #
        # # # # #
        global nloadsteps = 10
        global ε = 1e5 # eller?

        # # # # # # # # #
        # Equillibrium  #
        # # # # # # # # #
        a, _, Fₑₓₜ, Fᵢₙₜ, K, traction = solver_C(dh, coord, Δ, nloadsteps)

        # # # # # # # # #
        # Sensitivities #
        # # # # # # # # #
        ∂rᵤ_∂x = similar(K)
        ∂rᵤ_∂x = drᵤ_dx_c(∂rᵤ_∂x, dh, mp, t, a, coord, enod, ε)
        dr_dd  = drψ(dr_dd, dh0, Ψ, λ, d, Γ_robin, coord₀)

        # # # # # # #
        # Objective #
        # # # # # # #
        # Max reaction force
        g     = -T' * Fᵢₙₜ
        ∂g_∂x = -T' * ∂rᵤ_∂x
        ∂g_∂u = -T' * K
        # Compliance
        # g            = -a[pdofs]' * Fᵢₙₜ[pdofs]
        # ∂g_∂x[fdofs] = -a[pdofs]' * ∂rᵤ_∂x[pdofs, fdofs]
        # ∂g_∂u[fdofs] = -a[pdofs]' * K[pdofs, fdofs]

        # # # # # # #
        # Adjoints  #
        # # # # # # #


        solveq!(λᵤ, K', ∂g_∂u, bcdof_o2, bcval_o2)
        solveq!(λψ, Kψ', ∂g_∂x' - ∂rᵤ_∂x' * λᵤ, bcdof_o2, bcval_o2)

        # # # # # # # # # # #
        # Full sensitivity  #
        # # # # # # # # # # #
        ∂g_∂d = (-transpose(λψ) * dr_dd)'

        # # # # # # # # # # #
        # Volume constraint #
        # # # # # # # # # # #
        g₁    = volume(dh, coord, enod) / Vₘₐₓ - 1.0
        ∂Ω_∂x = volume_sens(dh, coord)
        solveq!(λᵥₒₗ, Kψ, ∂Ω_∂x, bcdof_o2, bcval_o2)
        ∂Ω∂d = Real.(-transpose(λᵥₒₗ) * dr_dd ./ Vₘₐₓ)
        #∂g_∂d[locked_d] .= 0.0 # fulfix?

        # # # # # # # # # # #
        # Volume constraint #
        # # # # # # # # # # #
        g₁    = volume(dh,coord,enod) / Vₘₐₓ - 1.0
        ∂Ω_∂x = volume_sens(dh,coord)
        solveq!(λᵥₒₗ, Kψ, ∂Ω_∂x, bcdof_o2, bcval_o2.*0);
        ∂Ω∂d  = Real.( -transpose(λᵥₒₗ)*dr_dd ./ Vₘₐₓ) ;
        #∂Ω∂d[locked_d] .= 0.0

        # # # # # # # # # # # #
        # Pressure constraint #
        # # # # # # # # # # # #
        # p = 2
        # X_ordered = getXfromCoord(coord)
        # g₂         = contact_pnorm_s(X_ordered, a, ε, p) / 0.5 - 1.0
        # ∂g₂_∂x     = ForwardDiff.gradient(x -> contact_pnorm_ordered_s(x, a, ε, p), getXinDofOrder(dh, X_ordered, coord))
        # ∂g₂_∂u     = ForwardDiff.gradient(u -> contact_pnorm_s(X_ordered, u, ε, p), a)

        # solveq!(λᵤ, K',  ∂g₂_∂u, bcdof_o, bcval_o)
        # solveq!(λψ, Kψ', ∂g₂_∂x - ∂rᵤ_∂x' * λᵤ, bcdof_o2, bcval_o2)
        # ∂g₂_∂d            = Real.( (-transpose(λψ) * dr_dd)' ./ 0.5 )'
        postprocess_opt(∂g_∂d, dh, "results/🛸-bild" * string(true_iteration))
        # # # # #
        # M M A #
        # # # # #
        d_new, ymma, zmma, lam, xsi, eta, mu, zet, S, low, upp = mmasub(m, n_mma, OptIter, d[:], xmin[:], xmax[:], xold1[:], xold2[:], g .* 100, ∂g_∂d .* 100, g₁ .* 100, ∂Ω∂d .* 100, low, upp, a0, am, C, d2)
        #d_new, ymma, zmma, lam, xsi, eta, mu, zet, S, low, upp = mmasub(m, n_mma, OptIter, d, xmin, xmax, xold1, xold2, -10 * g, -10 * ∂g_∂d, hcat([g₁; g₂]), vcat([∂Ω∂d; ∂g₂_∂d]), low, upp, a0, am, C, d2)
        xold2  = xold1
        xold1  = d
        d      = d_new
        change = norm(d .- xold1)

        # # # # # # # # # #
        # Postprocessing  #
        # # # # # # # # # #
        v_hist[true_iteration] = g₁
        p_hist[true_iteration] = g₂
        g_hist[true_iteration] = g

        #historia[true_iteration,:] = [∂g_∂d[677] ∂g_∂d[678] coord[273,1] coord[273,2]]

        #The residual vector of the KKT conditions is calculated:
        #residu,kktnorm,residumax = kktcheck(m,n,X,ymma,zmma,lam,xsi,eta,mu,zet,S, xmin,xmax,∂g_∂d,[0.0],zeros(size(d)),a0,a,C,d2);
        kktnorm = change
        println("Iter: ", true_iteration, " Norm of change: ", kktnorm, " Objective: ", g)
        if mod(OptIter,1) == 0
            coord = getCoord(getX(dh0), dh0)
            postprocess_opt(Ψ, dh0, "results/Current design" * string(true_iteration))
            postprocess_opt(d, dh0, "results/design_variables" * string(true_iteration))
        end

        println("Objective: ", g_hist[1:true_iteration], " Constraint: ", v_hist[1:true_iteration] , p_hist[1:true_iteration])


        p2 = plot(1:true_iteration,[p_hist[1:true_iteration],g_hist[1:true_iteration]],label = ["Constraint" "Objective"])
        display(p2)

        #p3 = plot(1:true_iteration,g_hist[1:true_iteration],legend=false, marker=3, reuse = false, lc =:darkgreen)
        #display(p3)

        #if true_iteration == 1
        #    g_ini = 0
        #    n     = 0
        #    xval  = 0
        #elseif true_iteration == 2
        #    @save "tva.jld2"
        #elseif true_iteration == 200
        #    @save "steg100.jld2"
        #    break
        #end
    end
    return g_hist, v_hist, OptIter, traction, historia
end

function plotTraction()
    traction = ExtractContactTraction(a, ε, coord)
    X_c = []
    tract = []
    for (key, val) ∈ traction
        append!(X_c, coord[key, 1])
        append!(tract, val)
    end
    ϵᵢⱼₖ = sortperm(X_c)
    tract = tract[ϵᵢⱼₖ]
    X_c = X_c[ϵᵢⱼₖ]
    return X_c, tract
end

# plot(coord[collect(n_robin),1], ∂g_∂d[free_d], seriestype=:scatter)
g_hist, v_hist, OptIter, traction, historia = Optimize(dh)
