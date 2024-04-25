using ForwardDiff #, Mortar2D
using Ferrite, FerriteGmsh#, FerriteMeshParser
using LinearSolve, SparseArrays # LinearSolvePardiso
using IterativeSolvers, IncompleteLU    # AlgebraicMultigrid
using Plots, Printf, JLD2, Statistics
#using LazySets: convex_hull
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
# FEM quantities
ip      = Lagrange{2,RefTetrahedron,1}()
qr      = QuadratureRule{2,RefTetrahedron}(1)
qr_face = QuadratureRule{1,RefTetrahedron}(1)
cv      = CellVectorValues(qr, ip)
fv      = FaceVectorValues(qr_face, ip)
# Create two grids
r₀      = 0.5
h       = 0.15 #0.2 #0.15 #0.075 # 0.075 #0.15 #0.1
Δx      = r₀  # * π/2  # 0.5
x₀      = 0.0
y₀      = 0.5
Δy      = 0.501 #1.001
grid1   = createHalfCircleMesh("circle", 0.0, 1.5, r₀, 0.75h)
grid2   = createBoxMeshRev("box_1",  0.0, y₀, Δx, Δy, 0.5h)
case    = "box"
        # - - - - - Eller?
        #y₀      = 0.501
        #grid2   = createHalfCircleMeshFlipped("circlef", x₀, y₀, r₀, h*0.4)
        #case    = "circle"
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

if case == "box"

    # ------------------ #
    # Create master sets #
    # ------------------ #
    addfaceset!(dh.grid, "Γ_slave", x -> x[2] ≈ y₀ + Δy)
    global Γs = getfaceset(dh.grid, "Γ_slave")
    addnodeset!(dh.grid, "nₛ", x -> x[2] ≈ y₀ + Δy)
    global nₛ = getnodeset(dh.grid, "nₛ")
    ###
    #addfaceset!(dh.grid, "Γ_master", x -> x[2] ≈ y₀ + Δy)
    #global Γm = getfaceset(dh.grid, "Γ_master")
    #addnodeset!(dh.grid, "nₘ", x -> x[2] ≈ y₀ + Δy)
    #global nₘ = getnodeset(dh.grid, "nₘ")


    # ------------------ #
    # Create left | sets #
    # ------------------ #
    addfaceset!(dh.grid, "Γ_left", x ->  x[1] ≈ 0.0)
    global Γ_left = getfaceset(dh.grid, "Γ_left")

    addnodeset!(dh.grid, "nₗ", x ->  x[1] ≈ 0.0)
    global n_left = getnodeset(dh.grid, "nₗ")
    # ------------------- #
    # Create right | sets #
    # ------------------- #
    addfaceset!(dh.grid, "Γ_right", x ->  x[1] ≈ Δx)
    global Γ_right = getfaceset(dh.grid, "Γ_right")

    addnodeset!(dh.grid, "nᵣ", x ->  x[1] ≈ Δx)
    global nᵣ = getnodeset(dh.grid, "nᵣ")

else
    # ------------------ #
    # Create slave sets  #
    # ------------------ #
    addfaceset!(dh.grid, "Γ_slave", x -> ((x[1] - x₀)^2 + (x[2] - y₀ )^2) ≈ r₀^2 )
    global Γs = getfaceset(dh.grid, "Γ_slave")

    addnodeset!(dh.grid, "nₛ", x -> ((x[1] - x₀)^2 + (x[2] - y₀ )^2) ≈ r₀^2 )
    global nₛ = getnodeset(dh.grid, "nₛ")
    # ------------------ #
    # Create left | sets #
    # ------------------ #
    addfaceset!(dh.grid, "Γ_left", x ->  x[1] ≈ 0.0)
    global Γ_left = getfaceset(dh.grid, "Γ_left")

    addnodeset!(dh.grid, "nₗ", x ->  x[1] ≈ 0.0)
    global n_left = getnodeset(dh.grid, "nₗ")
end

# ----------------- #
# Create slave sets #
# ----------------- #
 addfaceset!(dh.grid, "Γ_master", x -> ((x[1] - 0.0 )^2 + (x[2] - 1.5)^2) ≈ r₀^2 )
 global Γm = getfaceset(dh.grid, "Γ_master")
 addnodeset!(dh.grid, "nₘ", x -> ((x[1] - 0.0)^2 + (x[2] - 1.5)^2) ≈ r₀^2 )
 global nₘ = getnodeset(dh.grid, "nₘ")

#addfaceset!(dh.grid, "Γ_slave", x -> ((x[1] - 0.0 )^2 + (x[2] - 1.5)^2) ≈ r₀^2 )
#global Γs = getfaceset(dh.grid, "Γ_slave")
#
#addnodeset!(dh.grid, "nₛ", x -> ((x[1] - 0.0)^2 + (x[2] - 1.5)^2) ≈ r₀^2 )
#global nₛ = getnodeset(dh.grid, "nₛ")


# Extract all nbr nodes and dofs
global contact_dofs = getContactDofs(nₛ, nₘ)
global contact_nods = getContactNods(nₛ, nₘ)
global order        = Dict{Int64,Int64}()
for (i, nod) ∈ enumerate(contact_nods)
    push!(order, nod => i)
end
global freec_dofs    = setdiff(1:dh.ndofs.x,contact_dofs)

# Define top nodeset for displacement controlled loading
addfaceset!(dh.grid, "Γ_top", x -> x[2] ≈ 1.5)
global Γ_top = getfaceset(dh.grid, "Γ_top")

addnodeset!(dh.grid, "n_top", x -> x[2] ≈ 1.5)
global n_top = getnodeset(dh.grid, "n_top")

if case == "box"
    # Define bottom nodeset subject to  u(X) = 0 ∀ X ∈ Γ_bot
    addfaceset!(dh.grid, "Γ_bot", x -> x[2] ≈ y₀)
    global Γ_bot = getfaceset(dh.grid, "Γ_bot")

    addnodeset!(dh.grid, "n_bot", x -> x[2] ≈ y₀)
    global n_bot = getnodeset(dh.grid, "n_bot")
else
    # Define bottom nodeset subject to  u(X) = 0 ∀ X ∈ Γ_bot
    addfaceset!(dh.grid, "Γ_bot", x -> x[2] ≈ y₀)
    global Γ_bot = getfaceset(dh.grid, "Γ_bot")

    addnodeset!(dh.grid, "n_bot", x -> x[2] ≈ y₀)
    global n_bot = getnodeset(dh.grid, "n_bot")
end

# Final preparations for contact
global register = getNodeDofs(dh)
global X = getX(dh)
global coord = getCoordfromX(X)

# # # # # # # # #
# Init fictious #
# # # # # # # # #
global coord₀ = deepcopy(coord)

global Γ_robin = union(
   #getfaceset(dh.grid, "Γ_top"), # test: bc-ränder
    #getfaceset(dh.grid, "Γ_bot"), # test: bc-ränder
    getfaceset(dh.grid, "Γ_slave"),
    getfaceset(dh.grid, "Γ_right"),
    getfaceset(dh.grid, "Γ_master")
)
global n_robin = union(
    #getnodeset(dh.grid, "n_top"), # test: bc-ränder
    #getnodeset(dh.grid, "n_bot"), # test: bc-ränder
    getnodeset(dh.grid, "nₛ"),
    getnodeset(dh.grid, "nᵣ"), # Intressanta resultat när endast n_r exkluderades
    getnodeset(dh.grid, "nₘ")
)


global free_d = []
for jnod in n_robin
    if in(jnod,n_bot) || in(jnod,n_top)
        #append!(free_d, register[jnod, 1] )
    else
        append!(free_d, register[jnod, 1] )
        append!(free_d, register[jnod, 2] )
    end
end
global locked_d = setdiff(1:dh.ndofs.x,free_d)

# ------------------ #
# To test asymptotes #
# ------------------ #
global asy_counter = zeros(dh.ndofs.x,400)

function Optimize(dh)
        #
        #
        #
        # Initialize tangents
        global K    = create_sparsity_pattern(dh)
        global Kψ   = create_sparsity_pattern(dh)
        global a    = zeros(dh.ndofs.x)
        global d    = zeros(dh.ndofs.x)
        global Ψ    = zeros(dh.ndofs.x)
        global Fᵢₙₜ= zeros(dh.ndofs.x)
        global Fₑₓₜ= zeros(dh.ndofs.x)
        # boundary conditions for contact analysis
        bcdof_top_o, _    = setBCY(0.0, dh, n_top)
        bcdof_bot_o, _    = setBCY(0.0, dh, n_bot)
        bcdof_left_o, _   = setBCX(0.0, dh, n_left)
        bcdof_o           = [bcdof_top_o; bcdof_bot_o; bcdof_left_o]
        ϵᵢⱼₖ             = sortperm(bcdof_o)
        global bcdof_o    = bcdof_o[ϵᵢⱼₖ]
        global bcval_o    = bcdof_o .* 0.0
        bcdof_top_o2, _   = setBCY(0.0, dh, n_top)
        bcdof_bot_o2, _   = setBCY(0.0, dh, n_bot)
        bcdof_left_o2, _  = setBCX(0.0, dh, n_left)
        bcdof_o2          = [bcdof_top_o2; bcdof_bot_o2; bcdof_left_o2]
        ϵᵢⱼₖ             = sortperm(bcdof_o)
        global bcdof_o2   = bcdof_o2[ϵᵢⱼₖ]
        global bcval_o2   = bcdof_o2 .* 0.0
        global dr_dd      = similar(K)
        global ∂g_∂x      = zeros(size(a)) # behövs inte om vi har lokal funktion?
        global ∂g_∂u      = zeros(size(d)) # behövs inte om vi har lokal funktion?
        global λᵤ         = similar(a)
        global λψ         = similar(a)
        global Δ          = -0.05
        global nloadsteps = 10
        global kktnorm    = 1.0

        global dr_dd = similar(K)
        # Material parameters
        global mp₀   = [1.0 5.0] # [0.1 5.0]#
        global mp    = [175 80.769230769230759]
        global t     = 1.0
        # Optimization parameters
        global m             = 1;
        global n_mma         = length(d);
        global epsimin       = 0.0000001;
        global xvalue        = d[:];
        global xold1         = xvalue;
        global xold2         = xvalue;
        global xmin          = zeros(n_mma)#-ones(n_mma)/1000;
        global xmax          = zeros(n_mma)# ones(n_mma)/1000;
        global C             = 1000*ones(m);
        global d2            = zeros(m);
        global a0            = 1;
        global am            = zeros(m);
        global outeriter     = 0;
        global kkttol        = 0.001;
        global changetol     = 0.001;
        global kktnorm       = kkttol + 10;
        global outit         = 0;
        global change        = 1;
        global xmax         .=  0.2 #  0.2 # 0.2 #* π/2
        global xmin         .= -0.2 # -0.2 # 0.2 #* π/2
        global low           =-ones(n_mma);
        global upp           = ones(n_mma);
        # Flytta allt nedan till init_opt?
        global dh0     = deepcopy(dh)
        global λψ      = similar(a)
        global λᵤ      = similar(a)
        global λᵥₒₗ   = similar(a)
        Vₘₐₓ          = 0.75 #0.5# 0.6 #0.75 # volume(dh, coord, enod) # 4(r₀ + Δy) * Δx # 0.5 #0.9 #1.1 * volume(dh, coord, enod)
        tol            = 1e-6
        OptIter        = 0
        global true_iteration = 0
        global coord₀
        v_hist         = zeros(1000)
        p_hist         = zeros(1000)
        g_hist         = zeros(1000)
        historia       = zeros(1000,4)
        global T       = zeros(size(a))
        global T[bcdof_bot_o[bcdof_bot_o .% 2 .==0]] .= -1.0
        global T[bcdof_top_o[bcdof_top_o .% 2 .==0]] .=  1.0
        g₁ = 0.0
        g₂ = 0.0
        # # # # # # # # # # # # # #
        #       Definitions       #
        # # # # # # # # # # # # # #
    #
    while kktnorm > tol || OptIter < 200
        # # # # # # # # # # # # # #
        OptIter += 1
        true_iteration +=1

        if true_iteration % 10 == 0 && true_iteration < 200
            dh0 = deepcopy(dh)
            d          = zeros(dh.ndofs.x)
            xold1      = d[:]
            xold2      = d[:]
            #if true_iteration > 50
            #    xmin .= -0.1
            #    xmax .=  0.1
            #end
            low        = xmin
            upp        = xmax
            OptIter    = 1
        end

            # # # # #
            # test  #
            # # # # #
            global nloadsteps = 20
            # 1e5 för h=0.015
            # 5e3 för h=0.03
            # 1e4 standard
            global μ = 1e4 #5e2 # 1e2 # 1e3/2 #1e4

            # # # # # # # # # # # # # #
            # Fictitious equillibrium #
            # # # # # # # # # # # # # #
            global coord₀ = getCoord(getX(dh0), dh0) # x₀
            Ψ, _, Kψ, _, λ = fictitious_solver_with_contact_half(d, dh0, coord₀, nloadsteps)

            # # # # # #
            # Filter  #
            # # # # # #
            global dh    = deepcopy(dh0)
            updateCoords!(dh, Ψ) # x₀ + Ψ = x
            global coord = getCoord(getX(dh), dh)

        # # # # #
        # test  #
        # # # # #
        global nloadsteps = 10
        global ε = 1e4 #1e4 # 1e5 #1e4 # 2?

        # # # # # # # # #
        # Equillibrium  #
        # # # # # # # # #
        a, _, Fₑₓₜ, Fᵢₙₜ, K, _ = solver_C_half(dh, coord, Δ, nloadsteps)

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
        g     = - T' * Fᵢₙₜ
        ∂g_∂x =  -T' * ∂rᵤ_∂x #
        ∂g_∂u =  -T' * K # ?

        # # # # # # #
        # Adjoints  #
        # # # # # # #
        solveq!(λᵤ, K',  ∂g_∂u, bcdof_o, bcval_o)
        solveq!(λψ, Kψ', ∂g_∂x' - ∂rᵤ_∂x' * λᵤ, bcdof_o2, bcval_o2)

        # # # # # # # # # # #
        # Full sensitivity  #
        # # # # # # # # # # #
        ∂g_∂d            = (-transpose(λψ) * dr_dd)'
        #∂g_∂d[locked_d] .= 0.0 # fulfix?

        # # # # # # # # # # #
        # Volume constraint #
        # # # # # # # # # # #
        g₁    = volume(dh,coord,enod) / Vₘₐₓ - 1.0
        ∂Ω_∂x = volume_sens(dh,coord)./ Vₘₐₓ
        solveq!(λᵥₒₗ, Kψ, ∂Ω_∂x, bcdof_o2, bcval_o2.*0);
        ∂Ω∂d  = Real.( -transpose(λᵥₒₗ)*dr_dd ) ;
        #∂Ω∂d[locked_d] .= 0.0

        # # # # # # # # # # # #
        # Pressure constraint #
        # # # # # # # # # # # #
        #=
        p = 2
        X_ordered = getXfromCoord(coord)
        g₂         = contact_pnorm_s(X_ordered, a, ε, p) / 0.5 - 1.0
        ∂g₂_∂x     = ForwardDiff.gradient(x -> contact_pnorm_ordered_s(x, a, ε, p), getXinDofOrder(dh, X_ordered, coord))
        ∂g₂_∂u     = ForwardDiff.gradient(u -> contact_pnorm_s(X_ordered, u, ε, p), a)

        solveq!(λᵤ, K',  ∂g₂_∂u, bcdof_o, bcval_o)
        solveq!(λψ, Kψ', ∂g₂_∂x - ∂rᵤ_∂x' * λᵤ, bcdof_o2, bcval_o2)
        ∂g₂_∂d            = Real.( (-transpose(λψ) * dr_dd)' ./ 0.5 )'
        =#

        # # # # #
        # M M A #
        # # # # #
        d_old   = d
        low_old = low
        upp_old = upp
        #d_new, ymma, zmma, lam, xsi, eta, mu, zet, S, low, upp = mmasub(m, n_mma, OptIter, d[:], xmin[:], xmax[:], xold1[:], xold2[:], g, ∂g_∂d, g₁.*100, ∂Ω∂d.*100, low, upp, a0, am, C, d2)
        d_new, ymma, zmma, lam, xsi, eta, mu, zet, S, low, upp = mmasub(m, n_mma, OptIter, d[:], xmin[:], xmax[:], xold1[:], xold2[:], g, ∂g_∂d, g₁.*10, ∂Ω∂d.*10, low, upp, a0, am, C, d2)
        α      = 1.0  # 0.4 # 0.1 #
        if true_iteration > 150
            α = 0.5 # kanske testa 0.5 / 0.6 / 0.75 ... kan också testa annat villkor för dämpning
        end
        d_new  = d_old   + α .* (d_new - d_old)
        low    = low_old + α .* (low   - low_old)
        upp    = upp_old + α .* (upp   - upp_old)
        xold2  = xold1
        xold1  = d
        d      = d_new
        change = norm(d .- xold1)
        kktnorm = change

        # # # # # # # # # #
        # Postprocessing  #
        # # # # # # # # # #
        v_hist[true_iteration] = g₁
        p_hist[true_iteration] = g₂
        g_hist[true_iteration] = g

        # Print results
        println("Iter: ", true_iteration, " Norm of change: ", kktnorm, " Objective: ", g)
        postprocess_opt(Ψ, dh0, "results/Current design" * string(true_iteration))
        postprocess_opt(d, dh0, "results/design_variables" * string(true_iteration))
        postprocess_opt(∂g_∂d,dh,"results/🛸" * string(true_iteration))
        println("Objective: ", g_hist[1:true_iteration], " Constraint: ", v_hist[1:true_iteration] , p_hist[1:true_iteration])
        p2 = plot(1:true_iteration,[v_hist[1:true_iteration].*10,g_hist[1:true_iteration]],label = ["Volume Constraint" "Objective"], marker = :cross)
        display(p2)
        GC.gc()
        @save "färdig_cyl.jld2" a dh dh0 OptIter v_hist g_hist d ε μ true_iteration
        #Profile.take_heap_snapshot( "snapshot.heapsnapshot" * string(true_iteration) )
    end
    return g_hist, v_hist, OptIter, historia
end
g_hist, v_hist, OptIter, historia = Optimize(dh)
@save "färdig_cyl.jld2" a dh dh0 OptIter v_hist g_hist d ε μ true_iteration

# Byta från att ladda Mortar2D från packages till lokal mapp


# Plot traction , can be moved to function...
traction = ExtractContactTraction(a, 1e5, coord)
#traction = τ_c
X_c = []
tract = []
for (key, val) ∈ traction
    append!(X_c, coord[key, 1])
    append!(tract, val)
end
ϵᵢⱼₖ = sortperm(X_c)
tract = tract[ϵᵢⱼₖ]
X_c = X_c[ϵᵢⱼₖ]
p = plot(X_c, tract, legend=false, marker=4, lc=:tomato, mc=:tomato)
display(p)


vtk_grid("gap📄" , dh) do vtkfile
    vtk_point_data(vtkfile, dh, a) # displacement field
end
