using Mortar2D, ForwardDiff
using Ferrite, FerriteGmsh#, FerriteMeshParser
using LinearSolve, SparseArrays # LinearSolvePardiso
using IterativeSolvers, IncompleteLU    # AlgebraicMultigrid
using Plots, Printf, JLD2, Statistics
using LazySets: convex_hull
include("..//mesh_reader.jl")
include("Contact//contact_help.jl")
include("assemLin.jl")
include("assemElemLin.jl")
include("..//material.jl")
include("..//fem.jl")
include("run_linear.jl")
include("sensitivitiesLin.jl")
include("..//mma.jl")
#include("objects.jl")
# FEM quantities
ip = Lagrange{2,RefTetrahedron,1}()
qr = QuadratureRule{2,RefTetrahedron}(1)
qr_face = QuadratureRule{1,RefTetrahedron}(1)
cv = CellVectorValues(qr, ip)
fv = FaceVectorValues(qr_face, ip)
#fem = FEM(
#    create_sparsity_pattern(dh),
#    create_sparsity_pattern(dh),
#    zeros(dh.ndofs.x),
#    zeros(dh.ndofs.x),
#    zeros(dh.ndofs.x),
#    zeros(dh.ndofs.x),
#    zeros(dh.ndofs.x),
#    zeros(dh.ndofs.x),
#    zeros(dh.ndofs.x),
#    zeros(dh.ndofs.x)
#)
#
r₀ = 0.5
# Create two grids
xₗ = 0.0
yₗ = 0.5
Δx = 1.0
Δy = 0.5
h  = 0.025
case = "box"
rounded = false
# najs för ~0.05
if rounded == true
    filename1 = "box_rounded"
    grid1 = createBoxMeshRounded_Flipped(filename1, 0.35,  2yₗ, Δy, h)
    Γ_1   = getBoundarySet(grid1);
else
    filename1 = "circle"
    grid1 = createCircleMesh(filename1,  Δx/2, yₗ + 2Δy, r₀, h)
end
filename2 = "box"
grid2 = createBoxMeshRev(filename2, xₗ, yₗ, Δx, 0.501, h)
## Merge into one grid
grid_tot = merge_grids(grid1, grid2; tol=1e-6)
grid1 = nothing
grid2 = nothing
# Create dofhandler with displacement field u
global dh = DofHandler(grid_tot)
add!(dh, :u, 2)
close!(dh)
# Exrtact full boundary
Γ_all   = Ferrite.__collect_boundary_faces(dh.grid);
addfaceset!(dh.grid,"Γ_all", Γ_all);
Γ_all  = getfaceset(dh.grid, "Γ_all");
#
n_all = getBoundarySet(dh.grid, Γ_all);
addnodeset!(dh.grid, "n_all", n_all);
# Extract CALFEM-style matrices
global coord, enod = getTopology(dh)
global register = index_nod_to_grid(dh, coord)
if case == "box" && rounded == false
    # ------------------ #
    # Create slave sets #
    # ------------------ #
    addfaceset!(dh.grid, "Γ_slave", x -> x[2] ≈ 1.001)
    global Γs = getfaceset(dh.grid, "Γ_slave")
    addnodeset!(dh.grid, "nₛ", x -> x[2] ≈ 1.001)
    global nₛ = getnodeset(dh.grid, "nₛ")
    # ------------------ #
    # Create left | sets #
    # ------------------ #
    addfaceset!(dh.grid, "Γ_left", x -> x[2] < 1.001 && x[1] ≈ xₗ)
    global Γ_left = getfaceset(dh.grid, "Γ_left")
    addnodeset!(dh.grid, "nₗ", x -> x[2] < 1.001 && x[1] ≈ xₗ)
    global n_left = getnodeset(dh.grid, "nₗ")
    # ------------------ #
    # Create right  sets #
    # ------------------ #
    addfaceset!(dh.grid, "Γ_right", x -> x[2] < 1.001 && x[1] ≈ xₗ + Δx)
    global Γ_right = getfaceset(dh.grid, "Γ_right")
    addnodeset!(dh.grid, "nᵣ", x -> x[2] < 1.001 && x[1] ≈ xₗ + Δx)
    global n_right = getnodeset(dh.grid, "nᵣ")
elseif rounded == true
    # ------------------ #
    # Create slave sets #
    # ------------------ #
    addfaceset!(dh.grid, "Γ_slave", x -> x[2] ≈ 1.001)
    global Γs = getfaceset(dh.grid, "Γ_slave")
    addnodeset!(dh.grid, "nₛ", x -> x[2] ≈ 1.001)
    global nₛ = getnodeset(dh.grid, "nₛ")
    # ------------------ #
    # Create left | sets #
    # ------------------ #
    addfaceset!(dh.grid, "Γ_left", x ->  x[1] ≈ xₗ)
    global Γ_left = getfaceset(dh.grid, "Γ_left")
    addnodeset!(dh.grid, "nₗ", x ->  x[1] ≈ xₗ)
    global n_left = getnodeset(dh.grid, "nₗ")
    # ------------------ #
    # Create right  sets #
    # ------------------ #
    addfaceset!(dh.grid, "Γ_right", x ->  x[1] ≈ xₗ + Δx)
    global Γ_right = getfaceset(dh.grid, "Γ_right")
    addnodeset!(dh.grid, "nᵣ", x ->  x[1] ≈ xₗ + Δx)
    global n_right = getnodeset(dh.grid, "nᵣ")
else
    # ------------------ #
    # Create master sets #
    # ------------------ #
    addfaceset!(dh.grid, "Γ_slave", x -> ((x[1] - r₀)^2 + (x[2] - 0.5001 )^2) ≈ r₀^2 )
    global Γs = getfaceset(dh.grid, "Γ_slave")

    addnodeset!(dh.grid, "nₛ", x -> ((x[1] - r₀)^2 + (x[2] - 0.5001 )^2) ≈ r₀^2 )
    global nₛ = getnodeset(dh.grid, "nₛ")
end
# ----------------- #
# Create slave sets #
# ----------------- #
if rounded == true
    addfaceset!(dh.grid,"Γ_master", x -> x ∈ Γ_1 );
    Γm = getfaceset(dh.grid, "Γ_master");
    #Γm = intersect(Γm, Γ_all);
    #
    nₘ = getBoundarySet(dh.grid,Γm)
    addnodeset!(dh.grid, "nₘ" ,nₘ)
else
    addfaceset!(dh.grid, "Γ_master", x -> ((x[1] - r₀)^2 + (x[2] - 1.5)^2) ≈ r₀^2 )
    global Γm = getfaceset(dh.grid, "Γ_master")
    addnodeset!(dh.grid, "nₘ", x -> ((x[1] - r₀)^2 + (x[2] - 1.5)^2) ≈ r₀^2 )
    global nₘ = getnodeset(dh.grid, "nₘ")
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
addfaceset!(dh.grid, "Γ_top", x -> x[2] ≈ 1.5)
global Γ_top = getfaceset(dh.grid, "Γ_top")
addnodeset!(dh.grid, "n_top", x -> x[2] ≈ 1.5)
global n_top = getnodeset(dh.grid, "n_top")
#
if case == "box"
    # Define bottom nodeset subject to  u(X) = 0 ∀ X ∈ Γ_bot
    addfaceset!(dh.grid, "Γ_bot", x -> x[2] ≈ yₗ)
    global Γ_bot = getfaceset(dh.grid, "Γ_bot")

    addnodeset!(dh.grid, "n_bot", x -> x[2] ≈ yₗ)
    global n_bot = getnodeset(dh.grid, "n_bot")
else
    # Define bottom nodeset subject to  u(X) = 0 ∀ X ∈ Γ_bot
    addfaceset!(dh.grid, "Γ_bot", x -> x[2] ≈ 0.5001)
    global Γ_bot = getfaceset(dh.grid, "Γ_bot")

    addnodeset!(dh.grid, "n_bot", x -> x[2] ≈ 0.5001)
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
global free_d = []
for jnod in n_robin
    if in(jnod,n_left) || in(jnod,n_right)
        append!(free_d, register[jnod, 1] )
        append!(free_d, register[jnod, 2] ) ## Fundera på om detta skall vara med
    else
        append!(free_d, register[jnod, 1] )
        append!(free_d, register[jnod, 2] )
    end
end
global locked_d = setdiff(1:dh.ndofs.x,free_d)
# Initialize tangents <- fem
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
# "fem"
#  boundary conditions for contact analysis
bcdof_top_o, _ = setBCXY(0.0, dh, n_top)
bcdof_bot_o, _ = setBCXY(0.0, dh, n_bot)
bcdof_o = [bcdof_top_o; bcdof_bot_o]
ϵᵢⱼₖ = sortperm(bcdof_o)
global bcdof_o = bcdof_o[ϵᵢⱼₖ]
global bcval_o = bcdof_o .* 0.0
# fictitious bcs
bcdof_top_o2, _ = setBCXY(0.0, dh, n_top)
bcdof_bot_o2, _ = setBCXY(0.0, dh, n_bot)
bcdof_o2 = [bcdof_top_o2; bcdof_bot_o2]
ϵᵢⱼₖ = sortperm(bcdof_o)
global bcdof_o2 = bcdof_o2[ϵᵢⱼₖ]
global bcval_o2 = bcdof_o2 .* 0.0
# sensitivities
global dr_dd      = similar(K)
global ∂rψ_∂d     = similar(K)
global ∂g_∂x      = zeros(size(a)) # behövs inte om vi har lokal funktion?
global ∂g_∂u      = zeros(size(d)) # behövs inte om vi har lokal funktion?
global ∂g₂_∂x     = zeros(size(a)) # behövs inte om vi har lokal funktion?
global ∂g₂_∂u     = zeros(size(d)) # behövs inte om vi har lokal funktion?
global λᵤ         = similar(a) # Intermediate, kanske bara behöver en adjoint?
global λψ         = similar(a)
global λᵥₒₗ       = similar(a)
# "fem"
global Δ          = -0.01
global nloadsteps = 10
# mma
include("initOptLin.jl")
#
function Optimize(dh)
    # Flytta allt nedan till init_opt?
        # - - - - - - -  #
        # Initialization #
        # - - - - - - -  #
        global dh0     = deepcopy(dh)
        global λψ      = similar(a)
        global λᵤ      = similar(a)
        global λᵥₒₗ   = similar(a)
        Vₘₐₓ          = 1.0#0.9
        tol            = 1e-3
        OptIter        = 0
        true_iteration = 0
        global coord₀
        v_hist         = zeros(1000)
        p_hist         = zeros(1000)
        g_hist         = zeros(1000)
        historia       = zeros(1000,4)
        global T       = zeros(size(a))
        global T[bcdof_bot_o[bcdof_bot_o .% 2 .==0]] .= -1.0
        global T[bcdof_top_o[bcdof_top_o .% 2 .==0]] .=  1.0
        g₁             = 0.0
        g₂             = 0.0
    while change > tol && OptIter < 200 #|| OptIter < 3
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
            #global traction
        # # # # # # # # # # # # # #
        OptIter += 1
        true_iteration +=1
        # detta ska 100% vara en rutin
        if true_iteration % 10 == 0 #|| OptIter == 1
            @save "innan_remeshh.jld2" h dh coord enod register Γs nₛ Γm nₘ contact_dofs contact_nods order freec_dofs free_d locked_d bcdof_o bcval_o d dh0 coord₀
            #@load "innan_remeshh.jld2"
            #break
            print("\n", " ⏳ --------  Remeshing -------- ⏳ ", "\n")
            reMeshGrids!(h, dh, coord, enod, register, Γs, nₛ, Γm, nₘ, contact_dofs, contact_nods, order, freec_dofs, free_d, locked_d, bcdof_o, bcval_o, d, dh0, coord₀)
            # Initialize tangents, också en rutin
            global K      = create_sparsity_pattern(dh) # behövs
            global Kψ     = create_sparsity_pattern(dh) # behövs
            global a      = zeros(dh.ndofs.x) # behövs
            global Ψ      = zeros(dh.ndofs.x) # behövs
            global Fᵢₙₜ  = zeros(dh.ndofs.x) # behövs?
            global rc     = zeros(dh.ndofs.x) # behövs?
            global Fₑₓₜ  = zeros(dh.ndofs.x) # behövs ?
            global a      = zeros(dh.ndofs.x) # behövs ?
            global d      = zeros(dh.ndofs.x)
            global Δa     = zeros(dh.ndofs.x) # behövs inte
            global res    = zeros(dh.ndofs.x) # behövs inte
            global ∂rᵤ_∂x = similar(K) # behövs inte om vi har lokal funktion?
            global dr_dd  = similar(K) # behövs inte om vi har lokal funktion?
            global ∂rψ_∂d = similar(K) # behövs inte om vi har lokal funktion?
            global ∂g_∂x  = zeros(size(a)) # behövs inte om vi har lokal funktion?
            global ∂g_∂u  = zeros(size(d)) # behövs inte om vi har lokal funktion?
            global ∂g₂_∂x = zeros(size(a)) # behövs inte om vi har lokal funktion?
            global ∂g₂_∂u = zeros(size(d)) # behövs inte om vi har lokal funktion?
            global ∂g₂_∂d = zeros(size(d)) # behövs inte om vi har lokal funktion?
            global ∂rᵤ_∂x = similar(K) # behövs inte om vi har lokal funktion?
            global λᵤ     = similar(a) # behövs inte om vi har lokal funktion?
            global λψ     = similar(a) # behövs inte om vi har lokal funktion?
            global λᵥₒₗ  = similar(a) # behövs inte om vi har lokal funktion?
            #
            include("initOptLin.jl")
            #
            #=
                global m = 1 # behöver inte skrivas över
                global n_mma = length(d) # behöver skrivas över
                global epsimin = 0.0000001 # behöver inte skrivas över
                global xval = d[:] # behöver skrivas över
                global xold1 = xval # behöver skrivas över
                global xold2 = xval # behöver skrivas över
                global xmin = -ones(n_mma) / 20 # behöver skrivas över
                global xmax =  ones(n_mma)  / 20 # behöver skrivas över
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
                global xmin .= -0.1
                global xmax .=  0.1
                global low           = -ones(n_mma);
                global upp           =  ones(n_mma);
            =#
            #  boundary conditions for contact analysis
            bcdof_top_o, _ = setBCXY(0.0, dh, n_top)
            bcdof_bot_o, _ = setBCXY(0.0, dh, n_bot)
            bcdof_o = [bcdof_top_o; bcdof_bot_o]
            ϵᵢⱼₖ = sortperm(bcdof_o)
            global bcdof_o = bcdof_o[ϵᵢⱼₖ]
            global bcval_o = bcdof_o .* 0.0
            # fictitious bcs
            bcdof_top_o2, _ = setBCXY(0.0, dh, n_top)
            bcdof_bot_o2, _ = setBCXY(0.0, dh, n_bot)
            bcdof_o2 = [bcdof_top_o2; bcdof_bot_o2]
            ϵᵢⱼₖ = sortperm(bcdof_o)
            global bcdof_o2 = bcdof_o2[ϵᵢⱼₖ]
            global bcval_o2 = bcdof_o2 .* 0.0
            #
            global T       = zeros(size(a))
            global T[bcdof_bot_o[bcdof_bot_o .% 2 .==0]] .= -1.0
            global T[bcdof_top_o[bcdof_top_o .% 2 .==0]] .=  1.0
            # # # # # # # # # # # # # # # #
            # Reset Optimization problem  #
            # # # # # # # # # # # # # # # #
            OptIter = 1
            global d    .= 0
            global xold1 = d[:]
            global xold2 = d[:]
        end
        #
        if OptIter % 10 == 0 && g₁ < 0.0 #  && g₂ < 0.0
            dh0 = deepcopy(dh)
            global d          = zeros(dh.ndofs.x)
            global xold1      = d[:]
            global xold2      = d[:]
            global low        = xmin
            global upp        = xmax
            OptIter           = 1
        end
        #
        # # # # #
        # test  #
        # # # # #
        global nloadsteps = 20
        global μ = 1e4 # var μ = 1e4
        #
        # # # # # # # # # # # # # #
        # Fictitious equillibrium #
        # # # # # # # # # # # # # #
        global coord₀ = getCoord(getX(dh0), dh0) # x₀
        Ψ, _, Kψ, _, λ = fictitious_solver_with_contact(d, dh0, coord₀, nloadsteps)
        # # # # # #
        # Filter  #
        # # # # # #
        #
        global dh    = deepcopy(dh0)
        updateCoords!(dh, Ψ) # x₀ + Ψ = x
        global coord = getCoord(getX(dh), dh)
        #
        # # # # #
        # test  #
        # # # # #
        global nloadsteps = 10
        global ε = 1e4 # eller?
        #
        # # # # # # # # #
        # Equillibrium  #
        # # # # # # # # #
        #
        a, _, Fₑₓₜ, Fᵢₙₜ, K, traction = solver_C(dh, coord, Δ, nloadsteps)
        # # # # # # # # #
        # Sensitivities #
        # # # # # # # # #
        #
        ∂rᵤ_∂x = similar(K)
        ∂rᵤ_∂x = drᵤ_dx_c(∂rᵤ_∂x, dh, mp, t, a, coord, enod, ε)
        dr_dd  = drψ(dr_dd, dh0, Ψ, λ, d, Γ_robin, coord₀)
        #
        # # # # # # #
        # Objective #
        # # # # # # #
        # Max reaction force
        #
        g     = - T' * Fᵢₙₜ
        ∂g_∂x =  -T' * ∂rᵤ_∂x #
        ∂g_∂u =  -T' * K # ?
        #
        # # # # # # #
        # Adjoints  #
        # # # # # # #
        #
        solveq!(λᵤ, K',  ∂g_∂u, bcdof_o, bcval_o)
        solveq!(λψ, Kψ', ∂g_∂x' - ∂rᵤ_∂x' * λᵤ, bcdof_o2, bcval_o2)
        #
        # # # # # # # # # # #
        # Full sensitivity  #
        # # # # # # # # # # #
        #
        ∂g_∂d            = (-transpose(λψ) * dr_dd)'
        #
        # # # # # # # # # # #
        # Volume constraint #
        # # # # # # # # # # #
        #
        g₁    = volume(dh,coord,enod) / Vₘₐₓ - 1.0
        ∂Ω_∂x = volume_sens(dh,coord)
        solveq!(λᵥₒₗ, Kψ, ∂Ω_∂x, bcdof_o2, bcval_o2.*0);
        ∂Ω∂d  = Real.( -transpose(λᵥₒₗ)*dr_dd ./ Vₘₐₓ) ;
        # # # # # # # # # # # #
        # Pressure constraint #
        # # # # # # # # # # # #
        #=
            p = 2
            X_ordered  = getXfromCoord(coord)
            g₂         = contact_pnorm_s(X_ordered, a, ε, p) / 10.0 - 1.0
            ∂g₂_∂x     = ForwardDiff.gradient(x -> contact_pnorm_ordered_s(x, a, ε, p), getXinDofOrder(dh, X_ordered, coord))
            ∂g₂_∂u     = ForwardDiff.gradient(u -> contact_pnorm_s(X_ordered, u, ε, p), a)
            solveq!(λᵤ, K',  ∂g₂_∂u, bcdof_o, bcval_o)
            solveq!(λψ, Kψ', ∂g₂_∂x - ∂rᵤ_∂x' * λᵤ, bcdof_o2, bcval_o2)
            ∂g₂_∂d            = Real.( (-transpose(λψ) * dr_dd)' ./ 10.0 )'
            # fulfix
            g₂      = -1.
            ∂g₂_∂d .= 0.
        =#
        # # # # # # # # # # #
        # Lås horisontellt  # // # Dålig lösning?
        # # # # # # # # # # #
        #
            # ∂g_∂d[1:2:end-1]  .= 0.0
            # ∂Ω∂d[1:2:end-1]   .= 0.0
            # ∂g₂_∂d[1:2:end-1] .= 0.0
        # # # # #
        # M M A #
        # # # # #
        #
        #d_new, ymma, zmma, lam, xsi, eta, mu, zet, S, low, upp = mmasub(m, n_mma, OptIter, d[:], xmin[:], xmax[:], xold1[:], xold2[:], g, ∂g_∂d, hcat([g₁.*100; g₂]), vcat([∂Ω∂d.*100; ∂g₂_∂d]), low, upp, a0, am, C, d2)
        d_new, ymma, zmma, lam, xsi, eta, mu, zet, S, low, upp = mmasub(m, n_mma, OptIter, d[:], xmin[:], xmax[:], xold1[:], xold2[:], g, ∂g_∂d, g₁.*100, ∂Ω∂d.*100, low, upp, a0, am, C, d2)
        xold2  = xold1
        xold1  = d
        d      = d_new
        change = norm(d .- xold1)
        #
        # # # # # # # # # #
        # Postprocessing  #
        # # # # # # # # # #
        #
        v_hist[true_iteration] = g₁
        p_hist[true_iteration] = g₂
        g_hist[true_iteration] = g
        println("Iter: ", true_iteration, " Norm of change: ", kktnorm, " Objective: ", g)
        println("Objective: ", g_hist[true_iteration], " Constraint: ", v_hist[true_iteration] , p_hist[true_iteration])
        # - - - - - -  #
        # Write to VTK #
        # - - - - - -  #
        coord = getCoord(getX(dh0), dh0)
        postprocess_opt(Ψ, dh0, "results/Current design" * string(true_iteration))
        postprocess_opt(d, dh0, "results/design_variables" * string(true_iteration))
        postprocess_opt(∂g_∂d,dh,"results/🛸" * string(true_iteration))
        # - - - - - - - - - - - - - - -  #
        # Plot objective and constraints #
        # - - - - - - - - - - - - - - -  #
        p2 = plot(1:true_iteration,[v_hist[1:true_iteration].*100,p_hist[1:true_iteration],g_hist[1:true_iteration]],label = ["Volume Constraint" "Uniform pressure Constraint" "Objective"])
        display(p2)
    end
#    jld2save("250_iter_circle.jld2")
    return g_hist, v_hist, OptIter, traction, historia
end

g_hist, v_hist, OptIter, traction, historia = Optimize(dh)
