using Mortar2D, ForwardDiff, Ferrite, FerriteGmsh, FerriteMeshParser
using LinearSolve, SparseArrays, IterativeSolvers, IncompleteLU
using SparseDiffTools, Plots, Printf, JLD2, Statistics

include("..//mesh_reader.jl")
include("Contact//contact_help.jl")
include("assemLin.jl")
include("assemElemLin.jl")
include("..//material.jl")
include("..//fem.jl")
include("run_linear.jl")
include("sensitivitiesLin.jl")
include("..//mma.jl")
# # # # # # #
# Geometry  #
# # # # # # #
xₗ = 0.0
xᵤ = 0.5
Δx = 1.0
yₗ = 0.0
yᵤ = 0.501
Δy = 0.5
# # # # # # # # # #
# Finite element  #
# # # # # # # # # #
ip      = Lagrange{2,RefTetrahedron,1}()
qr      = QuadratureRule{2,RefTetrahedron}(1)
qr_face = QuadratureRule{1,RefTetrahedron}(1)
cv      = CellVectorValues(qr, ip)
fv      = FaceVectorValues(qr_face, ip)

# # # # # # # # #
# Create grids  #
# # # # # # # # #
grid1 = createBoxMeshRev("box_1", xₗ, yₗ, Δx, Δy, 0.051)
grid2 = createBoxMeshRev("box_2", xᵤ, yᵤ, Δx, Δy, 0.05)
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

# ----------------- #
# Create slave sets #
# ----------------- #
#addfaceset!(dh.grid,"Γ_slave", x-> ( (x[2] ≈ yᵤ) || (( x[2]≥yᵤ) && (x[1] ≈ xᵤ) ) ))
addfaceset!(dh.grid,"Γ_slave", x-> ( (x[2] ≈ yᵤ) ))
global Γs = getfaceset(dh.grid, "Γ_slave")

#addnodeset!(dh.grid, "nₛ", x -> ( ( x[2] ≈ yᵤ) || ( (x[2] ≥ yᵤ) && (x[1] ≈ xᵤ) ) ) )
addnodeset!(dh.grid, "nₛ", x -> ( ( x[2] ≈ yᵤ) ))
global nₛ = getnodeset(dh.grid, "nₛ")

# ------------------ #
# Create master sets #
# ------------------ #
#addfaceset!(dh.grid, "Γ_master", x -> ( ( x[2] ≈ yₗ + Δy) || ( (x[2] ≤ yₗ + Δy) && (x[1] ≈ xₗ) ) ) )
addfaceset!(dh.grid, "Γ_master", x -> ( ( x[2] ≈ yₗ + Δy) ))
global Γm = getfaceset(dh.grid, "Γ_master")

#addnodeset!(dh.grid, "nₘ", x ->  ( (x[2] ≈ yₗ + Δy ) || ( (x[2] ≤ yₗ + Δy) && (x[1] ≈ xₗ + Δx) ) ) )
addnodeset!(dh.grid, "nₘ", x ->  ( (x[2] ≈ yₗ + Δy ) ))
global nₘ = getnodeset(dh.grid, "nₘ" )

# ------------------ #
# Create right  sets #
# ------------------ #
addfaceset!(dh.grid, "Γ_right", x -> x[1] ≈ xᵤ + Δx)
global Γ_right = getfaceset(dh.grid, "Γ_right")

addnodeset!(dh.grid, "nᵣ", x -> x[1] ≈ xᵤ + Δx)
global n_right = getnodeset(dh.grid, "nᵣ")

# ------------------ #
# Create left | sets #
# ------------------ #
addfaceset!(dh.grid, "Γ_left", x -> x[1] ≈ xₗ)
global Γ_left = getfaceset(dh.grid, "Γ_left")

addnodeset!(dh.grid, "nₗ", x -> x[1] ≈ xₗ)
global n_left = getnodeset(dh.grid, "nₗ")

# --------------- #
# Create top sets #
# --------------- #
addfaceset!(dh.grid, "Γ_top", x -> x[2] ≈ yᵤ + Δy)
global Γ_top = getfaceset(dh.grid, "Γ_top")

addnodeset!(dh.grid, "n_top", x -> x[2] ≈ yᵤ + Δy)
global n_top = getnodeset(dh.grid, "n_top")

# --------------- #
# Create bot sets #
# --------------- #
addfaceset!(dh.grid, "Γ_bot", x -> x[2] ≈ yₗ)
global Γ_bot = getfaceset(dh.grid, "Γ_bot")

addnodeset!(dh.grid, "n_bot", x -> x[2] ≈ yₗ)
global n_bot = getnodeset(dh.grid, "n_bot")

# Extract all nbr nodes and dofs
global contact_dofs = getContactDofs(nₛ, nₘ)
global contact_nods = getContactNods(nₛ, nₘ)
global order = Dict{Int64,Int64}()
for (i, nod) ∈ enumerate(contact_nods)
    push!(order, nod => i)
end
global freec_dofs = setdiff(1:dh.ndofs.x, contact_dofs)

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
    getfaceset(dh.grid, "Γ_bot"),
    getfaceset(dh.grid, "Γ_top"),
    getfaceset(dh.grid, "Γ_master")
)
global n_robin = union(
    getnodeset(dh.grid, "nₛ"),
    getnodeset(dh.grid, "n_bot"),
    getnodeset(dh.grid, "n_top"),
    getnodeset(dh.grid, "nₘ")
)

global free_d = []
for jnod in n_robin
    if in(jnod, n_left) || in(jnod, n_right)
        append!(free_d, register[jnod, 1])
        append!(free_d, register[jnod, 2]) ## Fundera på om detta skall vara med
    else
        append!(free_d, register[jnod, 1])
        append!(free_d, register[jnod, 2])
    end
end
global locked_d = setdiff(1:dh.ndofs.x, free_d)

# Initialize tangents
global K = create_sparsity_pattern(dh)
global Kψ = create_sparsity_pattern(dh)
global a = zeros(dh.ndofs.x)
global d = zeros(dh.ndofs.x)
global Ψ = zeros(dh.ndofs.x)
global Fᵢₙₜ = zeros(dh.ndofs.x)
global rc = zeros(dh.ndofs.x)
global Fₑₓₜ = zeros(dh.ndofs.x)
global a = zeros(dh.ndofs.x)
global Δa = zeros(dh.ndofs.x)
global res = zeros(dh.ndofs.x)
global dr_dd = similar(K)
global ∂rψ_∂d = similar(K)
global ∂g_∂x = zeros(size(a)) # behövs inte om vi har lokal funktion?
global ∂g_∂u = zeros(size(d)) # behövs inte om vi har lokal funktion?
global λᵤ = similar(a)
global λψ = similar(a)

global Δ = -0.15
global nloadsteps = 10
global asy_counter = zeros(dh.ndofs.x, 400)
include("initOptLinHook.jl")

# ------------------- #
# Boundary conditions #
# ------------------- #
bcdof_left, _    = setBCXY_both(0.0, dh, n_left)
bcdof_right,_    = setBCXY_both(0.0, dh, n_right)
bcdofs_opt       = [bcdof_left; bcdof_right];
ϵᵢⱼₖ            = sortperm(bcdofs_opt)
global bcdofs_opt = bcdofs_opt[ϵᵢⱼₖ]
global bcval_opt = bcdofs_opt .* 0.0

function Optimize(dh)
    # Flytta allt nedan till init_opt?
        global dh0   = deepcopy(dh)
        global λψ    = similar(a)
        global λᵤ    = similar(a)
        global λᵥₒₗ  = similar(a)
        Vₘₐₓ         = 1.2  #
        tol            = 1e-3
        OptIter        = 0
        global true_iteration = 0
        global coord₀
        v_hist         = zeros(1000)
        g_hist         = zeros(1000)
        historia       = zeros(200,4)
        global T       = zeros(size(a))
        global T[bcdof_right[bcdof_right .% 2 .==0]] .= 1.0
        g₁             = 0.0
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
            global ∂g_Ωd
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
            global pdofs    = bcdofs_opt
            global fdofs    = setdiff(1:length(a), pdofs)
            global locked_d = setdiff(1:length(a),free_d)
            global low
            global upp
            global traction
        # # # # # # # # # # # # # #
        OptIter += 1
        true_iteration +=1

        # # # # #
        # test  #
        # # # # #
        global nloadsteps = 10
        global μ = 1e3 # var μ = 1e4

        if OptIter % 25 == 0 && g₁ < 0
            dh0 = deepcopy(dh)
            global d          = zeros(dh.ndofs.x)
            global xold1      = d[:]
            global xold2      = d[:]
            global low        = xmin
            global upp        = xmax
            OptIter           = 1
        end

        # # # # # # # # # # # # # #
        # Fictitious equillibrium #
        # # # # # # # # # # # # # #
        global coord₀ = getCoord(getX(dh0), dh0) # x₀
        Ψ, _, Kψ, _, λ = fictitious_solver_with_contact_hook(d, dh0, coord₀, nloadsteps)

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
        global ε = 1e5 # eller?

        # # # # # # # # #
        # Equillibrium  #
        # # # # # # # # #
        a, _, Fₑₓₜ, Fᵢₙₜ, K, traction = solver_C_hook(dh, coord, Δ, nloadsteps)

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
        ∂g_∂x = -T' * ∂rᵤ_∂x ## stämmer? innehåller kontaktkänslighet men dessa träffar bara kontaktdofs som inte är kopplade till bc.
        ∂g_∂u = -T' * K
        # Compliance
        # g            = -a[pdofs]' * Fᵢₙₜ[pdofs]
        # ∂g_∂x[fdofs] = -a[pdofs]' * ∂rᵤ_∂x[pdofs, fdofs]
        # ∂g_∂u[fdofs] = -a[pdofs]' * K[pdofs, fdofs]

        # # # # # # #
        # Adjoints  #
        # # # # # # #
        solveq!(λᵤ, K',  ∂g_∂u, bcdofs_opt, bcval_opt)
        solveq!(λψ, Kψ', ∂g_∂x' - ∂rᵤ_∂x' * λᵤ, bcdofs_opt, bcval_opt)

        # # # # # # # # # # #
        # Full sensitivity  #
        # # # # # # # # # # #
        ∂g_∂d            = (-transpose(λψ) * dr_dd)'

        # # # # # # # # # # #
        # Volume constraint #
        # # # # # # # # # # #
        g₁    = volume(dh,coord,enod) / Vₘₐₓ - 1.0
        ∂Ω_∂x = volume_sens(dh,coord)
        solveq!(λᵥₒₗ, Kψ, ∂Ω_∂x, bcdofs_opt, bcval_opt);
        ∂Ω∂d  = Real.( -transpose(λᵥₒₗ)*dr_dd ./ Vₘₐₓ) ;

        # # # # #
        # M M A #
        # # # # #
        #d_new, ymma, zmma, lam, xsi, eta, mu, zet, S, low, upp = mmasub(m, n_mma, OptIter, d, xmin, xmax, xold1, xold2, g , ∂g_∂d , g₁, ∂Ω∂d, low, upp, a0, am, C, d2)
        d_new, ymma, zmma, lam, xsi, eta, mu, zet, S, low, upp = mmasub(m, n_mma, OptIter, d[:], xmin[:], xmax[:], xold1[:], xold2[:], g.*100, ∂g_∂d.*100, g₁.*100, ∂Ω∂d.*100, low, upp, a0, am, C, d2)
        xold2  = xold1
        xold1  = d
        d      = d_new
        change = norm(d .- xold1)

        # # # # # # # # # #
        # Postprocessing  #
        # # # # # # # # # #
        v_hist[true_iteration] = g₁
        g_hist[true_iteration] = g


        #The residual vector of the KKT conditions is calculated:
        #residu,kktnorm,residumax = kktcheck(m,n,X,ymma,zmma,lam,xsi,eta,mu,zet,S, xmin,xmax,∂g_∂d,[0.0],zeros(size(d)),a0,a,C,d2);
        kktnorm = change
        println("Iter: ", true_iteration, " Norm of change: ", kktnorm, " Objective: ", g)
        if mod(OptIter,1) == 0
            coord = getCoord(getX(dh0), dh0)
            postprocess_opt(Ψ, dh0, "results/Current design" * string(true_iteration))
            postprocess_opt(d, dh0, "results/design_variables" * string(true_iteration))
        end

        println("Objective: ", g_hist[1:true_iteration], " Constraint: ", v_hist[1:true_iteration])


        p2 = plot(1:true_iteration,[v_hist[1:true_iteration], g_hist[1:true_iteration]],label = ["Volume Constraint" "Objective"])
        display(p2)

        #p3 = plot(1:true_iteration,g_hist[1:true_iteration],legend=false, marker=3, reuse = false, lc =:darkgreen)
        #display(p3)

        if true_iteration == 1
            g_ini = 0
            n     = 0
            xval  = 0
            #@save "stegett.jld2"
        elseif true_iteration == 2
            #@save "tva.jld2"
        elseif true_iteration == 200
            #@save "steg200.jld2"
            #break
        end
    end
    jld2save("färdig.jld2",a,dh,dh0,Opiter,v_hist,g_hist,d)
    return g_hist, v_hist, OptIter, traction, historia
end

# plot(coord[collect(n_robin),1], ∂g_∂d[free_d], seriestype=:scatter)
g_hist, v_hist, OptIter, traction, historia = Optimize(dh)
