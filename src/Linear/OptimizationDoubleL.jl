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


xl = 0.0
yl = 0.0
xr = -.75
yr = 1.45
Δx = 1.25
Δy = 1.0
th = 0.25
r1 = 0.05
r2 = 0.05

# # # # # # # # # #
# Finite element  #
# # # # # # # # # #
ip = Lagrange{2,RefTetrahedron,1}()
qr = QuadratureRule{2,RefTetrahedron}(1)
qr_face = QuadratureRule{1,RefTetrahedron}(1)
cv = CellVectorValues(qr, ip)
fv = FaceVectorValues(qr_face, ip)

# # # # # # # # #
# Create grids  #
# # # # # # # # #
grid1 = createLMesh("mesh_1", xl, yl, Δx, Δy, th, r1, r2, 0.05);
Γ_1   = getBoundarySet(grid1);
grid2 = createLMeshRev("mesh_2", xr, yr, Δx, Δy, th, r1, r2, 0.05);
Γ_2   = getBoundarySet(grid2);
grid_tot = merge_grids(grid1, grid2; tol=1e-6);
grid1 = nothing;
grid2 = nothing;
# Create dofhandler with displacement field u
global dh = DofHandler(grid_tot);
add!(dh, :u, 2);
close!(dh);

# Extract CALFEM-style matrices
global coord, enod = getTopology(dh);
global register = index_nod_to_grid(dh, coord);


# Exrtact full boundary
Γ_all   = Ferrite.__collect_boundary_faces(dh.grid);
addfaceset!(dh.grid,"Γ_all", Γ_all);
Γ_all  = getfaceset(dh.grid, "Γ_all");

n_all = getBoundarySet(dh.grid, Γ_all);
addnodeset!(dh.grid, "n_all", n_all);

Γ_all_dofs = Vector{Int64}()
# --------------
# Master
addfaceset!(dh.grid, "Γ_master", x -> x ∈ Γ_2 );
Γm = getfaceset(dh.grid, "Γ_master");
Γm =  intersect(Γm, Γ_all);

nₘ = getBoundarySet(dh.grid,Γm);
addnodeset!(dh.grid,"nₘ" ,nₘ);

# ---------------
# Slave
addfaceset!(dh.grid,"Γ_slave", x -> x ∈ Γ_1 );
Γs = getfaceset(dh.grid, "Γ_slave");
Γs =  intersect(Γs, Γ_all);

nₛ = getBoundarySet(dh.grid,Γs)
addnodeset!(dh.grid, "nₛ" ,nₛ)


# ---------------
# Displacement bc boundary u(x) = Δ ∀ x ∈ Γ_Δ
addfaceset!(dh.grid, "Γ_right", x -> x[1] ≈ xl + Δx)
Γ_right = getfaceset(dh.grid, "Γ_right")

addnodeset!(dh.grid, "n_right", x -> x[1] ≈ xl + Δx)
n_right = getnodeset(dh.grid, "n_right")

# --------------
# Displacement bc boundary u(x) = 0 ∀ x ∈ Γ_0
addfaceset!(dh.grid, "Γ_left", x -> x[1] ≈ xr )
Γ_left = getfaceset(dh.grid, "Γ_left")

addnodeset!(dh.grid, "n_left", x -> x[1] ≈ xr )
n_left = getnodeset(dh.grid, "n_left")

# --------------
# bottom
addfaceset!(dh.grid, "Γ_bot", x -> x[2] ≈ yl)
Γ_bot = getfaceset(dh.grid, "Γ_bot")

addnodeset!(dh.grid, "n_bot", x -> x[2] ≈ yl)
n_bot = getnodeset(dh.grid, "n_bot")

# --------------
# Top
addfaceset!(dh.grid, "Γ_top", x -> x[2] ≈ yr)
Γ_top = getfaceset(dh.grid, "Γ_top")

addnodeset!(dh.grid, "n_top", x -> x[2] ≈ yr)
n_top = getnodeset(dh.grid, "n_top")

# ---------------
# Design boundaries
#Γ_robin = setdiff(Γ_all, union(Γ_left, Γ_right, Γ_bot, Γ_top))
Γ_robin = setdiff(Γ_all, union(Γ_left, Γ_right))
addfaceset!(dh.grid, "Γ_robin", Γ_robin)

n_robin = getBoundarySet(dh.grid,Γ_robin)
addnodeset!(dh.grid, "n_robin", n_robin)


# # # # # # # # # # # # #
# Collect contact dofs  #
# # # # # # # # # # # # #
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

# # # # # # # # # # # #
# Collect design dofs #
# # # # # # # # # # # #
global free_d = []
for jnod in n_robin
    append!(free_d, register[jnod, 1])
    append!(free_d, register[jnod, 2])
end
##global locked_d = setdiff(1:dh.ndofs.x, free_d)
#global locked_d = Vector{Int64}()
#for n ∈ n_left
#    push!(locked_d,register[n,1])
#    push!(locked_d,register[n,2])
#end
#for n ∈ n_right
#    push!(locked_d,register[n,1])
#    push!(locked_d,register[n,2])
#end

# Initialize tangents
global K      = create_sparsity_pattern(dh)
global Kψ     = create_sparsity_pattern(dh)
global a      = zeros(dh.ndofs.x)
global d      = zeros(dh.ndofs.x)
global Ψ      = zeros(dh.ndofs.x)
global Fᵢₙₜ  = zeros(dh.ndofs.x)
global rc     = zeros(dh.ndofs.x)
global Fₑₓₜ  = zeros(dh.ndofs.x)
global a      = zeros(dh.ndofs.x)
global Δa     = zeros(dh.ndofs.x)
global res    = zeros(dh.ndofs.x)
global dr_dd  = similar(K)
global ∂rψ_∂d = similar(K)
global ∂g_∂x  = zeros(size(a)) # behövs inte om vi har lokal funktion?
global ∂g_∂u  = zeros(size(d)) # behövs inte om vi har lokal funktion?
global ∂g₂_∂x = zeros(size(a)) # behövs inte om vi har lokal funktion?
global ∂g₂_∂u = zeros(size(d)) # behövs inte om vi har lokal funktion?
global λᵤ     = similar(a)
global λψ     = similar(a)

global Δ = 0.05
global nloadsteps = 10

include("initOptLinHook.jl")

# ------------------- #
# Boundary conditions #
# ------------------- #
bcdof_left, _    = setBCXY_X(0.0, dh, n_left)
bcdof_right, _   = setBCXY_X(0.0, dh, n_right)
bcdof_bot, _     = setBCY(0.0, dh, n_bot)
bcdof_top, _     = setBCY(0.0, dh, n_top)

bcdof_bot, _     = Vector{Int64}(), Vector{Float64}()
bcdof_top, _     = Vector{Int64}(), Vector{Float64}()

bcdofs_opt       = [bcdof_left; bcdof_right; bcdof_bot; bcdof_top];
ϵᵢⱼₖ            = sortperm(bcdofs_opt)
global bcdofs_opt= bcdofs_opt[ϵᵢⱼₖ]
global bcval_opt = bcdofs_opt .* 0.0

function Optimize(dh)
    # Flytta allt nedan till init_opt?
        global dh0     = deepcopy(dh)
        global λψ      = similar(a)
        global λᵤ      = similar(a)
        global λᵥₒₗ   = similar(a)
        Vₘₐₓ          = 1.2  #
        tol            = 1e-3
        OptIter        = 0
        true_iteration = 0
        global coord₀
        v_hist         = zeros(1000)
        g_hist         = zeros(1000)
        historia       = zeros(200,4)
        global T       = zeros(size(a))
        global T[bcdof_left[isodd.(bcdof_left)]] .= 1.0
        g₁             = 0.0
        g₂             = 0.0
    #
    while kktnorm > tol || OptIter < 200

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
            global pdofs    = bcdofs_opt
            global fdofs    = setdiff(1:length(a), pdofs)
            #global locked_d = setdiff(1:length(a),free_d)
            global low
            global upp
            global traction
        # # # # # # # # # # # # # #
        OptIter        += 1
        true_iteration += 1

        # # # # #
        # test  #
        # # # # #
        global nloadsteps = 20
        global μ          = 1e3 # var μ = 1e4

        if OptIter % 10 == 0 && g₁ < 0
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
        global coord₀  = getCoord(getX(dh0), dh0) # x₀
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
        global ε          = 1e4 # eller?

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

        # # # # # # # # # # # #
        # Pressure constraint #
        # # # # # # # # # # # #
        #p = 2
        #X_ordered = getXfromCoord(coord)
        #g₂ = contact_pnorm_s(X_ordered, a, ε, p) / 0.5 - 1.0
        #∂g₂_∂x = ForwardDiff.gradient(x -> contact_pnorm_ordered_s(x, a, ε, p), getXinDofOrder(dh, X_ordered, coord))
        #∂g₂_∂u = ForwardDiff.gradient(u -> contact_pnorm_s(X_ordered, u, ε, p), a)
#
        #solveq!(λᵤ, K', ∂g₂_∂u, bcdofs_opt, bcval_opt)
        #solveq!(λψ, Kψ', ∂g₂_∂x - ∂rᵤ_∂x' * λᵤ, bcdofs_opt, bcval_opt)
        #∂g₂_∂d = Real.((-transpose(λψ) * dr_dd)' ./ 0.5)'

        # # # # #
        # M M A #
        # # # # #
        d_new, ymma, zmma, lam, xsi, eta, mu, zet, S, low, upp = mmasub(m, n_mma, OptIter, d, xmin, xmax, xold1, xold2, g.*10, ∂g_∂d.*10, g₁.*100, ∂Ω∂d.*100, low, upp, a0, am, C, d2)
        #d_new, ymma, zmma, lam, xsi, eta, mu, zet, S, low, upp = mmasub(m, n_mma, OptIter, d, xmin, xmax, xold1, xold2, g .* 100, ∂g_∂d .* 100, hcat([g₁; g₂]), vcat([∂Ω∂d; ∂g₂_∂d]), low, upp, a0, am, C, d2)
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
            postprocess_opt(Ψ, dh0, "results/Current design"   * string(true_iteration))
            postprocess_opt(d, dh0, "results/design_variables" * string(true_iteration))

            postprocess_opt(∂g_∂d, dh0, "test")
        end

        println("Objective: ", g_hist[1:true_iteration], " Constraint: ", v_hist[1:true_iteration])


        p2 = plot(1:true_iteration,[v_hist[1:true_iteration], g_hist[1:true_iteration] .*100],label = ["Volume Constraint" "Objective"])
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
