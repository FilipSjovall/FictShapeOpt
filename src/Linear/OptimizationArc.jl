using Mortar2D, ForwardDiff, Ferrite, FerriteGmsh, FerriteMeshParser
using LinearSolve, SparseArrays, IterativeSolvers, IncompleteLU
using SparseDiffTools, Plots, Printf, JLD2, Statistics, AlgebraicMultigrid
#
plotlyjs()
#
include("..//mesh_reader.jl")
include("Contact//contact_help.jl")
include("assemLin.jl")
include("assemElemLin.jl")
include("..//material.jl")
include("..//fem.jl")
include("run_linear.jl")
include("sensitivitiesLin.jl")
include("..//mma.jl")
# ------------------- #
# Geometry parameters #
# ------------------- #
th = 0.2
x₀ = 0.0
y₀ = 0.0
Δx = 1.0
Δy = 0.25
# grid size
h = 0.025
# # # # # # # # # #
# Finite element  #
# # # # # # # # # #
ip = Lagrange{2,RefTetrahedron,1}()
qr = QuadratureRule{2,RefTetrahedron}(3)
qr_face = QuadratureRule{1,RefTetrahedron}(2)
cv = CellVectorValues(qr, ip)
fv = FaceVectorValues(qr_face, ip)
# # # # # # # # #
# Create grids  #
# # # # # # # # #
grid1 = createArcMesh("mesh_1", x₀, y₀, Δx, Δy, th, h);
Γ_1 = getBoundarySet(grid1);

# ------------------------------------------- #
# Create dofhandler with displacement field u #
# ------------------------------------------- #
global dh = DofHandler(grid1);
add!(dh, :u, 2);
close!(dh);
grid1 = nothing;
# Extract CALFEM-style matrices
global coord, enod = getTopology(dh);
global register = index_nod_to_grid(dh, coord);
# Exrtact full boundary
Γ_all = Ferrite.__collect_boundary_faces(dh.grid);
addfaceset!(dh.grid, "Γ_all", Γ_all);
Γ_all = getfaceset(dh.grid, "Γ_all");
#
n_all = getBoundarySet(dh.grid, Γ_all);
addnodeset!(dh.grid, "n_all", n_all);
#
Γ_all_dofs = Vector{Int64}()
#
# ---------------
# Displacement bc boundary u(x) = Δ ∀ x ∈ Γ_Δ
addfaceset!(dh.grid, "Γ_right", x -> x[1] ≈ x₀ + Δx)
Γ_right = getfaceset(dh.grid, "Γ_right")
#
addnodeset!(dh.grid, "n_right", x -> x[1] ≈ x₀ + Δx)
n_right = getnodeset(dh.grid, "n_right")
# addfaceset!(dh.grid, "Γ_right", x -> (x[1] ≈ x₀ + Δx && x[2] ≈ x[2] + Δy + th))
# Γ_right = getfaceset(dh.grid, "Γ_right")
# #
# addnodeset!(dh.grid, "n_right", x -> (x[1] ≈ x₀ + Δx && x[2] ≈ x[2] + Δy + th))
# n_right = getnodeset(dh.grid, "n_right")

# -------------------------------------------- #
# Displacement bc boundary u(x) = 0 ∀ x ∈ Γ_0 #
# ------------------------------------------- #
addfaceset!(dh.grid, "Γ_left", x -> x[1] ≈ x₀)
Γ_left = getfaceset(dh.grid, "Γ_left")
#
addnodeset!(dh.grid, "n_left", x -> x[1] ≈ x₀)
n_left = getnodeset(dh.grid, "n_left")
#addfaceset!(dh.grid, "Γ_left", x -> (x[1] ≈ x₀ && x[2] ≈ y₀))
#Γ_left = getfaceset(dh.grid, "Γ_left")
##
#addnodeset!(dh.grid, "n_left", x -> (x[1] ≈ x₀ && x[2] ≈ y₀))
#n_left = getnodeset(dh.grid, "n_left")


# ----------------- #
# Design boundaries #
# ----------------- #
# Γ_robin = setdiff(Γ_all, union(Γ_left, Γ_right, Γm, Γs))
#Γ_robin = setdiff(Γ_all, union(Γ_left, Γ_right, Γ_bot, Γ_top))
Γ_robin = setdiff(Γ_all, union(Γ_left, Γ_right))
addfaceset!(dh.grid, "Γ_robin", Γ_robin)

n_robin = getBoundarySet(dh.grid, Γ_robin)
addnodeset!(dh.grid, "n_robin", n_robin)
# # # # # # # # # # # # #
# Collect contact dofs  #
# # # # # # # # # # # # #
global contact_dofs = Vector{Int64}()
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
global ∂g₂_∂x = zeros(size(a)) # behövs inte om vi har lokal funktion?
global ∂g₂_∂u = zeros(size(d)) # behövs inte om vi har lokal funktion?
global λᵤ = similar(a)
global λψ = similar(a)
global Δ = -0.25
global nloadsteps = 10
global a_hist = zeros(dh.ndofs.x, nloadsteps)
global Ψ_hist = zeros(dh.ndofs.x, nloadsteps)
global d_hist = zeros(dh.ndofs.x, nloadsteps)
global F_tar = [-0.02, -0.04, -0.06, -0.08, -0.1, -0.12, -0.14, -0.16, -0.18, -0.20] * 8
global F_tar[5:end] .= F_tar[5]
#global F_tar[10] = F_tar[10] .* 1.5
global F_d = zeros(10)
global F₀ = zeros(10)
global g = 0.0
# # # # # # # # # # # # # # # #
# Init optimization variables #
# # # # # # # # # # # # # # # #
include("initHookNoContact.jl")
# ------------------- #
# Boundary conditions #
# ------------------- #
bcdof_left, _ = setBCXY_X(0.0, dh, n_left)
bcdof_right, _ = setBCXY_Y(0.0, dh, n_right)
#bcdof_bot, _   = setBCY(0.0, dh, n_bot)
#bcdof_top, _   = setBCY(0.0, dh, n_top)

#bcdof_contact, _ = setBCXY_both(0.0, dh, union(nₘ,nₛ))
bcdof_contact, _ = Vector{Int64}(), Vector{Float64}()

bcdof_bot, _ = Vector{Int64}(), Vector{Float64}()
bcdof_top, _ = Vector{Int64}(), Vector{Float64}()

#bcdofs_opt = [bcdof_left; bcdof_right; bcdof_bot; bcdof_top];
#bcdofs_opt = [bcdof_left; bcdof_right; bcdof_bot; bcdof_top; bcdof_contact];
bcdofs_opt = [bcdof_left; bcdof_right; bcdof_contact];
ϵᵢⱼₖ = sortperm(bcdofs_opt)
global bcdofs_opt = bcdofs_opt[ϵᵢⱼₖ]
global bcval_opt = bcdofs_opt .* 0.0
global asy_counter = zeros(dh.ndofs.x, 300)

global low_hist = zeros(length(d), 300)
global upp_hist = zeros(length(d), 300)
global d_hist2 = zeros(length(d), 300)
# -------------------- #
# Optimization program #
# -------------------- #
function Optimize(dh)
    # Flytta allt nedan till init_opt?
    global dh0 = deepcopy(dh)
    global λψ = similar(a)
    global λᵤ = similar(a)
    global λᵥₒₗ = similar(a)
    Vₘₐₓ = volume(dh, coord, enod) * 1.0  # 2.0
    tol = 1e-3
    global OptIter = 0
    global true_iteration = 0
    global coord₀
    v_hist = zeros(1000)
    g_hist = zeros(1000)
    #historia = zeros(200, 4)
    global T = zeros(size(a))
    global T[bcdof_left[isodd.(bcdof_left)]]   .= 1.0
    #global T[bcdof_right[isodd.(bcdof_right)]] .= -1.0
    g₁ = 0.0
    #
    while kktnorm > tol || OptIter < 2

        global d
        global Ψ
        global a
        global Fₑₓₜ
        global K
        global Kψ
        global ∂rᵤ_∂x
        global dr_dd
        global ∂rψ_∂d
        global ∂g_∂d = zeros(dh.ndofs.x)
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
        global pdofs = bcdofs_opt
        global fdofs = setdiff(1:length(a), pdofs)
        global low
        global upp
        global traction
        # # # # # # # # # # # # # #
        global OptIter += 1
        global true_iteration += 1

        # # # # #
        # Reset #
        # # # # #
        if OptIter % 50 == 0 && true_iteration < 150
            dh0 = deepcopy(dh)
            global d = zeros(dh.ndofs.x)
            global xold1 = d[:]
            global xold2 = d[:]
            global low = xmin
            global upp = xmax
            OptIter = 1
        end

        # # # # # # # # # # # # # #
        # Fictitious equillibrium #
        # # # # # # # # # # # # # #
        global nloadsteps = 10 #10
        global μ = 2e3 #1e3
        global coord₀ = getCoord(getX(dh0), dh0) # x₀
        Ψ, _, Kψ, _, λ = fictitious_solver_hook(d, dh0, coord₀, nloadsteps)

        # # # # # # # # # # # # # # #
        # Apply filter: x₀ + Ψ = x  #
        # # # # # # # # # # # # # # #
        global dh = deepcopy(dh0)
        updateCoords!(dh, Ψ) #
        global coord = getCoord(getX(dh), dh)

        # # # # # # # # #
        # Equillibrium  #
        # # # # # # # # #
        global nloadsteps = 10
        global ε = 5e2 # funkade ok med 1e4
        a, _, Fₑₓₜ, Fᵢₙₜ, K, a_hist = solver_arc(dh, coord, Δ, nloadsteps)


        global g = 0.0
        assemGlobal!(Kψ, FΨ, dh0, mp₀, t, Ψ, coord₀, enod, λ, d, Γ_robin)
        dr_dd = drψ(dr_dd, dh0, Ψ, λ, d, Γ_robin, coord₀)

        for n = 1:nloadsteps
            a = a_hist[:, n]
            #Ψ = Ψ_hist[:,n]
            λ = 1.0 #(1/nloadsteps)*n
            assemGlobal!(K, Fᵢₙₜ, dh, mp, t, a, coord, enod)
            if n ∈ [1, 4, 6, 10]
                # # # # # # # # #
                # Sensitivities #
                # # # # # # # # #
                ∂rᵤ_∂x = similar(K)
                ∂rᵤ_∂x = drᵤ_dx(∂rᵤ_∂x, dh, mp, t, a, coord, enod)


                # # # # # # #
                # Objective #
                # # # # # # #
                # Max reaction force
                g += 0.5 * (-T' * Fᵢₙₜ - F_tar[n])^2
                ∂g_∂x = -T' * ∂rᵤ_∂x
                ∂g_∂u = -T' * K
                # Compliance
                # g            = -a[pdofs]' * Fᵢₙₜ[pdofs]
                # ∂g_∂x[fdofs] = -a[pdofs]' * ∂rᵤ_∂x[pdofs, fdofs]
                # ∂g_∂u[fdofs] = -a[pdofs]' * K[pdofs, fdofs]

                # # # # # # #
                # Adjoints  #
                # # # # # # #
                solveq!(λᵤ, K', ∂g_∂u, bcdofs_opt, bcval_opt)
                solveq!(λψ, Kψ', ∂g_∂x' - ∂rᵤ_∂x' * λᵤ, bcdofs_opt, bcval_opt)

                # # # # # # # # # # #
                # Full sensitivity  #
                # # # # # # # # # # #
                ∂g_∂d += (-T' * Fᵢₙₜ - F_tar[n]) * (-transpose(λψ) * dr_dd)'
            end
            F_d[n] = -T' * Fᵢₙₜ
        end

        if true_iteration == 1
            global F₀ = deepcopy(F_d)
        end
        # # # # # # # # # # #
        # Volume constraint #
        # # # # # # # # # # #
        g₁ = volume(dh, coord, enod) / Vₘₐₓ - 1.0
        ∂Ω_∂x = volume_sens(dh, coord)
        solveq!(λᵥₒₗ, Kψ, ∂Ω_∂x, bcdofs_opt, bcval_opt)
        ∂Ω∂d = Real.(-transpose(λᵥₒₗ) * dr_dd ./ Vₘₐₓ)

        # # # # #
        # M M A #
        # # # # #
        # d_old   = d
        # low_old = low
        # upp_old = upp
        # d_new, ymma, zmma, lam, xsi, eta, mu, zet, S, low, upp = mmasub(m, n_mma, OptIter, d[:], xmin[:], xmax[:], xold1[:], xold2[:], g .* 1000, ∂g_∂d .* 1000, g₁ .* 100, ∂Ω∂d .* 100, low, upp, a0, am, C, d2)
        # #d_new, ymma, zmma, lam, xsi, eta, mu, zet, S, low, upp = mmasub(m, n_mma, OptIter, d, xmin, xmax, xold1, xold2, g .* 100, ∂g_∂d .* 100, hcat([g₁; g₂]), vcat([∂Ω∂d; ∂g₂_∂d]), low, upp, a0, am, C, d2)
        # # ----------------- #
        # # Test - new update #
        # # ----------------- #
        # α      = 1.0 # 0.4 # 0.1 #
        # d_new  = d_old   + α .* (d_new - d_old)
        # low    = low_old + α .* (low   - low_old)
        # upp    = upp_old + α .* (upp   - upp_old)
        # # ----------------- #
        # xold2  = xold1
        # xold1  = d
        # d      = d_new
        # change = norm(d .- xold1)
        d_old = d[free_d]
        low_old = low
        upp_old = upp
        d_new, ymma, zmma, lam, xsi, eta, mu, zet, S, low, upp = mmasub(m, n_mma, OptIter, d[free_d], xmin[:], xmax[:], xold1[:], xold2[:], g .* 1000, ∂g_∂d[free_d] .* 1000, g₁ .* 100, ∂Ω∂d[free_d]' .* 100, low, upp, a0, am, C, d2)
        # ----------------- #
        # Test - new update #
        # ----------------- #
        α = 0.5 # 0.4 # 0.1 #
        d_new = d_old + α .* (d_new - d_old)
        low = low_old + α .* (low - low_old)
        upp = upp_old + α .* (upp - upp_old)
        # ----------------- #
        xold2 = xold1
        xold1 = d[free_d]
        d[free_d] = d_new
        change = norm(d[free_d] .- xold1)

        # # # # # # # # # #
        # Postprocessing  #
        # # # # # # # # # #
        v_hist[true_iteration] = g₁
        g_hist[true_iteration] = g

        #The residual vector of the KKT conditions is calculated:
        #residu,kktnorm,residumax = kktcheck(m,n,X,ymma,zmma,lam,xsi,eta,mu,zet,S, xmin,xmax,∂g_∂d,[0.0],zeros(size(d)),a0,a,C,d2);
        kktnorm = change
        println("Iter: ", true_iteration, " Norm of change: ", kktnorm, " Objective: ", g)
        postprocess_opt(Ψ, dh0, "results/Current design" * string(true_iteration))
        postprocess_opt(d, dh0, "results/design_variables" * string(true_iteration))
        #postprocess_opt(∂g_∂d, dh, "results/🛸" * string(true_iteration))
        println("Objective: ", g_hist[1:true_iteration], " Constraint: ", v_hist[1:true_iteration])
        # append?
        p2 = plot(1:true_iteration, v_hist[1:true_iteration], label="Volume", background_color=RGB(0.2, 0.2, 0.2), legend=:outerleft, grid=false, lc=:orange)
        p3 = plot(1:true_iteration, g_hist[1:true_iteration] .* 1000, label="Objective", background_color=RGB(0.2, 0.2, 0.2), legend=:outerleft, grid=false, lc=:purple)
        #display(p2)
        p4 = plot(0.0:0.01:0.1, vcat(0.0, abs.(F_d)), label="Design", background_color=RGB(0.2, 0.2, 0.2), legend=:outerleft, grid=false) # ,
        plot!(0.0:0.01:0.1, vcat(0.0, abs.(F₀)), label="Initial")
        plot!([0.0 0.01 0.04 0.06 0.1]', hcat(0.0, abs.([F_tar[1] F_tar[4] F_tar[6] F_tar[10]]))', label="Target", marker=:x, color=:crimson, linestyle=:dash)
        #display(p3)
        p = plot(p2, p3, p4, layout=(3, 1), size=(800, 600))
        display(p)
        # For investigative purpose
        low_hist[free_d, true_iteration] = low
        upp_hist[free_d, true_iteration] = upp
        d_hist2[free_d, true_iteration] = d[free_d]
        @save "asymptoter.jld2" low_hist upp_hist d_hist
    end
    #jld2save("färdig.jld2",a,dh,dh0,Opiter,v_hist,g_hist,d)
    return g_hist, v_hist, OptIter, traction
end

g_hist, v_hist, OptIter, traction = Optimize(dh)


# # # # # #
# Plot 3D #
# # # # # #

#
# Plotta bara "free_d"
design_indices = free_d
opt_iters = range(1, 100)
#
#
# Create x and y coordinates based on the matrix size
x = 1:size(low_hist[design_indices, opt_iters], 2)
y = 1:size(low_hist[design_indices, opt_iters], 1)

# Create a meshgrid from x and y
xgrid = repeat(x, 1, length(y))
ygrid = repeat(y, 1, length(x))'

# Create the 3D surface plot
plot(xgrid, ygrid, low_hist[design_indices, opt_iters]', st=:surface, xlabel="Design iterations", ylabel="Degree of freedom", zlabel="Asymptote value", title="Evolution of asymptotes")
plot!(xgrid, ygrid, upp_hist[design_indices, opt_iters]', st=:surface)


F = zeros(10)

for n = 1:nloadsteps
    a = a_hist[:, n]
    assemGlobal!(K, Fᵢₙₜ, rc, dh, mp, t, a, coord, enod, ε)
    assemGlobal!(Kψ, FΨ, dh0, mp₀, t, Ψ, coord₀, enod, λ, d, Γ_robin, μ)
    F[n] = -T' * Fᵢₙₜ
end
plot(0.01:0.01:0.1, abs.(F))





# For abstract
