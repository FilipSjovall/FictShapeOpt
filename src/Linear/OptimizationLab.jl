using Pkg
Pkg.activate()
# kolla Pkg.status() vid problem / jämför med att bara starta julia i en terminal
using ForwardDiff, Ferrite, FerriteGmsh, FerriteMeshParser
using LinearSolve, SparseArrays, IterativeSolvers, IncompleteLU
using SparseDiffTools, Plots, Printf, JLD2, Statistics, AlgebraicMultigrid
# kan behöva köra export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib
plotlyjs()
#
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
# ------------------- #
# Geometry parameters #
# ------------------- #
# - Block - #
th = 0.1
x₁ = 0.0
y₁ = 0.25001#1
Δx = 0.5
Δy = 0.1
# - Seal - #
x₀ = 0.0
y₀ = 0.0
B  = 0.15
b  = 0.1 #+ 0.0125*2
Δl = (Δx - B)  #0.05
H  = 0.15
r  = 0.025 #0.0125
r2 = 0.025
# för vertikal sida på gasket skall B/2 - b/2 - r = 0 gälla.
# grid size
h = 0.05 # 0.05
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
grid1 = createQuarterLabyrinthMeshRounded("mesh_1", x₀, y₀, th, B, b, Δl, H, r, h);
#grid1 = createQuarterLabyrinthMeshRoundedCavity("mesh_1", x₀, y₀, th, B, b, Δl, H, r, r2, h);
Γ_1 = getBoundarySet(grid1);
grid2 = createBoxMeshRev("mesh_2", x₁, y₁, Δx, Δy, h);
#grid2 = createBoxMesh("mesh_2", x₁, y₁, Δx, Δy, h);
Γ_2 = getBoundarySet(grid2);
grid_tot = merge_grids(grid1, grid2; tol=1e-8);
grid1 = nothing;
grid2 = nothing;
# ------------------------------------------- #
# Create dofhandler with displacement field u #
# ------------------------------------------- #
global dh = DofHandler(grid_tot);
add!(dh, :u, 2);
close!(dh);
# Extract CALFEM-style matrices
global coord, enod = getTopology(dh);
global register = index_nod_to_grid(dh, coord);
# - - - - - - - - #
# Create cellsets #
# - - - - - - - - #
top_mesh = addcellset!(dh.grid, "top mesh", x -> x[2] > th + H)
bot_mesh = addcellset!(dh.grid, "bot mesh", x -> x[2] ≤ th + H)
# Exrtact full boundary
Γ_all = Ferrite.__collect_boundary_faces(dh.grid);
addfaceset!(dh.grid, "Γ_all", Γ_all);
Γ_all = getfaceset(dh.grid, "Γ_all");
#
n_all = getBoundarySet(dh.grid, Γ_all);
addnodeset!(dh.grid, "n_all", n_all);
#
Γ_all_dofs = Vector{Int64}()
# ------ #
# Master #
# ------ #
addfaceset!(dh.grid, "Γ_master", x -> x[2] ≈ y₁);
Γm = getfaceset(dh.grid, "Γ_master");
Γm = intersect(Γm, Γ_all);
#
nₘ = getBoundarySet(dh.grid, Γm);
addnodeset!(dh.grid, "nₘ", nₘ);
#
# ----- #
# Slave #
# ----- #
addfaceset!(dh.grid, "Γ_slave", x ->  x ∈ Γ_1 );
Γs = getfaceset(dh.grid, "Γ_slave");
Γs = intersect(Γs, Γ_all);
#
global nₛ = getBoundarySet(dh.grid, Γs)
addnodeset!(dh.grid, "nₛ", nₛ)

# ------ #
# bottom #
# ------ #
addfaceset!(dh.grid, "Γ_bot", x -> x[2] ≈ y₀)
Γ_bot = getfaceset(dh.grid, "Γ_bot")

addnodeset!(dh.grid, "n_bot", x -> x[2] ≈ y₀)
n_bot = getnodeset(dh.grid, "n_bot")
# --- #
# Top #
# --- #
addfaceset!(dh.grid, "Γ_top", x -> x[2] ≈ y₁ + Δy)
Γ_top = getfaceset(dh.grid, "Γ_top")

addnodeset!(dh.grid, "n_top", x -> x[2] ≈ y₁ + Δy)
n_top = getnodeset(dh.grid, "n_top")


# --- #
# mid #
# --- #
# addnodeset!(dh.grid, "n_mid", x -> (x[2] ≈ H + th) && (x[1] ≈ 0.25))
# n_mid = getnodeset(dh.grid, "n_mid")

# ----------------------------- #
# left and right sides of block #
# ----------------------------- #
#addnodeset!(dh.grid, "n_lr", x -> ( x[2]≥ y₁ && (x[1] ≈ x₁ || x[1]≈ x₁ + Δx) ))
addnodeset!(dh.grid, "n_lr", x -> ( x[2]≥ y₁ && x[1] ≈ x₁  ))
nₗᵣ = getnodeset(dh.grid, "n_lr")

addfaceset!(dh.grid, "Γ_lr", x -> ( x[2]≥ y₁ && x[1] ≈ x₁  ))
Γ_lr = getfaceset(dh.grid, "Γ_lr")

# ------------------------------ #
# Middle nodes on top and bottom #
# ------------------------------ #
addnodeset!(dh.grid,"n_sym", x->x[1] ≈ 0.5)
n_sym = getnodeset(dh.grid, "n_sym")

addfaceset!(dh.grid, "Γ_sym", x->x[1] ≈ 0.5)
Γ_sym = getfaceset(dh.grid, "Γ_sym")

# -------------- #
# Left side both #
# -------------- #
addnodeset!(dh.grid,"n_left", x->x[1] ≈ x₀)
n_left = getnodeset(dh.grid, "n_left")

addfaceset!(dh.grid, "Γ_left", x->x[1] ≈ x₀)
Γ_left = getfaceset(dh.grid, "Γ_left")

# ----------------- #
# Design boundaries #
# ----------------- #
#Γ_robin = setdiff(Γ_all, union(Γ_top, Γ_bot, Γm, Γ_sym, Γ_lr))
#Γ_robin = setdiff(Γ_all, union(Γ_top, Γ_bot, Γm, Γ_sym, Γ_lr))
Γ_robin = setdiff(Γ_all, union(Γ_top, Γ_bot, Γ_sym, Γ_lr))
addfaceset!(dh.grid, "Γ_robin", Γ_robin)

n_robin = getBoundarySet(dh.grid, Γ_robin)
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

# Initialize tangents
global K  = create_sparsity_pattern(dh)
global Kψ = create_sparsity_pattern(dh)
global a  = zeros(dh.ndofs.x)
global d  = zeros(dh.ndofs.x)
global Ψ  = zeros(dh.ndofs.x)
global Fᵢₙₜ  = zeros(dh.ndofs.x)
global rc     = zeros(dh.ndofs.x)
global Fₑₓₜ  = zeros(dh.ndofs.x)
global a      = zeros(dh.ndofs.x)
global Δa     = zeros(dh.ndofs.x)
global res    = zeros(dh.ndofs.x)
global dr_dd  = similar(K)
global ∂rψ_∂d = similar(K)
global ∂g_∂x  = zeros(size(a))
global ∂g_∂u  = zeros(size(d))
global ∂g₂_∂x = zeros(size(a))
global ∂g₂_∂u = zeros(size(d))
global ∂g₃_∂d = zeros(size(d))
global λᵤ = similar(a)
global λψ = similar(a)
global Δ  = -0.025
global nloadsteps = 10
global g  = 0.0
global g₂ = 0.0
global g₃ = 0.0
global α  = 1.0
# # # # # # # # # # # # # # # #
# Init optimization variables #
# # # # # # # # # # # # # # # #
include("initLab.jl")

# ------------------- #
# Boundary conditions #
# ------------------- #
bcdof_bot, _ = setBCY(0.0, dh, n_bot)
bcdof_top, _ = setBCY(Δ, dh, n_top)

bcdof_right, _ = setBCX(0.0, dh, n_sym)
#bcdof_mid, _ = setBCX(0.0, dh, n_mid)
#bcdof_right = Vector{Int64}()
bcdof_left, _ = setBCX(0.0, dh, n_left)

# - - - - - - - - #
# Lås master dofs #
# - - - - - - - - #
bcdof_contact, _ = setBCXY_both(0.0, dh, nₘ) # union(n,n,n) om flera set skall slås samman
# bcdofs_opt = [bcdof_bot; bcdof_top; bcdof_right; bcdof_left];
bcdofs_opt = [bcdof_bot; bcdof_top; bcdof_right; bcdof_left];

ϵᵢⱼₖ      = sortperm(bcdofs_opt)
global bcdofs_opt  = bcdofs_opt[ϵᵢⱼₖ]
global bcval_opt   = bcdofs_opt .* 0.0
global asy_counter = zeros(dh.ndofs.x, 500)

# global low_hist = zeros(length(d), 300)
# global upp_hist = zeros(length(d), 300)
# global d_hist2  = zeros(length(d), 300)

# -------------------- #
# Optimization program #
# -------------------- #
function Optimize(dh)
    global dh0  = deepcopy(dh)
    global λψ   = similar(a)
    global λᵤ   = similar(a)
    global λᵥₒₗ = similar(a)
    Vₘₐₓ = volume(dh, coord, enod) * 1.0 # "volfrac"
    tol   = 1e-3
    global OptIter = 0
    global true_iteration = 0
    global coord₀
    v_hist = zeros(1000)
    g_hist = zeros(1000)
    au_hist = zeros(1000)
    al_hist = zeros(1000)
    global T = zeros(size(a))
    global T[bcdof_bot[iseven.(bcdof_bot)]] .= -1.0
    #global T[bcdof_top[iseven.(bcdof_top)]] .=  1.0
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
        global ∂g_∂d = zeros(dh.ndofs.x,1)
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
        global low
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
        if OptIter % 15 == 0 && true_iteration < 31
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
        global μ = 1e4 #1e3
        global coord₀ = getCoord(getX(dh0), dh0) # x₀
        Ψ, _, Kψ, _, λ = fictitious_solver_with_contact_lab(d, dh0, coord₀, nloadsteps)

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
        global ε = 1e5
        a, _, Fₑₓₜ, Fᵢₙₜ, K = solver_Lab(dh, coord, Δ, nloadsteps)

        # - - - - - - - #
        # Sensitivities #
        # - - - - - - - #
        ∂rᵤ_∂x = similar(K)
        ∂rᵤ_∂x = drᵤ_dx_c(∂rᵤ_∂x, dh, t, a, coord, enod, ε, mp₁, mp₂)
        dr_dd = drψ(dr_dd, dh0, Ψ, λ, d, Γ_robin, coord₀)
        # # # # # # #
        # Objective #
        # # # # # # #
        # g     = -T' * Fᵢₙₜ
        # ∂g_∂x = -T' * ∂rᵤ_∂x
        # ∂g_∂u = -T' * K
        p = 2
        X_ordered = getXfromCoord(coord)
        g     = -contact_pressure(X_ordered, a, ε, p)
        ∂g_∂x = -ForwardDiff.gradient(x -> contact_pressure_ordered(x, a, ε, p), getXinDofOrder(dh, X_ordered, coord))
        ∂g_∂u = -ForwardDiff.gradient(u -> contact_pressure(X_ordered, u, ε, p), a)
        # # # # # # #
        # Adjoints  #
        # # # # # # #
        solveq!(λᵤ, K', ∂g_∂u, bcdofs, bcvals)
        solveq!(λψ, Kψ', ∂g_∂x - ∂rᵤ_∂x' * λᵤ, bcdofs_opt, bcval_opt)
        # # # # # # # # # # #
        # Full sensitivity  #
        # # # # # # # # # # #
        ∂g_∂d =  (-transpose(λψ) * dr_dd)'
        # # # # # # # # # # #
        # Volume constraint #
        # # # # # # # # # # #
        g₁ = volume(dh, coord, enod) / Vₘₐₓ - 1.0
        ∂Ω_∂x = volume_sens(dh, coord)
        solveq!(λᵥₒₗ, Kψ, ∂Ω_∂x, bcdofs_opt, bcval_opt)
        ∂Ω∂d = Real.(-transpose(λᵥₒₗ) * dr_dd ./ Vₘₐₓ)
        # # # # # # # # # # # #
        # Area constraint #
        # # # # # # # # # # # #
        γ_max = 0.10
        γ_min = 0.04
        γc = contact_area(X_ordered, a, ε)

        # g     = -T' * Fᵢₙₜ
        # ∂g_∂x = -T' * ∂rᵤ_∂x
        # ∂g_∂u = -T' * K

        #g₂ = γc / γ_max - 1.0
        Fₘₐₓ  = 100.0
        g₂     = -T' * Fᵢₙₜ ./ Fₘₐₓ - 1.0
        ∂g₂_∂x = -T' * ∂rᵤ_∂x ./ Fₘₐₓ
        ∂g₂_∂u = -T' * K ./ Fₘₐₓ
        solveq!(λᵤ, K', ∂g₂_∂u, bcdofs, bcvals.*0)
        solveq!(λψ, Kψ', ∂g₂_∂x' - ∂rᵤ_∂x' * λᵤ, bcdofs_opt, bcdofs_opt.*0)
        ∂g₂_∂d = Real.((-transpose(λψ) * dr_dd)' )'

        ∂g₃_∂x = ForwardDiff.gradient(x -> contact_area_ordered(x, a, ε), getXinDofOrder(dh, X_ordered, coord))
        ∂g₃_∂u = ForwardDiff.gradient(u -> contact_area(X_ordered, u, ε), a)
        g₃ = 1.0 - γc / γ_min
        solveq!(λᵤ, K',  -∂g₃_∂u./γ_min, bcdofs, bcvals.*0)
        solveq!(λψ, Kψ', -∂g₃_∂x./γ_min - ∂rᵤ_∂x' * λᵤ, bcdofs_opt, bcdofs_opt.*0)
        ∂g₃_∂d = Real.((-transpose(λψ) * dr_dd)' )'
        # # # # #
        # M M A #
        # # # # #
        d_old = d[free_d]
        low_old = low
        upp_old = upp
        d_new, ymma, zmma, lam, xsi, eta, mu, zet, S, low, upp = mmasub(m, n_mma, OptIter, d[free_d], xmin[:], xmax[:],
                                                                        xold1[:], xold2[:], g ./ 1, ∂g_∂d[free_d] ./ 1,
                                                                        vcat(g₁ .* 1e3, g₂*1e3, g₃*1e3),
                                                                        hcat(∂Ω∂d[free_d] .* 1e3,
                                                                        ∂g₂_∂d[free_d]*1e3,
                                                                        ∂g₃_∂d[free_d]*1e3)',
                                                                        low, upp, a0, am, C, d2)
        #d_new, ymma, zmma, lam, xsi, eta, mu, zet, S, low, upp = mmasub(m, n_mma, OptIter, d[free_d], xmin[:], xmax[:], xold1[:], xold2[:], g ./ 1e3, ∂g_∂d[free_d] ./ 1e3, vcat(g₁ .* 1e2, g₂*100), hcat(∂Ω∂d[free_d] .* 1e2, ∂g₂_∂d[free_d]*100)', low, upp, a0, am, C, d2)
        #d_new, ymma, zmma, lam, xsi, eta, mu, zet, S, low, upp = mmasub(m, n_mma, OptIter, d[free_d], xmin[:], xmax[:], xold1[:], xold2[:], g .* 10, ∂g_∂d[free_d] .* 10, g₁ .* 1e2, ∂Ω∂d[free_d]' .* 1e2, low, upp, a0, am, C, d2)
        # ----------------- #
        # Test - new update #
        # ----------------- #
        if true_iteration % 100 == 0
            global α = α / 2
        end
        d_new = d_old   + α .* (d_new - d_old)
        low   = low_old + α .* (low - low_old)
        upp   = upp_old + α .* (upp - upp_old)
        # ----------------- #
        xold2     = xold1
        xold1     = d[free_d]
        d[free_d] = d_new
        change    = norm(d[free_d] .- xold1)

        # # # # # # # # # #
        # Postprocessing  #
        # # # # # # # # # #
        g_hist[true_iteration]  = g
        v_hist[true_iteration]  = g₁
        au_hist[true_iteration] = g₂
        al_hist[true_iteration] = g₃
        kktnorm = change
        println("Iter: ", true_iteration, " Norm of change: ", kktnorm, " Objective: ", g)
        println("Objective: ", g_hist[1:true_iteration])
        println("Volume constraint: ", v_hist[1:true_iteration])
        println("Area constraint", " γ_min ≤ γ ≤ γ_max:  ", γ_min, " ≤ ", γc ," ≤ ", γ_max )
        # ------------ #
        # write to vtu #
        # ------------ #
        postprocess_opt(Ψ, dh0, "results/Current design" * string(true_iteration))
        postprocess_opt(d, dh0, "results/design_variables" * string(true_iteration))
        postprocess_opt(∂g_∂d, dh, "results/🛸" * string(true_iteration))
        # ---- #
        # plot #
        # ---- #
        red_condition_v  = [y > 0 ? :red : :green for y in v_hist[1:true_iteration]]
        red_condition_au = [y > 0 ? :red : :orange for y in au_hist[1:true_iteration]]
        red_condition_al = [y > 0 ? :red : :yellow for y in al_hist[1:true_iteration]]

        p2 = plot(1:true_iteration, [v_hist[1:true_iteration]*10 au_hist[1:true_iteration] al_hist[1:true_iteration]],
                  label=["Volume" "γ_max" "γ_min"],
                  linecolor=hcat(red_condition_v, red_condition_au, red_condition_al),
                  background_color=RGB(0.2, 0.2, 0.2),
                  legend=:outerleft, grid=false)
        hspan!(p2,[-2,0], color = :green, alpha = 0.2, labels = "👌");
        hspan!(p2,[2,0], color = :red, alpha = 0.2, labels = "🤚");
        #p2 = plot(1:true_iteration, v_hist[1:true_iteration]*10 ,
        #           label="Volume" ,
        #           background_color=RGB(0.2, 0.2, 0.2),
        #           legend=:outerleft, grid=false)
        p3 = plot(1:true_iteration, g_hist[1:true_iteration]/1e4, label="Objective",
                  background_color=RGB(0.2, 0.2, 0.2), legend=:outerleft, lc=:purple, grid=false)
        X_c,tract = plotTraction()
        p4 = plot(X_c, tract, label="λ" , marker=4, lc=:tomato, mc=:tomato, grid=false, legend=:outerleft)
        p = plot(p2, p3, p4, layout=(3, 1), size=(600, 600))
        display(p)
        # For investigative purpose
        #low_hist[free_d, true_iteration] = low
        #upp_hist[free_d, true_iteration] = upp
        #d_hist2[free_d, true_iteration]  = d[free_d]
        #@save "asymptoter.jld2" low_hist upp_hist d_hist2
        GC.gc() # Collect garbage
    end
    return g_hist, v_hist, OptIter
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


g_hist, v_hist, OptIter = Optimize(dh)
