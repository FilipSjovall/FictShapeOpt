using LinearSolve, LinearSolvePardiso, SparseArrays, 
      #StaticArrays, CairoMakie,
      IterativeSolvers, AlgebraicMultigrid, IncompleteLU    

#ENV["PATH"]




function load_files()
      include("..//mesh_reader.jl")

      include("..//material.jl")
  
      include("..//fem.jl")
  
      include("assemElemLin.jl")
  
      include("assemLin.jl")
  
      include("sensitivitiesLin.jl")

      include("run_linear.jl")

      include("..//mma.jl")

      include("initLin.jl")

      include("initOptLin.jl")
end

load_files()


function Optimize(dh)
    
    # Flytta allt nedan till init_opt?
        dh0 = deepcopy(dh)
        λψ   = similar(a)
        λᵤ   = similar(a)
        λᵥₒₗ = similar(a)
        Vₘₐₓ = 1.2
        l    = similar(a)
        l   .= 0.5
        tol  = 1e-6
        OptIter = 0
        coord₀  = coord
        global g_hist = zeros(200)
        global v_hist = zeros(200)
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
            global pdofs = bcdof
            global fdofs = setdiff(1:length(a), pdofs)
            global locked_d = setdiff(1:length(a),free_d)
            global low
            global upp
        # # # # # # # # # # # # # #
        OptIter += 1
        
        # # # # # # # # # # # # # #
        # Fictitious equillibrium #
        # # # # # # # # # # # # # #
        coord₀ = getCoord(getX(dh0), dh0) # x₀ 
        Ψ, _, Kψ, _, λ = fictitious_solver(d, dh0, coord₀) # Döp om till "~coord0"

        # # # # # # 
        # Filter  #
        # # # # # # 
        dh    = deepcopy(dh0) 
        updateCoords!(dh, Ψ) # x₀ + Ψ = x
        coord = getCoord(getX(dh), dh)

        # # # # # # # # #
        # Equillibrium  #
        # # # # # # # # #
        a, _, Fₑₓₜ, _, K = solver(dh,coord)
        
        # # # # # # # 
        # Objective #
        # # # # # # #
        #g = -a[pdofs]' * Fᵢₙₜ[pdofs]
        g = a' * Fₑₓₜ

        # # # # # # # # # 
        # Sensitivities #
        # # # # # # # # # 
        # ∂g_∂u[fdofs] = -a[pdofs]' * K[pdofs, fdofs]
        ∂g_∂u        = Fₑₓₜ
        ∂rᵤ_∂x       = drᵤ_dx(∂rᵤ_∂x, dh, mp, t, a, coord, enod)
        dr_dd        = drψ(dr_dd, dh0, Ψ, fv, λ, d, Γ_robin)

        # # # # # # #
        # Adjoints  #
        # # # # # # #
        solveq!(λᵤ, K', ∂g_∂u, bcdof, bcval * 0)  # var Fₑₓₜ;
        ∂g_∂x[fdofs] = a[pdofs]' * ∂rᵤ_∂x[pdofs, fdofs] # ???
        solveq!(λψ, Kψ', ∂g_∂x - ∂rᵤ_∂x' * λᵤ, bcdof, bcval * 0)

        # # # # # # # # # # #
        # Full sensitivity  #
        # # # # # # # # # # #
        ∂g_∂d = (-transpose(λψ) * dr_dd)'
        ∂g_∂d[locked_d] .= 0.0 # fulfix?

        # # # # # # # # # # #
        # Volume constraint #
        # # # # # # # # # # #
        g₁    = volume(dh,coord) / Vₘₐₓ - 1
        ∂Ω_∂x = volume_sens(dh,coord)
        solveq!(λᵥₒₗ, Kψ, ∂Ω_∂x, bcdof, bcval.*0);
        ∂Ω∂d  = -transpose(λᵥₒₗ)*dr_dd ./ Vₘₐₓ; 
        ∂Ω∂d[locked_d] .= 0.0

        # # # # #
        # M M A # 
        # # # # # 
        X, ymma, zmma, lam, xsi, eta, mu, zet, S, low, upp = mmasub(m, n, OptIter, d, xmin, xmax, xold1, xold2, 10 * g, 10 * ∂g_∂d, g₁, ∂Ω∂d', low, upp, a0, am, C, d2)
        xold2 = xold1
        xold1 = d
        d = X
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
        if mod(OptIter,5) == 0
            coord = getCoord(getX(dh0), dh0)
            postprocess_opt(Ψ, dh0, "results/Shape" * string(OptIter))
            coord = getCoord(getX(dh), dh)
            postprocess_opt(a, dh, "results/Deformation" * string(OptIter))
        end 
        println("Objective: ", g_hist, " Constraint: ", g₁)
        if OptIter == 50
            break
        end
    end

    # # # # # # # # #
    # Plot history  #
    # # # # # # # # #
    fig = Figure()
    ax1, l1 = lines(fig[1, 1], 1..OptIter, g_hist[1:OptIter], color = :red)
    ax2, l2 = lines(fig[2, 1], 1..OptIter, v_hist[1:OptIter], color = :blue)
    Legend(fig[1:2, 2], [l1, l2], ["Objective", "Constraint"])
    fig
    return g_hist, v_hist, OptIter
end


using Makie
OptIter=37
fig = Figure()
ax1, l1 = lines(fig[1, 1], 1..OptIter, g_hist[1:OptIter], color = :red)
ax2, l2 = lines(fig[2, 1], 1..OptIter, v_hist[1:OptIter], color = :blue)
Legend(fig[1:2, 2], [l1, l2], ["Objective", "Constraint"])
fig
