using LinearSolve, LinearSolvePardiso, SparseArrays, StaticArrays
using FerriteMeshParser, Ferrite, IterativeSolvers, AlgebraicMultigrid, IncompleteLU

include("initialize.jl")
include("mma.jl")
init_hyper()

dh0 = deepcopy(dh);


function Optimize(dh)
    dh0 = deepcopy(dh)
    #init_MMA();
    # Flytta allt nedan till init_opt?
    λψ   = similar(a)
    λᵤ   = similar(a)
    λᵥₒₗ = similar(a)
    Vₘₐₓ = 0.6
    l    = similar(a)
    l   .= 0.5
    tol  = 1e-4
    OptIter = 0
    coord₀  = coord
    #
    while kktnorm > tol || OptIter < 3 #&& OptIter < 5

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
        global low
        global upp
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


        OptIter += 1
        global g_ini
        global g_hist = []
        global v_hist = []
        #A() = A(Float64[],[]).
        global pdofs = bcdof
        global fdofs = setdiff(1:length(a), pdofs)
        # # # # # # # # # # # # # #
        # Fictitious equillibrium #
        # # # # # # # # # # # # # #
        coord = getCoord(getX(dh0), dh0) # x₀ 
        Ψ, _, Kψ, _, λ = fictitious_solver(d, dh0, coord₀) # Döp om till "~coord0"

        # # # # # # 
        # Filter  #
        # # # # # # 
        dh    = deepcopy(dh0) 
        updateCoords!(dh, Ψ) # x₀ + Ψ = x
        coord = getCoord(getX(dh), dh)

        # **TEST**
        #updateCoords!(dh0,Ψ); # x₀ + Ψ = x

        # # # # # # # # #
        # Equillibrium  #
        # # # # # # # # #
        a, _, _, Fᵢₙₜ, K = solver(dh,coord)

        # # # # # # # 
        # Objective #
        # # # # # # #
        #g       =  compliance(l,a);
        g = -a[pdofs]' * Fᵢₙₜ[pdofs]

        # # # # # # # # # 
        # Sensitivities #
        # # # # # # # # # 
        ∂g_∂x = zeros(size(a))
        ∂g_∂u = zeros(size(d))


        ∂g_∂u[fdofs] = -a[pdofs]' * K[pdofs, fdofs]
        #println("norm coord outside ", norm(coord) )
        ∂rᵤ_∂x = drᵤ_dx(∂rᵤ_∂x, dh, mp, t, a, coord, enod)
        dr_dd = drψ(dr_dd, dh0, Ψ, fv, λ, d, ΓN)
        #∂g_∂d   = -transpose(λψ)*dr_dd; # gör till funktion?


        # # # # # # #
        # Adjoints  #
        # # # # # # #
        solveq!(λᵤ, K', ∂g_∂u, bcdof, bcval * 0)  # var Fₑₓₜ;
        ∂g_∂x[fdofs] = a[pdofs]' * ∂rᵤ_∂x[pdofs, fdofs]
        solveq!(λψ, Kψ', ∂g_∂x - ∂rᵤ_∂x' * λᵤ, bcdof, bcval * 0)

        # # # # # # # # # # #
        # Full sensitivity  #
        # # # # # # # # # # #
        ∂g_∂d = (-transpose(λψ) * dr_dd)'


        #
        # Volume constraint
        #
        g₁    = volume(dh) / Vₘₐₓ - 1
        ∂Ω_∂x = volume_sens(dh,coord)
        solveq!(λᵥₒₗ, Kψ, ∂Ω_∂x, bcdof, bcval.*0);
        ∂Ω∂d  = -transpose(λᵥₒₗ)*dr_dd ./ Vₘₐₓ; 

        # # # # 
        # MMA # 
        # # # #  
        X, ymma, zmma, lam, xsi, eta, mu, zet, S, low, upp = mmasub(m, n, OptIter, d, xmin, xmax, xold1, xold2, 10 * g, 10 * ∂g_∂d, g₁, ∂Ω∂d', low, upp, a0, am, C, d2)
        xold2 = xold1
        xold1 = d
        d = X
        #display(d[free_d])
        #Xg=lower_bound+(upper_bound-lower_bound).*X;
        change = norm(d .- xold1)
        if OptIter == 1
            g_ini = g
        end

        append!(v_hist, g₁)
        append!(g_hist, g / g_ini)
        #The residual vector of the KKT conditions is calculated:
        #residu,kktnorm,residumax = kktcheck(m,n,X,ymma,zmma,lam,xsi,eta,mu,zet,S, xmin,xmax,∂g_∂d,[0.0],zeros(size(d)),a0,a,C,d2);
        kktnorm = change
        println("Iter: ", OptIter, " Norm of change: ", kktnorm, " Objective: ", g)
        coord = getCoord(getX(dh0), dh0)
        postprocess_opt(Ψ, dh0, "results/Shape" * string(OptIter))
        coord = getCoord(getX(dh), dh)
        postprocess_opt(a, dh, "results/Deformation" * string(OptIter))
        println("Objective: ", g_hist, " Constraint: ", g₁)
    end
    return g_hist, v_hist
end


