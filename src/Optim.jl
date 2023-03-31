using LinearSolve, LinearSolvePardiso, SparseArrays,  StaticArrays
using FerriteMeshParse,Ferrite, IterativeSolvers, AlgebraicMultigrid, IncompleteLU

include("initialize.jl")
include("mma.jl")
init_hyper()


function Optimize(dh)
    dh0 = deepcopy(dh);
    init_MMA();
    # Flytta allt nedan till init_opt?
    λψ = similar(a);
    λᵤ = similar(a);
    l  = similar(a)
    l .= 0.5
    tol= 1e-4
    OptIter = 0
    #
    while kktnorm > tol || OptIter < 1
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
        OptIter +=1
        # # # # # # # # # # # # # #
        # Fictitious equillibrium #
        # # # # # # # # # # # # # #
        coord = getCoord(getX(dh0),dh0); # x₀ 
        Ψ, _, Kψ,_               = fictitious_solver(d,dh0); # Döp om till "~coord0"
        # # # # # # 
        # Filter  #
        # # # # # # 
        updateCoords!(dh,Ψ); # x₀ + Ψ = x
        coord = getCoord(getX(dh),dh);
        # # # # # # # # #
        # Equillibrium  #
        # # # # # # # # #
        a, _, _, _, K       = solver(dh);
        # # # # # # # # # 
        # Sensitivities #
        # # # # # # # # # 
        ∂rᵤ_∂x = drᵤ_dx(∂rᵤ_∂x,dh,mp,t,a,coord,enod);
        dr_dd  = drψ(dr_dd,dh0,Ψ,fv,λ,d,ΓN);
        # # # # # # # 
        # Objective #
        # # # # # # #
        g       =  compliance(l,a);
        ∂g_∂d   = -transpose(λψ)*dr_dd; # gör till funktion?
        ∂g_∂x   =  zeros(size(a));
        ∂g_∂u   =  l # gör till funktion ? 
        # # # # # # #
        # Adjoints  #
        # # # # # # #
        solveq!(λᵤ, K', l, bcdof, bcval*0);  # var Fₑₓₜ;
        solveq!(λψ, Kψ', -transpose(λᵤ)*∂rᵤ_∂x, bcdof, bcval*0);
        # # # # # # # # # # #
        # Full sensitivity  #
        # # # # # # # # # # #
        ∂g_∂d   = -transpose(λψ)*dr_dd;
        # # # # 
        # MMA # 
        # # # #  
        X,ymma,zmma,lam,xsi,eta,mu,zet,S,low,upp=mmasub(m,n,OptIter,d,xmin,xmax,xold1,xold2, g,∂g_∂d,0,zeros(size(d)),low,upp,a0,am,C,d2);
        xold2 = xold1;
        xold1 = d;
        d  = X;
        Xg=lower_bound+(upper_bound-lower_bound).*X;
        change=norm(d .-xold1);
        #The residual vector of the KKT conditions is calculated:
        residu,kktnorm,residumax = kktcheck(m,n,X,ymma,zmma,lam,xsi,eta,mu,zet,S, xmin,xmax,df0dx,fval,dfdx[:],a0,a,C,d2);
        println("Iter: ",OptIter, "KKT-norm: ",kktnorm, "Objective: ", f0val)
    end
end 