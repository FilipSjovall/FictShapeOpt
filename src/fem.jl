#   Module for FEM utilities
#
#   Solveq - Try different solvers
#
function solveq!(a,K,f,bcdof,bcval)
    nd       = size(K,1)
    pdofs    = bcdof
    fdofs    = setdiff(1:nd,pdofs)
    A = K[fdofs,fdofs] 
    b = f[fdofs] - K[fdofs,pdofs]*bcval
    
    
    prob = LinearProblem(A, b)
    
    # Samling lösare
    #for alg in (
    #    MKLPardisoFactorize(),
    #    MKLPardisoIterate(),
    #    UMFPACKFactorization(),
    #    KLUFactorization())
    #
    #    @time a[fdofs] = solve(prob, alg).u
    #end

    a[fdofs] = solve(prob, MKLPardisoFactorize()).u

    # Algebraic multigrid
    #ml = ruge_stuben(A) # Construct a Ruge-Stuben solver
    #pl = aspreconditioner(ml)
    #@time a[fdofs] = solve(prob, KrylovJL_GMRES(), Pl = pl).u
    
    # Incomplete LU
    
    #pl = ilu(A, τ = 0.01) # τ needs to be tuned per problem
    #@time a[fdofs] = solve(prob, KrylovJL_GMRES(), Pl = pl).u

    #a[fdofs] =  A\b


    # Naive attempt
    #a[fdofs] = K[fdofs,fdofs] \ (f[fdofs] - K[fdofs,pdofs]*bcval)
    # Attempt at using Pardiso through LinearSolve
    #a[fdofs] = solve(prob, MKLPardisoFactorize(); cache_kwargs...).u
    # Attempt at solve using Pardiso - needs linking?
    #solve!(ps,a,K[fdofs,fdofs],(f[fdofs] - K[fdofs,pdofs]*bcval) )
    
    
    a[pdofs] = bcval
end