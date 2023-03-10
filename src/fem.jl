#   Module for FEM utilities
#
#   Solveq - Try different solvers
#
function solveq!(a,K,f,bcdof,bcval)
    nd       = size(K,1)
    pdofs    = bcdof
    fdofs    = setdiff(1:nd,pdofs)
    
    
    prob = LinearProblem(K[fdofs,fdofs], f[fdofs] - K[fdofs,pdofs]*bcval)
    
    # Solvers
    #    MKLPardisoFactorize(),
    #    MKLPardisoIterate(),
    #    UMFPACKFactorization(),
    #    KLUFactorization())
    a[fdofs] = solve(prob, MKLPardisoIterate()).u
    
    # Algebraic multigrid
    #ml = ruge_stuben(K[fdofs,fdofs]) # Construct a Ruge-Stuben solver
    #pl = aspreconditioner(ml)
    #@time a[fdofs] = solve(prob, KrylovJL_GMRES(), Pl = pl).u
    
    # Incomplete LU
    #pl = ilu(K[fdofs,fdofs], τ = 0.01) # τ needs to be tuned per problem
    #@time a[fdofs] = solve(prob, KrylovJL_GMRES(), Pl = pl).u

    # Conjugate gradient
    #IterativeSolvers.cg!(a[fdofs], K[fdofs,fdofs], f[fdofs] - K[fdofs,pdofs]*bcval; maxiter=1000)

    # Naive attempt
    a[pdofs] = bcval
    
end