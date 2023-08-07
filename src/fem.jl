#   Module for FEM utilities
#
#   Solveq - Try different solvers
#
function solveq!(x,K,f,bcdof,bcval)
    nd       = size(K,1)
    pdofs    = bcdof
    fdofs    = setdiff(1:nd,pdofs)

    prob = LinearProblem(K[fdofs,fdofs], f[fdofs] - K[fdofs,pdofs]*bcval)

    # Solvers
    #    MKLPardisoFactorize(),
    #    MKLPardisoIterate(),
    #    UMFPACKFactorization(),
    #    KLUFactorization())
    #x[fdofs] = solve(prob, LUFactorization()).u
    x[fdofs] = solve(prob, UMFPACKFactorization()).u

    # Algebraic multigrid
    #ml = ruge_stuben(K[fdofs,fdofs]) # Construct x Ruge-Stuben solver
    #pl = aspreconditioner(ml)
    #@time x[fdofs] = solve(prob, KrylovJL_GMRES(), Pl = pl).u

    # Incomplete LU
    #pl = ilu(K[fdofs,fdofs], τ = 0.01) # τ needs to be tuned per problem
    #@time x[fdofs] = solve(prob, KrylovJL_GMRES(), Pl = pl).u

    # Conjugate gradient
    #IterativeSolvers.cg!(x[fdofs], K[fdofs,fdofs], f[fdofs] - K[fdofs,pdofs]*bcval; maxiter=1000)

    # Naive attempt
    x[pdofs] = bcval

end

function externalForce!(F_ext,r,F_int,F)
    F_ext[freeDofs] = r[freeDofs] - F_int(freeDofs)
    F_ext[dirDofs]  = F[dirDofs]
end
