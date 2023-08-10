#   Module for FEM utilities
#
#   Solveq - Try different solvers
#
function solveq!(x,K,f,bcdofs_in,bcval_in)
    nd          = size(K,1)
    pdofs_in    = bcdofs_in
    fdofs_in    = setdiff(1:nd,pdofs_in)

    prob = LinearProblem(K[fdofs_in,fdofs_in], f[fdofs_in] - K[fdofs_in,pdofs_in]*bcval_in)

    # Solvers
    #    MKLPardisoFactorize(),
    #    MKLPardisoIterate(),
    #    UMFPACKFactorization(),
    #    KLUFactorization())
    #x[fdofs_in] = solve(prob, LUFactorization()).u
    x[fdofs_in] = solve(prob, UMFPACKFactorization()).u

    # Algebraic multigrid
    #ml = ruge_stuben(K[fdofs_in,fdofs_in]) # Construct x Ruge-Stuben solver
    #pl = aspreconditioner(ml)
    #@time x[fdofs_in] = solve(prob, KrylovJL_GMRES(), Pl = pl).u

    # Incomplete LU
    #pl = ilu(K[fdofs_in,fdofs_in], τ = 0.01) # τ needs to be tuned per problem
    #@time x[fdofs_in] = solve(prob, KrylovJL_GMRES(), Pl = pl).u

    # Conjugate gradient
    #IterativeSolvers.cg!(x[fdofs_in], K[fdofs_in,fdofs_in], f[fdofs_in] - K[fdofs_in,pdofs_in]*bcval_in; maxiter=1000)

    # Naive attempt
    x[pdofs_in] = bcval_in

end

function externalForce!(F_ext,r,F_int,F)
    F_ext[freeDofs] = r[freeDofs] - F_int(freeDofs)
    F_ext[dirDofs]  = F[dirDofs]
end
