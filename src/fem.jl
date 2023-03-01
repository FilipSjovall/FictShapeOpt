#   Module for FEM utilities
#
#   Solveq - Try different solvers
#
function solveq!(a,K,f,bcdof,bcval)
    nd       = size(K,1)
    pdofs    = bcdof
    fdofs    = setdiff(1:nd,pdofs)
    a[fdofs] = K[fdofs,fdofs] \ (f[fdofs] - K[fdofs,pdofs]*bcval)
    a[pdofs] = bcval
end