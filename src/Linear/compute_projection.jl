using Symbolics
using LinearAlgebra
@variables N1 N2 xs1[1:2] xs2[1:2] xm[1:2] ξ n₁[1:2] n₂[1:2]

N1 = 1 - ξ
N2 = ξ

expr1 = [N1.*xs1[1] + N2.*xs2[1] - xm[1];N1.*xs1[2] + N2.*xs2[2] - xm[2];0]

expr2 = [N1 .* n₁[1] + N2 .* n₂[1]; N1 .* n₁[2] + N2 .* n₂[2]; 0]

expr3 = cross(expr1,expr2)

s = simplify(expr3)
#Symbolics.solve_for(expr3,ξ)

s[3]