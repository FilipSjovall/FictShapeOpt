from sympy import *
init_printing()
from sympy.vector import CoordSys3D
N = CoordSys3D('N')
# Define the symbolic variables and the system of nonlinear equations
N1, N2 = symbols('N1 N2')
xs1, xs2, xm = symbols('xs1 xs2 xm', cls=Function)
n1, n2 = symbols('n1 n2', cls=Function)
xi = Symbol('xi')

N1 = 1 - xi
N2 = xi

expr1 = Eq(N1*xs1(xi) + N2*xs2(xi) - xm(xi))
expr2 = Eq(N1*n1(xi) + N2*n2(xi))
expr3 = Eq(expr1.rhs.cross(expr2.rhs), [0,0,0])

# Solve the system of equations
soln = solve(expr3, xi)

# Print the solution
print(soln)