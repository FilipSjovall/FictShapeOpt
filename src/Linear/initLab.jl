# Solution vector
kktnorm = 1.0
dr_dd = similar(K) # ?
# Material parameters
# N/m² = N/(10⁶ mm² )
mp₁ = [185.19 75.76].*1e3
mp₂ = [5.556 0.34]
mp₀ = [0.01 2.0]
#
t = 1.0
# Optimization parameters
global m = 2;
global n_mma = length(free_d);
global epsimin = 0.0000001;
#global xvalue    = d[:];
global xvalue = d[free_d];
global xold1 = xvalue;
global xold2 = xvalue;
global C = 1000 * ones(m);
global d2 = zeros(m);
global a0 = 1;
global am = zeros(m);
global outeriter = 0;
global kkttol = 0.001;
global changetol = 0.001;
global outit = 0;
global change = 1;
global xmax =  0.01 * ones(n_mma);
global xmin = -0.01 * ones(n_mma);

global low = xmin;
global upp = xmax;
