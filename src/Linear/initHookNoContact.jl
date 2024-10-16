# Solution vector
kktnorm = 1.0
dr_dd = similar(K) # ?
# Material parameters
mp  = [175 80.769230769230759]
mp₀ = [0.001 1.0] # " [K₀ G₀] "
#
t = 1.0
# Optimization parameters
global m = 1;
global n_mma = length(free_d);
global epsimin = 0.0000001;
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
global xmax =  .5 * ones(n_mma);
global xmin = -.5 * ones(n_mma);

global low = xmin;
global upp = xmax;
