# Solution vector
kktnorm = 1.0
dr_dd = similar(K) # ?
# Material parameters
# mp₀ = [1.0 1.0]
mp₀= [1. 1.]
mp = [175 80.769230769230759]
#mp = [1.5 0.0003]
#
#mp₀ = mp
#
t = 1.0
# Optimization parameters
global m         = 1;
global n_mma     = length(d);
global epsimin   = 0.0000001;
global xvalue    = d[:];
global xold1     = xvalue;
global xold2     = xvalue;
global C         = 1000 * ones(m);
global d2        = zeros(m);
global a0        = 1;
global am        = zeros(m);
global outeriter = 0;
global kkttol    = 0.001;
global changetol = 0.001;
global outit     = 0;
global change    = 1;
global xmax      =  0.005 * ones(n_mma);
global xmin      = -0.005 * ones(n_mma);
#global xmin[contact_dofs] .= -0.0 # behöver skrivas över
#global xmax[contact_dofs] .=  0.0 # behöver skrivas över
global low = xmin;
global upp = xmax;
