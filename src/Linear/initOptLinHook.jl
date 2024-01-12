# Solution vector
kktnorm = 1.0
dr_dd = similar(K) # ?
# Material parameters
mp  = [175 80.769230769230759]
mp₀ = [1. 5.].*100 #* 3 # .*500
#mp₀ = [1. 1.]
#mp₀ = mp
#
t = 1.0
# Optimization parameters
global m         = 1;
#global n_mma     = length(d); # length(free_d) ?
global n_mma = length(free_d);
global epsimin   = 0.0000001;
#global xvalue    = d[:];
global xvalue = d[free_d];
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
global xmax      =  1.0*ones(n_mma);
global xmin      = -1.0*ones(n_mma);

m_indices = findall(x -> x in collect(nₘ), free_d)
s_indices = findall(x -> x in collect(nₛ), free_d)


global xmax[m_indices].=   .5
global xmin[m_indices].=  -.5
global xmax[s_indices].=   .5
global xmin[s_indices].=  -.5

global low = xmin;
global upp = xmax;
