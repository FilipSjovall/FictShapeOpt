# Solution vector
a = zeros(size(coord,1)*2)
# Optimization norm
kktnorm = 1.0
# Boundary values and associated dofs
bcdof,bcval = setBCLin(0.0,dh)
# Sensitivities
∂g_∂x = zeros(size(a)) 
∂g_∂u = zeros(size(d)) 
∂rᵤ_∂x= similar(K) 
dr_dd = similar(K) 
# Material parameters
mp₀   = [1.0 1.0]
mp    = [175 80.769230769230759]
# Optimization parameters
global m             = 1;
global n             = length(d);
global epsimin       = 0.0000001;
global xval          = d[:];
global xold1         = xval;
global xold2         = xval;
global xmin          = -ones(n)/20;
global xmax          = ones(n)/20;
global C             = 1000*ones(m);
global d2            = zeros(m);
global a0            = 1;
global am            = zeros(m);
global outeriter     = 0;
global kkttol        = 0.001;
global changetol     = 0.001;
global kktnorm       = kkttol + 10;
global outit         = 0;
global change        = 1;
global xmin[free_d] .= -0.8
global xmax[free_d] .=  0.8
global low           = xmin;
global upp           = xmax;