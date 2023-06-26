# Solution vector
K    = create_sparsity_pattern(dh)
# Optimization norm
kktnorm = 1.0
# Boundary values and associated dofs
bcdof,bcval = setBCLin(0.0,dh)

dr_dd = similar(K)
# Material parameters
mp₀   = [1.0 1.0]
mp    = [175 80.769230769230759]
# Optimization parameters
global m             = 1;
global n_mma         = length(d);
global epsimin       = 0.0000001;
global xval          = d[:];
global xold1         = xval;
global xold2         = xval;
global xmin          = -ones(n_mma)/20;
global xmax          =  ones(n_mma)/20;
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
global xmin[contact_dofs] .= -0.05
global xmax[contact_dofs] .=  0.05
global xmin[contact_dofs[findall(x -> x % 2 == 0, contact_dofs)]] .= -0.2
global xmax[contact_dofs[findall(x -> x % 2 == 0, contact_dofs)]] .=  0.2
global low           = xmin;
global upp           = xmax;
