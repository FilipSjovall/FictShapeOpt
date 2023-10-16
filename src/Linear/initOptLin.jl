# Optimization norm
kktnorm = 1.0
# Boundary values and associated dofs
bcdof,bcval = setBCLin(0.0,dh)

dr_dd = similar(K)
# Material parameters
#mp₀   = [1.0 1.0]
mp₀   = [0.0 5.0]
mp    = [175 80.769230769230759]

t = 1.0
# Optimization parameters
#global m             = 2;
global m             = 1;
global n_mma         = length(d);
global epsimin       = 0.0000001;
global xvalue        = d[:];
global xold1         = xvalue;
global xold2         = xvalue;
global xmin          = zeros(n_mma)#-ones(n_mma)/1000;
global xmax          = zeros(n_mma)# ones(n_mma)/1000;
#global xmin          = zeros(n_mma);
#global xmax          = zeros(n_mma);
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

#global xmin[free_d] .= -.0 # behöver skrivas över
#global xmax[free_d] .=  .0 # behöver skrivas över
#global xmin[free_d[findall(x -> x % 2 == 0, free_d)]] .= -.5 # behöver skrivas över
#global xmax[free_d[findall(x -> x % 2 == 0, free_d)]] .=  .5 # behöver skrivas över

# Gränser sfär / cylinder
#global xmin[register[collect(nₘ),1]] .= -.0
#global xmax[register[collect(nₘ),1]] .=  .0

global xmax .=  0.2
global xmin .= -0.2

#global xmin[1:2:end-1] .= -0.01
#global xmax[1:2:end-1] .=  0.01



# global xmax[register[collect(nₘ), 1]].=   .5
# global xmin[register[collect(nₘ), 1]].=  -.5
# global xmax[register[collect(nₘ), 2]] .=  .5
# global xmin[register[collect(nₘ), 2]] .= -.5

global low           = -ones(n_mma);
global upp           =  ones(n_mma);
#
