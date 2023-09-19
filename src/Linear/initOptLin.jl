# Solution vector
ip = Lagrange{2,RefTetrahedron,1}()
qr = QuadratureRule{2,RefTetrahedron}(1)
qr_face = QuadratureRule{1,RefTetrahedron}(1)
cv = CellVectorValues(qr, ip)
fv = FaceVectorValues(qr_face, ip)
# Optimization norm
kktnorm = 1.0
# Boundary values and associated dofs
bcdof,bcval = setBCLin(0.0,dh)

dr_dd = similar(K)
# Material parameters
#mp₀   = [1.0 1.0]
mp₀   = [1.0 5.0]
mp    = [175 80.769230769230759]

t = 1.0
# Optimization parameters
global m             = 2;
global n_mma         = length(d);
global epsimin       = 0.0000001;
global xvalue        = d[:];
global xold1         = xvalue;
global xold2         = xvalue;
global xmin          = -ones(n_mma)/100;
global xmax          =  ones(n_mma)/100;
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

#global xmin[contact_dofs] .= -0.01 # behöver skrivas över
#global xmax[contact_dofs] .=  0.01 # behöver skrivas över
#global xmin[contact_dofs[findall(x -> x % 2 == 0, contact_dofs)]] .= -0.01 # behöver skrivas över
#global xmax[contact_dofs[findall(x -> x % 2 == 0, contact_dofs)]] .=  0.01 # behöver skrivas över

#global xmin[free_d] .= -.05 # behöver skrivas över
#global xmax[free_d] .=  .05 # behöver skrivas över
global xmin[free_d[findall(x -> x % 2 == 0, free_d)]] .= -.05 # behöver skrivas över
global xmax[free_d[findall(x -> x % 2 == 0, free_d)]] .=  .05 # behöver skrivas över

# Gränser sfär / cylinder
#global xmin[register[collect(nₘ),1]] .= -.1
#global xmax[register[collect(nₘ),1]] .=  .1
global xmin[register[collect(nₘ),2]] .= -.5
global xmax[register[collect(nₘ),2]] .=  .5

global low           = xmin;
global upp           = xmax;
#
