# Solution vector
kktnorm = 1.0
dr_dd = similar(K) #
# - - - - - - - - - - - #
# Material parameters   #
# N/m² = N/(10⁶ mm² )   #
# - - - - - - - - - - - #
# https://www.azom.com/properties.aspx?ArticleID=920
mp₁ = [180 80].*1e3     # [K G]
mp₂ = [2.5 0.1].*1e3    #
#mp₂ = [154 1.5]    #
#mp₀ = [0.01 5.0]        #
mp₀ = [1.0 5.0]        #
# - - - - - - - - - - - #
#
t = 1.0
# Optimization parameters
global m = 3;
global n_mma   = length(free_d);
global epsimin = 0.0000001;
#global xvalue    = d[:];
global xvalue = d[free_d];
global xold1  = xvalue;
global xold2  = xvalue;
global C  = 1000 * ones(m);
global d2 = zeros(m);
global a0 = 1;
global am = zeros(m);
global outeriter = 0;
global kkttol    = 0.001;
global changetol = 0.001;
global outit  = 0;
global change = 1;
global xmax =  0.1 * ones(n_mma) ;
global xmin = -0.1 * ones(n_mma) ;

# global xmax[1:2:end-1] .= 0.0001
# global xmin[1:2:end-1] .= 0.0001

global low = xmin;
global upp = xmax;
