ϵ  = 1e-6

a,dh = solver();
mp = [1.0 1.0]
t  = 1.0
gp = 1
#dofs = edof[1,:]
dofs = Ferrite.celldofs(dh,1)
nods = enod[1][2:end]


x_glob = reshape(coord,(449*2))
ed = a[dofs]
xe = x_glob[dofs]

## Test drᵤ_dx
dX = init_∂X();
load_files()
@benchmark dr = dr_GP(coord[nods,:],ed,gp,mp,t)
index1,index2 = [1,1]

ke = zeros(12,12)
fe = zeros(12,2)

for pert in 1:2
    coord[nods[index1],index2] = coord[nods[index1],index2] + ϵ * (-real(1*im^(pert)))
    println((-real(1*im^(pert))))
    ke,fe[:,pert] = assemGP(coord[nods,:],ed,gp,mp,t)
end

numsens = (fe[:,2] - fe[:,1])/ϵ
asens   = dr[:,1]
numsens./asens



## Test f_int vs K

## Test f_int vs K in fictious domain

## Test objective function and sensitivity

## Test global sensitivity with adjoint
