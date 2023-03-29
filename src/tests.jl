ϵ  = 1e-6
load_files()
a,dh = solver();
mp = [1.0 1.0]
t  = 1.0
gp = 1

dofs = Ferrite.celldofs(dh,1)
nods = enod[1][2:end]


x_glob = reshape(coord,(449*2))
ed = a[dofs]
xe = x_glob[dofs]
addfaceset!(dh.grid, "Γ₁", x -> norm(x[1]) ≈ 0.5)
addfaceset!(dh.grid, "Γ₂", x -> norm(x[2]) ≈ 0.5)
#addfaceset!(dh.grid, "Γ₃", x -> norm(x[2]) ≈ 1.0)

ΓN = union(
    getfaceset(dh.grid, "Γ₁"),
    getfaceset(dh.grid, "Γ₂"),
    #getfaceset(grid, "Γ₃"),
)

ip = Lagrange{2, RefTetrahedron, 2}()
qr = QuadratureRule{2, RefTetrahedron}(2)
qr_face = QuadratureRule{1, RefTetrahedron}(2)

cv = CellVectorValues(qr, ip)
fv = FaceVectorValues(qr_face, ip)
######
##### Test drᵤ_dx
######
    dX = init_∂X();
    dr = dr_GP(coord[nods,:],ed,gp,mp,t)
    index1,index2 = [1,2]

    ke = zeros(12,12)
    fe = zeros(12,2)

    for pert in 1:2
        coord[nods[index1],index2] = coord[nods[index1],index2] + ϵ * (-real(1*im^(pert)))
        println((-real(1*im^(pert))))
        ke,fe[:,pert] = assemGP(coord[nods,:],ed,gp,mp,t)
    end
    numsens = (fe[:,2] - fe[:,1])/ϵ
    asens   = dr[:,2]
    numsens./asens
######
##### Test f_int vs K
    indexet = 7

    ke = zeros(12,12)
    fe = zeros(12,2)
    for pert in 1:2
        ed[indexet] = ed[indexet] + ϵ * (-real(1*im^(pert)))
        ke, fe[:,pert] = assemElem(coord[nods,:],ed,mp,t)
    end
    numsens = (fe[:,2] - fe[:,1])/ϵ
    asens   = ke[:,indexet]
    numsens./asens
#####
## Test f_int vs K in fictious domain
    elnum = 35

    ke = zeros(12,12)
    fe = zeros(12,2)
    ke2= zeros(12,12)
    fe2 = zeros(12)
    de = 0.5*ones(12)

    cell = CellIterator(dh.grid)
    dofs = Ferrite.celldofs(dh,elnum)
    nods = enod[elnum][2:end]
    cΕll    = Ferrite.reinit!(cell.cc,elnum)
    λ    = 0.5

    indexet = 8
    for pert in 1:2
        fe[:,pert] .= 0.0
        ke= zeros(12,12)
        ed[indexet] = ed[indexet] + ϵ * (-real(1*im^(pert)))
        #ke, fe[:,pert]  = assemElem(coord[nods,:],ed,mp,t)
        ke,fe[:,pert] = RobinIntegral(ke,fe[:,pert],cΕll,ΓN,fv,ed,λ,de)
    end
    numsens = (fe[:,2] - fe[:,1])/ϵ
    asens   = ke[:,indexet] 
    numsens./asens
####
## Test r_fictitious w.r.t. d
    elnum = 35

    ke = zeros(12,12)
    fe = zeros(12,2)
    dfe= zeros(12,12)
    de = 0.5*ones(12)

    cell = CellIterator(dh.grid)
    dofs = Ferrite.celldofs(dh,elnum)
    nods = enod[elnum][2:end]
    cΕll = Ferrite.reinit!(cell.cc,elnum)
    λ    = 0.5

    indexet = 8
    for pert in 1:2
        ke = zeros(12,12)
        dfe= zeros(12,12)
        fe[:,pert] .= 0.0
        de[indexet] = de[indexet] + ϵ * (-real(1*im^(pert)))
        ke,fe[:,pert] = RobinIntegral(ke,fe[:,pert],cΕll,ΓN,fv,ed,λ,de)
        dfe = d_RobinIntegral(dfe,cΕll,ΓN,fv,ed,λ,de)
    end
    numsens = (fe[:,2] - fe[:,1])/ϵ
    asens   = dfe[:,indexet] 
    numsens./asens
#####
## Test objective function and sensitivity
    a, dh, Fₑₓₜ, Fᵢₙₜ, K = solver(dh);
    C = zeros(2)
    indexet = 13
    for pert in 1:2
        if pert == 1
            a[indexet] = a[indexet] + ϵ 
        else
            a[indexet] = a[indexet] - ϵ 
        end
        C[pert] = compliance(Fₑₓₜ,a)
    end
    numsens = ( C[1] - C[2] ) / ϵ
    asens   = Fₑₓₜ[indexet]
    numsens./asens
#####
## Test global dr_dd
    ## bla bla
    Fψ = similar(Fᵢₙₜ)
    C = zeros(898,2)
    indexet = 294
    dr_dd = drψ(dr_dd,dh0,Ψ,fv,λ,d,ΓN);
    ϵ = 0.1
    for pert in 1:2
        if pert == 1
            d[indexet] = d[indexet] + ϵ 
        else
            d[indexet] = d[indexet] - ϵ 
        end
        coord = getCoord(getX(dh0),dh0)
        assemGlobal!(Kψ,Fψ,dh0,mp,t,Ψ,coord,enod,fv,λ,d,ΓN)
        C[:,pert]                   = Fψ;
    end
    numsens = (C[:,1]-C[:,2])/ϵ
    asens   = dr_dd[:,indexet]
    kvot = numsens./asens
    filter_kvot = filter(x-> abs(x)<10,kvot)
## Test global dr_dx
    ϵ = 1e-2
    
    X = getX(dh)
    incr = zeros(898)
    C = zeros(898,2)
    for pert in 1:2
        if pert == 1
            X[indexet]   +=   ϵ
            incr[indexet] =   ϵ 
            updateCoords!(dh,incr)
        else
            X[indexet]   -=   ϵ
            incr[indexet] = - ϵ 
            updateCoords!(dh,incr)
        end
        coord = getCoord(X,dh) # borde flyttas in i solver..
        assemGlobal!(K,Fᵢₙₜ,dh,mp,t,a,coord,enod)
        #a, _, Fₑₓₜ, Fᵢₙₜ, K       = solver(dh);
        #C[pert]                   = compliance(Fₑₓₜ,a)
        C[:,pert]                   = Fᵢₙₜ;
    end
    ∂rᵤ_∂x = drᵤ_dx(∂rᵤ_∂x,dh,mp,t,a,coord,enod);
    numsens = (C[:,1]-C[:,2])/ϵ
    asens   = ∂rᵤ_∂x[:,indexet]
    kvot    = numsens./asens
    filter_kvot = filter(x-> abs(x)<10,kvot)
## Test global sensitivity with adjoint
    load_files()
    filename = "mesh2.txt"
    coord, enod, edof = readAscii(filename);
    grid = get_ferrite_grid("data/mesh2.inp")
    dh = DofHandler(grid)
    add!(dh, :u, 2)
    close!(dh)
    addfaceset!(dh.grid, "Γ₁", x -> norm(x[1]) ≈ 0.5)
    addfaceset!(dh.grid, "Γ₂", x -> norm(x[2]) ≈ 0.5)
    ΓN = union(
            getfaceset(grid, "Γ₁"),
            getfaceset(grid, "Γ₂"),
        )
    ip = Lagrange{2, RefTetrahedron, 2}()
    qr = QuadratureRule{2, RefTetrahedron}(2)
    qr_face = QuadratureRule{1, RefTetrahedron}(2)
    cv = CellVectorValues(qr, ip)
    fv = FaceVectorValues(qr_face, ip)

    bcdof,bcval = setBC(0.0,dh)

    dh0 = deepcopy(dh)
    d  = ones(898)*0.5

    Ψ  = similar(d)
    a  = similar(d)
    Fₑₓₜ= similar(d)
    K  = create_sparsity_pattern(dh) 
    Kψ = similar(K)
    C  = zeros(2)
    ∂rᵤ_∂x = similar(K)
    dr_dd = similar(K)
    ∂rψ_∂d = similar(K)
    
    mp     = [175 80.769230769230759]
    t      = 1.0 

    indexet = 294
    ϵ       = 1e-6

    
    l  = similar(a)
    l .= 0.5

    for pert in 1:2
        if pert == 1
            # perturbera d
            #dh = dh0
            dh.grid.nodes = deepcopy(dh0.grid.nodes)
            d[indexet] = d[indexet] + ϵ 
        else
            # perturbera d och resetta dh
            dh.grid.nodes = deepcopy(dh0.grid.nodes)
            #dh = dh0
            d[indexet] = d[indexet] - ϵ 
        end

        # Check that grid is updated correctly
        #println(dh.grid.nodes[1],dh0.grid.nodes[1])
        coord = getCoord(getX(dh0),dh0)
        Ψ, _, Kψ,Fψ               = fictitious_solver(d,dh0);
        # update coords
        updateCoords!(dh,Ψ)
        
        # Check that grid is updated correctly
        #println(dh.grid.nodes[1],dh0.grid.nodes[1])
        
        coord = getCoord(getX(dh),dh)
        a, _, Fₑₓₜ, Fᵢₙₜ, K       = solver(dh);
        #C[pert]                   = compliance(Fₑₓₜ,a)
        C[pert]                   = transpose(l)*a;
    end

    ∂g_∂u  = Fₑₓₜ
    
    λ      = 0.2

    dX = init_∂X();
    ∂rᵤ_∂x = drᵤ_dx(∂rᵤ_∂x,dh,mp,t,a,coord,enod);
    dr_dd = drψ(dr_dd,dh0,Ψ,fv,λ,d,ΓN);
    #∂rψ_∂d = drψ(dr_dd,dh0,Ψ,fv,λ,d,ΓN);
    λψ = similar(a);
    λᵤ = similar(a);
    solveq!(λᵤ, K', l, bcdof, bcval*0);  # var Fₑₓₜ;
    solveq!(λψ, Kψ', -transpose(λᵤ)*∂rᵤ_∂x, bcdof, bcval*0);
    ∂g_∂d   = -transpose(λψ)*dr_dd;
    numsens = (C[1]-C[2])/ϵ
    asens   = ∂g_∂d[indexet]
    numsens/asens
####
    # 
    # call adjoint function