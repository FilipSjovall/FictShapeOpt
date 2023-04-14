ϵ  = 1e-8
load_files()
a,dh = solver();
mp = [1.0 1.0]
t  = 1.0
gp = 1

dofs = Ferrite.celldofs(dh,2)
nods = enod[2][2:end]


x_glob = reshape(coord,(length(dh0.grid.nodes)*2))
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
    #dr = dr_GP(coord[nods,:],ed,gp,mp,t)
    dr = assem_dr(coord[nods,:],ed,mp,t)

    
    ke = zeros(12,12)
    fe = zeros(12,2)

    index1,index2 = [1,2]
    for pert in 1:2
        coord[nods[index1],index2] = coord[nods[index1],index2] + ϵ * (-real(1*im^(pert)))
        #ke,fe[:,pert] = assemGP(coord[nods,:],ed,gp,mp,t)
        
        ke, fe[:,pert] = assemElem(coord[nods,:],ed,mp,t)
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
        if pert == 1
            ed[indexet] = ed[indexet] + ϵ 
        else
            ed[indexet] = ed[indexet] - ϵ 
        end
        ke, fe[:,pert] = assemElem(coord[nods,:],ed,mp,t)
    end
    numsens = (fe[:,2] - fe[:,1])/ϵ
    asens   = ke[:,indexet]
    numsens./asens
#####
## Test f_int vs K in fictious domain
    elnum = 92


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

    indexet = 2
    for pert in 1:2
        fe[:,pert] .= 0.0
        ke= zeros(12,12)
        if pert == 1
            ed[indexet] = ed[indexet] + ϵ 
        else
            ed[indexet] = ed[indexet] - ϵ 
        end
        # * (-real(1*im^(pert)))
        #ke, fe[:,pert]  = assemElem(coord[nods,:],ed,mp,t)
        #ke,fe[:,pert] = RobinIntegral(ke,fe[:,pert],cΕll,ΓN,fv,ed,λ,de,coord[nods,:])
        ke[f1dofs,f1dofs],fe[f1dofs,pert]   = Robin(coord[nods[f1],:],ed[f1dofs],de[f1dofs],λ)
    end
    numsens = (fe[:,1] - fe[:,2])/ϵ
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
    λ    = 1

    indexet = 2
    for pert in 1:2
        ke = zeros(12,12)
        dfe= zeros(12,12)
        fe[:,pert] .= 0.0
        if pert == 1
            de[indexet] = de[indexet] + ϵ 
        else
            de[indexet] = de[indexet] - ϵ 
        end
        dfe[f1dofs,f1dofs],fe[f1dofs,pert]   = Robin(coord[nods[f1],:],ed[f1dofs],de[f1dofs],λ)
        dfe = -dfe
    end
    numsens = (fe[:,1] - fe[:,2])/ϵ
    asens   = dfe[:,indexet] 
    numsens./asens
#####
## Test objective function and sensitivity
    a, dh, Fₑₓₜ, Fᵢₙₜ, K = solver(dh);
    C = zeros(2)
    indexet = 868
    for pert in 1:2
        if pert == 1
            a[indexet] = a[indexet] + ϵ 
        else
            a[indexet] = a[indexet] - ϵ 
        end
        assemGlobal!(K,Fᵢₙₜ,dh,mp,t,a,coord,enod)
        C[pert] = -a[pdofs]'*Fᵢₙₜ[pdofs]
    end
    numsens = ( C[1] - C[2] ) / ϵ;
    ∂g_∂u = zeros(size(a));
    ∂g_∂u[fdofs]  = -a[pdofs]'*K[pdofs,fdofs];
    ∂g_∂u[pdofs]  = -I(length(pdofs))'*Fᵢₙₜ[pdofs]  - (a[pdofs]'*K[pdofs,pdofs])';
    asens = ∂g_∂u[indexet];
    numsens./asens
#####
## Test global dr_dd
    ## bla bla
    Fψ = similar(Fᵢₙₜ)
    C = zeros(898,2)
    indexet = 294
    dr_dd = drψ(dr_dd,dh0,Ψ,fv,λ,d,ΓN);
    ϵ = 1e-6
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
    numsens = (C[:,1] - C[:,2])/ϵ
    asens   = dr_dd[:,indexet]
    kvot = numsens./asens
    filter_kvot = filter(x-> abs(x)<10,kvot)
## Test global dr_dx
    ϵ = 1e-6
    X = getX(dh)
    incr = zeros(898)
    C = zeros(898,2)
    indexet = 6
    for pert in 1:2
        if pert == 1
            X[indexet]   +=   ϵ
            #incr[indexet] =   ϵ 
            #updateCoords!(dh,incr)
        else
            X[indexet]   -=   ϵ
            #incr[indexet] = - ϵ 
            #updateCoords!(dh,incr)
        end
        coord = getCoord(X,dh) # borde flyttas in i solver..
        Fᵢₙₜ .=0
        assemGlobal!(K,Fᵢₙₜ,dh,mp,t,a,coord,enod)
        #a, _, Fₑₓₜ, Fᵢₙₜ, K       = solver(dh);
        #C[pert]                   = compliance(Fₑₓₜ,a)
        C[:,pert]                   = Fᵢₙₜ;
        #C[pert] = -a[pdofs]'*Fᵢₙₜ[pdofs]
    end
    ∂rᵤ_∂x = drᵤ_dx(∂rᵤ_∂x,dh,mp,t,a,coord,enod);
    numsens = (C[:,1]-C[:,2])./ϵ   
    asens   = ∂rᵤ_∂x[:,indexet]
    asens   = ∂rᵤ_∂x[:,2]; ## dof 6 motsvarar nod 2
    kvot    = numsens./asens
    filter_kvot = filter(x-> abs(x)<10,kvot)


## Test dg_dx
    indexet = 6
    C = zeros(2)
    for pert in 1:2
        if pert == 1
            X[indexet]   +=   ϵ
            #incr[indexet] =   ϵ 
            #updateCoords!(dh,incr)
        else
            X[indexet]   -=   ϵ
            #incr[indexet] = - ϵ 
            #updateCoords!(dh,incr)
        end
        coord = getCoord(X,dh) # borde flyttas in i solver..
        Fᵢₙₜ .=0
        assemGlobal!(K,Fᵢₙₜ,dh,mp,t,a,coord,enod)
        #a, _, Fₑₓₜ, Fᵢₙₜ, K       = solver(dh);
        #C[pert]                   = compliance(Fₑₓₜ,a)
        #C[:,pert]                   = Fᵢₙₜ;
        C[pert] = -a[pdofs]'*Fᵢₙₜ[pdofs]
    end
    ∂rᵤ_∂x = drᵤ_dx(∂rᵤ_∂x,dh,mp,t,a,coord,enod);
    numsens = (C[1]-C[2])./ϵ   
    asens   = -a[pdofs]'*∂rᵤ_∂x[pdofs,2]
    #asens   = ∂rᵤ_∂x[:,2]; ## dof 6 motsvarar nod 2
    kvot    = numsens./asens
    


    # Create conversion chart nods < - > dofs
    function index_nod_to_grid(dh,coord)
        X = getX(dh)
        X_nods = reshape_to_nodes(dh, X, :u)[1:2,:] 
        index_register = zeros(Int,length(dh.grid.nodes),2)
        for ii in 1:449
            temp1 = coord[ii,:]
            for jj in 1:449
                temp2 = X_nods[:,jj]
                if temp1 == temp2
                    index_register[ii,:] = [ii,jj]
                end
            end
        end
        return index_register
    end

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
    λᵤ = similar(a)
    λψ = similar(a)

    mp     = [175 80.769230769230759]
    t      = 1.0 

    indexet = 465
    ϵ       = 1e-6


    test = zeros(2)
    for pert in 1:2
        if pert == 1
            # perturbera d
            #dh = dh0
            dh.grid.nodes = deepcopy(dh0.grid.nodes)
            d[indexet] = d[indexet] + ϵ 
        else
            # perturbera d och resetta dh
            #dh = dh0
            dh.grid.nodes = deepcopy(dh0.grid.nodes)
            d[indexet] = d[indexet] - ϵ 
        end

        # Check that grid is updated correctly
        coord = getCoord(getX(dh0),dh0)
        Ψ, _, Kψ, _, λ               = fictitious_solver(d,dh0); # 
        # Update coords
        updateCoords!(dh,Ψ)
        coord = getCoord(getX(dh),dh)

        a, _, _, Fᵢₙₜ, K          = solver(dh);
        test[pert]                   = -a[pdofs]'*Fᵢₙₜ[pdofs];
    end

    ∂g_∂u = zeros(size(d))
    ∂g_∂u[fdofs]  = -a[pdofs]'*K[pdofs,fdofs];
    #∂g_∂u[pdofs]  = -I(length(pdofs))'*Fᵢₙₜ[pdofs]  - (a[pdofs]'*K[pdofs,pdofs])';

    ∂rᵤ_∂x = drᵤ_dx(∂rᵤ_∂x,dh,mp,t,a,coord,enod);

    dr_dd = drψ(dr_dd,dh0,Ψ,fv,λ,d,ΓN);

    ∂g_∂x = zeros(size(d))
    ∂g_∂x[fdofs]  = -a[pdofs]'*∂rᵤ_∂x[pdofs,fdofs]
    #∂g_∂x[pdofs]  = -a[pdofs]'*∂rᵤ_∂x[pdofs,pdofs]

    solveq!(λᵤ, K, ∂g_∂u, bcdof, bcval.*0);  # var Fₑₓₜ;
    solveq!(λψ, Kψ,∂g_∂x - ∂rᵤ_∂x'*λᵤ, bcdof, bcval.*0);

    ∂g_∂d   = -transpose(λψ)*dr_dd;

    numsens = (test[1] - test[2])/ϵ
    asens   = ∂g_∂d[indexet]

    numsens/asens


    ## Volume constraint direct sensitivity w.r.t. x
    indexet = 6
    test    = zeros(2)
    X       = getX(dh)
    ϵ       = 1e-6
    for pert in 1:2
        if pert == 1
            X[indexet]   +=   ϵ
        else
            X[indexet]   -=   ϵ
        end
        coord = getCoord(X,dh) # borde flyttas in i solver..
        test[pert] = volume(dh)
    end
    numsens = (test[1] - test[2]) / ϵ
    asens   = volume_sens(dh,coord)

    kvot    = numsens / asens[2]

    ## Volume constraint - sensitivity w.r.t d via adjoint sensitivity analysis
    indexet = 296
    test    = zeros(2)
    X       = getX(dh)
    ϵ       = 1e-6
    for pert in 1:2
        if pert == 1
            # perturbera d
            #dh = dh0
            dh.grid.nodes = deepcopy(dh0.grid.nodes)
            d[indexet] = d[indexet] + ϵ 
        else
            # perturbera d och resetta dh
            #dh = dh0
            dh.grid.nodes = deepcopy(dh0.grid.nodes)
            d[indexet] = d[indexet] - ϵ 
        end

        # Check that grid is updated correctly
        coord = getCoord(getX(dh0),dh0)
        Ψ, _, Kψ, _, λ               = fictitious_solver(d,dh0,coord); # 
        # Update coords
        updateCoords!(dh,Ψ)
        coord = getCoord(getX(dh),dh)

        a, _, _, Fᵢₙₜ, K          = solver(dh,coord);
        test[pert] = volume(dh)
    end

    dr_dd = drψ(dr_dd,dh0,Ψ,fv,λ,d,ΓN);

    ∂Ω_∂x = volume_sens(dh,coord)
    
    λᵥₒₗ    = similar(a)

    solveq!(λᵥₒₗ, Kψ,∂Ω_∂x, bcdof, bcval.*0);

    ∂Ω∂d   = -transpose(λψ)*dr_dd;

    numsens = (test[1] - test[2])/ϵ
    asens   = ∂Ω∂d[indexet]

    numsens/asens