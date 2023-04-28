ϵ  = 1e-6
load_files()
Ψ, _, Kψ, Fψ, λ = fictitious_solver(d, dh0, coord₀) # Döp om till "~coord0"
a, _, Fₑₓₜ, Fᵢₙₜ, K = solver(dh,coord)


dofs = Ferrite.celldofs(dh,1)
nods = enod[1][2:end]


x_glob = reshape(coord,(length(dh0.grid.nodes)*2))
ed = a[dofs]
xe = x_glob[dofs]
######
##### Test drᵤ_dx
######
    dX = init_∂X();
    #dr = dr_GP(coord[nods,:],ed,gp,mp,t)
    dr = assem_dr(coord[nods,:],ed,mp,t)
    ke = zeros(6,6)
    fe = zeros(6,2)

    #index1,index2 = [1,2]
    ia = 1
    indexet = nods[ia]
    for pert in 1:2
        coord[indexet,1] = coord[indexet,1] + ϵ * (-real(1*im^(pert)))
        #ke,fe[:,pert] = assemGP(coord[nods,:],ed,gp,mp,t)
        
        ke, fe[:,pert] = assemElem(coord[nods,:],ed,mp,t)
    end
    numsens = (fe[:,2] - fe[:,1])/ϵ
    asens   = dr[:,ia]
    numsens./asens
######
##### Test f_int vs K
    indexet = 3
    ke = zeros(6,6)
    fe = zeros(6,6)
    for pert in 1:2 
        println(" vad är fel ? ")
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


    ke = zeros(6,6)
    fe = zeros(6,2)
    ke2= zeros(6,6)
    fe2 = zeros(6)
    de = 0.5*ones(6)

    cell = CellIterator(dh.grid)
    dofs = Ferrite.celldofs(dh,elnum)
    nods = enod[elnum][2:end]
    cΕll = Ferrite.reinit!(cell.cc,elnum)
    λ    = 0.5
    hej  = CellCache(dh.grid)
    ie   = 103
    Ferrite.reinit!(hej,ie)
    indexet = 2
    for pert in 1:2
        fe[:,pert] .= 0.0
        ke= zeros(6,6)
        if pert == 1
            ed[indexet] = ed[indexet] + ϵ 
        else
            ed[indexet] = ed[indexet] - ϵ 
        end
        # * (-real(1*im^(pert)))
        #ke, fe[:,pert]  = assemElem(coord[nods,:],ed,mp,t)
        #ke,fe[:,pert] = RobinIntegral(ke,fe[:,pert],cΕll,ΓN,fv,ed,λ,de,coord[nods,:])
        for face in 1:nfaces(hej)
            if (cellid(hej), face) in Γ_robin
                face_nods      = [Ferrite.facedof_indices(ip)[face][1]; Ferrite.facedof_indices(ip)[face][2]]
                face_dofs      = [face_nods[1]*2-1; face_nods[1]*2; face_nods[2]*2-1; face_nods[2]*2]
                X              = coord[enod[ie][face_nods.+1] ,:]
                ke[face_dofs,face_dofs],fe[face_dofs,pert]   = Robin(X,ed[face_dofs],de[face_dofs],λ)
            end    
        end
    end
    numsens = (fe[:,1] - fe[:,2])/ϵ
    asens   = ke[:,indexet] 
    numsens./asens
####
## Test r_fictitious w.r.t. d
    ie = 103

    ke = zeros(6,6)
    fe = zeros(6,2)
    dfe= zeros(6,6)
    de = 0.5*ones(6)

    cell = CellIterator(dh.grid)
    dofs = Ferrite.celldofs(dh,ie)
    nods = enod[ie][2:end]
    cΕll = Ferrite.reinit!(cell.cc,ie)
    λ    = 0.5
    hej  = CellCache(dh.grid)
    Ferrite.reinit!(hej,ie)

    λ    = 1

    indexet = 2
    for pert in 1:2
        ke = zeros(6,6)
        dfe= zeros(6,6)
        fe[:,pert] .= 0.0
        if pert == 1
            de[indexet] = de[indexet] + ϵ 
        else
            de[indexet] = de[indexet] - ϵ 
        end
        for face in 1:nfaces(hej)
            if (cellid(hej), face) in Γ_robin
                face_nods      = [Ferrite.facedof_indices(ip)[face][1]; Ferrite.facedof_indices(ip)[face][2]]
                face_dofs      = [face_nods[1]*2-1; face_nods[1]*2; face_nods[2]*2-1; face_nods[2]*2]
                X              = coord[enod[ie][face_nods.+1] ,:]
                dfe[face_dofs,face_dofs],fe[face_dofs,pert]   = Robin(X,ed[face_dofs],de[face_dofs],λ)
            end    
        end
        dfe = -dfe
    end
    numsens = (fe[:,1] - fe[:,2])/ϵ
    asens   = dfe[:,indexet] 
    numsens./asens
#####
## Test objective function and sensitivity
    a, _, Fₑₓₜ, Fᵢₙₜ, K = solver(dh,coord)
    C = zeros(2)
    indexet = 227
    for pert in 1:2
        Fₑₓₜ .=0
        if pert == 1
            a[indexet] = a[indexet] + ϵ 
        else
            a[indexet] = a[indexet] - ϵ 
        end
        
        τ        = [0.1;0.1].*n
        assemGlobal!(K,Fᵢₙₜ,dh,mp,t,a,coord,enod,Γt,τ)
        assemGlobal!(Fₑₓₜ,dh,t,a,coord,enod,Γt,τ)
        #Fₑₓₜ[bcdof] = - Fᵢₙₜ[bcdof]
        C[pert]    = a'*Fₑₓₜ
    end
    numsens = ( C[1] - C[2] ) / ϵ;
    ∂g_∂u = zeros(size(a));
    ∂g_∂u = Fₑₓₜ
    #∂g_∂u[fdofs]  = -a[pdofs]'*K[pdofs,fdofs];
    #∂g_∂u[pdofs]  = -I(length(pdofs))'*Fᵢₙₜ[pdofs]  - (a[pdofs]'*K[pdofs,pdofs])';
    asens = ∂g_∂u[indexet];
    numsens./asens
#####
## Test global dr_dd
    ## bla bla
    load_files()
    Ψ, _, Kψ, _, λ = fictitious_solver(d, dh0, coord₀) 
    λ       = 1.0
    Fψ      = similar(Fᵢₙₜ)
    test    = zeros(284,2)
    indexet = 147
    dr_dd   = drψ(dr_dd,dh0,Ψ,λ,d,Γ_robin,coord₀);
    ϵ       = 1e-6
    for pert in 1:2
        if pert == 1
            d[indexet] = d[indexet] + ϵ 
        else
            d[indexet] = d[indexet] - ϵ 
        end
        coord₀ = getCoord(getX(dh0),dh0)
        assemGlobal!(Kψ,Fψ,dh0,mp₀,t,Ψ,coord₀,enod,λ,d,Γ_robin)
        #assemGlobal!(Kψ,Fᵢₙₜ,dh0,mp₀,t,Ψ,coord₀,enod,λ,d,Γ_robin)
        test[:,pert]                   = Fψ;
    end
    numsens = (test[:,1] - test[:,2])/ϵ
    asens   = dr_dd[:,indexet]
    kvot = numsens./asens
    filter_kvot = filter(x-> abs(x)<10,kvot)
#######
## Test global dr_dx
    Fᵢₙₜ = zeros(ndof)
    dr = similar(K)
    ϵ = 1e-8
    #X = getX(dh)
    incr = zeros(284)
    test = zeros(284,2)
    
    indexet = 20#72
    for pert in 1:2
        if pert == 1
            coord[indexet,1] += ϵ 
            #X[indexet]   +=   ϵ
            #incr[indexet] =   ϵ 
            #updateCoords!(dh,incr)
        else
            coord[indexet,1] -= ϵ 
           # X[indexet]   -=   ϵ
            #incr[indexet] = - ϵ 
            #updateCoords!(dh,incr)
        end
        #coord = getCoord(X,dh) # borde flyttas in i solver..
        #Fᵢₙₜ .=0
        τ        = [0.1;0.1].*10
        assemGlobal!(K,Fᵢₙₜ,dh,mp,t,a,coord,enod,Γt,τ)
        test[:,pert]                   = Fᵢₙₜ;
    end
    
    numsens = ( test[:,1] - test[:,2] )./ϵ;   
    ∂rᵤ_∂x  = drᵤ_dx(∂rᵤ_∂x,dh,mp,t,a,coord,enod,τ, Γt);
    asens   = ∂rᵤ_∂x[:,173];#∂rᵤ_∂x[:,1]; ## dof 6 motsvarar nod 2
    kvot    = numsens./asens;
    filter_kvot = filter(x-> abs(x)<1000,kvot)
#####
## Test dg_dx
    ϵ  = 1e-6
    indexet = 18
    test = zeros(2)
    for pert in 1:2
        if pert == 1
            coord[indexet,2] += ϵ 
            #X[indexet]   +=   ϵ
            #incr[indexet] =   ϵ 
            #updateCoords!(dh,incr)
        else
            coord[indexet,2] -= ϵ 
            #X[indexet]   -=   ϵ
            #incr[indexet] = - ϵ 
            #updateCoords!(dh,incr)
        end
        #coord = getCoord(X,dh) # borde flyttas in i solver..
        Fᵢₙₜ = zeros(ndof)
        Fₑₓₜ .=0
        τ        = [0.1;0.1].*n

        assemGlobal!(K,Fᵢₙₜ,dh,mp,t,a,coord,enod,Γt,τ)
        assemGlobal!(Fₑₓₜ,dh,t,a,coord,enod,Γt,τ)

        #Fₑₓₜ[bcdof] = - Fᵢₙₜ[bcdof]
        test[pert]  =   a'*Fₑₓₜ
    end
    #∂rᵤ_∂x = drᵤ_dx(∂rᵤ_∂x,dh,mp,t,a,coord,enod);
    dFₑₓₜ_dx  = dFext_dx(dFₑₓₜ_dx,dh,mp,t,a,coord,enod,τ,Γt);
    numsens   = (test[1]-test[2])./ϵ   
    dF[fdofs] = a[fdofs]'*dFₑₓₜ_dx[fdofs,fdofs]
    asens     = dF[134]
    #asens    = ∂rᵤ_∂x[:,2]; ## dof 6 motsvarar nod 2
    kvot      = numsens./asens
    


    # Create conversion chart nods < - > dofs
    function index_nod_to_grid(dh,coord)
        X = getX(dh)
        X_nods = reshape_to_nodes(dh, X, :u)[1:2,:] 
        index_register = zeros(Int,length(dh.grid.nodes),2)
        for ii in 1:142
            temp1 = coord[ii,:]
            for jj in 1:142
                temp2 = X_nods[:,jj]
                if temp1 == temp2
                    index_register[ii,:] = [ii,jj]
                end
            end
        end
        return index_register
    end
#####
## Test global sensitivity with adjoint
    load_files()

    #bcdof,bcval = setBCLin(0.0,dh)

    dh0        = deepcopy(dh)
    d          = zeros(284)
    d[free_d] .= 0.1

    locked_d = setdiff(1:ndof,free_d)

    Ψ        = similar(d)
    a        = similar(d)
    Fₑₓₜ     = similar(d)
    K        = create_sparsity_pattern(dh) 
    Kψ       = similar(K)
    dFₑₓₜ_dx = similar(K)
    C        = zeros(2)
    ∂rᵤ_∂x   = similar(K)
    dr_dd    = similar(K)
    ∂rψ_∂d   = similar(K)
    λᵤ       = similar(a)
    λψ       = similar(a)
    fdofs    = setdiff(1:ndof,bcdof)
    indexet  = 264
    ϵ        = 1e-6
    test     = zeros(2)
    for pert in 1:2
        if pert == 1
            # perturbera d
            dh = deepcopy(dh0)
            #dh.grid.nodes = deepcopy(dh0.grid.nodes)
            d[indexet] = d[indexet] + ϵ 
        else
            # perturbera d och resetta dh
            dh = deepcopy(dh0)
            #dh.grid.nodes = deepcopy(dh0.grid.nodes)
            d[indexet] = d[indexet] - ϵ 
        end

        # Check that grid is updated correctly

        Ψ, _, Kψ, _, λ = fictitious_solver(d, dh0, coord₀)
        updateCoords!(dh,Ψ)
        coord = getCoord(getX(dh),dh)
        
        a, _, Fₑₓₜ, Fᵢₙₜ, K             = solver(dh,coord)
        test[pert]                      = a'*Fₑₓₜ;
    end
    if 1 == 1
        ∂g_∂u         = zeros(size(d));

        ∂g_∂u         = Fₑₓₜ;


        dFₑₓₜ_dx      = dFext_dx(dFₑₓₜ_dx,dh,mp,t,a,coord,enod,τ, Γt);

        ∂rᵤ_∂x        = drᵤ_dx(∂rᵤ_∂x,dh,mp,t,a,coord,enod,τ, Γt);

        dr_dd         = drψ(dr_dd,dh0,Ψ,λ,d,Γ_robin,coord₀);
        #dr_dd[locked_d,:] .=0
        #dr_dd[:,locked_d] .=0 
        ∂g_∂x         = zeros(size(d));

        ∂g_∂x[fdofs]  = a[fdofs]'*dFₑₓₜ_dx[fdofs,fdofs];

        solveq!(λᵤ, K',  ∂g_∂u, bcdof, bcval.*0);  # var Fₑₓₜ;

        solveq!(λψ, Kψ', ∂g_∂x - ∂rᵤ_∂x'*λᵤ, bcdof, bcval.*0);

        ∂g_∂d         = -transpose(λψ)*dr_dd;

        asens         =  ∂g_∂d[indexet]

        numsens = (test[1] - test[2])/ϵ

        numsens/asens
    end
    

## Volume constraint direct sensitivity w.r.t. x
    indexet = 14
    test    = zeros(2)
    X       = getX(dh)
    ϵ       = 1e-6
    for pert in 1:2
        if pert == 1
            X[indexet]   +=   ϵ
        else
            X[indexet]   -=   ϵ
        end
        coord      = getCoord(X,dh) # borde flyttas in i solver..
        test[pert] = volume(dh,coord,enod)
    end
    numsens = (test[1] - test[2]) / ϵ
    asens   = volume_sens(dh,coord)

    kvot    = numsens / asens[114]

## Volume constraint - sensitivity w.r.t d via adjoint sensitivity analysis
    load_files()
    indexet = 264
    test    = zeros(2)
    X       = getX(dh)
    ϵ       = 1e-6
    for pert in 1:2
        if pert == 1
            # perturbera d
            dh = deepcopy(dh0)
            #dh.grid.nodes = deepcopy(dh0.grid.nodes)
            d[indexet] = d[indexet] + ϵ 
        else
            # perturbera d och resetta dh
            dh = deepcopy(dh0)
            #dh.grid.nodes = deepcopy(dh0.grid.nodes)
            d[indexet] = d[indexet] - ϵ 
        end

        # Check that grid is updated correctly
        #coord = getCoord(getX(dh0),dh0)
        Ψ, _, Kψ, _, λ               = fictitious_solver(d,dh0,coord₀); # 
        # Update coords
        updateCoords!(dh,Ψ)
        coord = getCoord(getX(dh),dh)

        a, _, _, Fᵢₙₜ, K          = solver(dh,coord);
        test[pert] = volume(dh,coord,enod)
    end
    dr_dd         = drψ(dr_dd,dh0,Ψ,λ,d,Γ_robin,coord₀);

    ∂Ω_∂x         = volume_sens(dh,coord);
    
    λᵥₒₗ          = similar(a);

    solveq!(λᵥₒₗ, Kψ,∂Ω_∂x, bcdof, bcval.*0);

    ∂Ω∂d   = -transpose(λᵥₒₗ)*dr_dd;

    numsens = (test[1] - test[2])/ϵ
    asens   = ∂Ω∂d[indexet]

    numsens/asens