function init_hyper()
    include("run.jl")
    load_files()

    #init_grid()

    #init_mats()

    # Choose objective and constraints
end


#function init_grid()
    load_files()
    filename = "mesh2.txt"
    coord, enod, edof = readAscii(filename);
    grid = get_ferrite_grid("data/Quadratic/mesh2.inp")
    dh = DofHandler(grid)
    add!(dh, :u, 2)
    close!(dh)
    addfaceset!(dh.grid, "Γ₁", x -> norm(x[1]) ≈ 0.0)
    addfaceset!(dh.grid, "Γ₂", x -> norm(x[2]) ≈ 0.0)
    ΓN = union(
            getfaceset(grid, "Γ₁"),
            getfaceset(grid, "Γ₂"),
        )
    Γ1 = getfaceset(grid, "Γ₁")
    Γ2 = getfaceset(grid, "Γ₂")
    ip = Lagrange{2, RefTetrahedron, 2}()
    qr = QuadratureRule{2, RefTetrahedron}(2)
    qr_face = QuadratureRule{1, RefTetrahedron}(2)
    cv = CellVectorValues(qr, ip)
    fv = FaceVectorValues(qr_face, ip)
    bcdof,bcval = setBC(0.0,dh)
#end


#function init_mats()
    #global d  = zeros(dh.ndofs.x)
    global d      = zeros(dh.ndofs.x)
    global Ψ      = similar(d)
    global a      = similar(d)
    global Fₑₓₜ   = similar(d)
    global K      = create_sparsity_pattern(dh) 
    global Kψ     = similar(K)
    global ∂rᵤ_∂x = similar(K)
    global dr_dd  = similar(K)
    global ∂rψ_∂d = similar(K)
    global mp     = [175 80.769230769230759]
    global mp₀    = [1.0 1.0]
    global t      = 1.0 
#end

#function init_MMA()
    global m        = 1;
    global n        = length(d);
    global epsimin  = 0.0000001;
    global xval     = d[:];
    global xold1    = xval;
    global xold2    = xval;
    global xmin     = -ones(n)/20;
    global xmax     = ones(n)/20;
    #global xmin     = zeros(n);
    #global xmax     = zeros(n);
    
    global C        = 1000*ones(m);
    global d2       = zeros(m);
    global a0       = 1;
    global am       = zeros(m);
    global outeriter= 0;
    global kkttol   = 0.001;
    global changetol= 0.001;
    global kktnorm  = kkttol + 10;
    global outit    = 0;
    global change   = 1;


    function index_nod_to_grid(dh,coord)
        coord = getCoord(getX(dh),dh);
        X = getX(dh)
        X_nods = reshape_to_nodes(dh, X, :u)[1:2,:] 
        index_register = zeros(Int,length(dh.grid.nodes),2)
        for ii in 1:length(coord[:,1])
            temp2 = X_nods[:,ii]
            for jj in 1:length(coord[:,1])
                temp1 = coord[jj,:]
                if temp1 == temp2
                    index_register[ii,:] = [ii,jj]
                end
            end
        end
        return index_register
    end

    register = index_nod_to_grid(dh,coord)

    free_d = []
    
    addnodeset!(dh.grid, "Γ₃", x -> norm(x[1]) == 0.5)
    addnodeset!(dh.grid, "Γ₄", x -> norm(x[2]) == 0.5)

    nodx = Ferrite.getnodeset(dh.grid,"Γ₃")
    nody = Ferrite.getnodeset(dh.grid,"Γ₄")

    for inod in nodx
       append!(free_d,register[inod,2]*2-1)
    end

    for jnod in nody
       append!(free_d,register[jnod,2]*2)
    end
    #Γ4 = getnodeset(grid, "Γ₄")
    global xmin[free_d] .= -0.8
    global xmax[free_d] .=  0.8
    d[free_d]    .=  0.0

    global low    = xmin;
    global upp    = xmax;
    #end

    addfaceset!(dh.grid, "Γ₃", x -> norm(x[1]) == 0.5)
    addfaceset!(dh.grid, "Γ₄", x -> norm(x[2]) == 0.5)
    Γt = union(
            getfaceset(grid, "Γ₃"),
            getfaceset(grid, "Γ₄"),
        )
    #Γt = getfaceset(grid,"Γₜ")

function postprocess_opt(Ψ, dh, str)
    begin
        vtk_grid(str, dh) do vtkfile
            vtk_point_data(vtkfile, dh, Ψ)
        end
    end
end

function postprocess(a,dh)
        begin
        vtk_grid("hyperelasticity_2", dh) do vtkfile
            vtk_point_data(vtkfile, dh, a)
        end
    end
end