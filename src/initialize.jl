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
#end


#function init_mats()
    #global d  = zeros(dh.ndofs.x)
    global d  = zeros(dh.ndofs.x)
    global Ψ  = similar(d)
    global a  = similar(d)
    global Fₑₓₜ= similar(d)
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
        for ii in 1:449
            
            temp2 = X_nods[:,ii]
            for jj in 1:449
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
    
    for cell in CellIterator(dh)
        for face in 1:nfaces(cell)
            if (cellid(cell), face) in ΓN
                face_nod = Ferrite.faces(dh.grid.cells[cellid(cell)])[1]
                println(face_nod)
                # Detta är fel dofs
                doftemp  = [register[face_nod[1],2]*2-1 register[face_nod[1],2]*2 ]
                append!(free_d,doftemp)
                doftemp  = [register[face_nod[2],2]*2-1 register[face_nod[2],2]*2 ]
                append!(free_d,doftemp)
                #append!(free_d,[face_nod[1]*2-1 face_nod[1]*2 face_nod[2]*2-1 face_nod[2]*2])
            end
        end
    end

    xmin[free_d] .= -0.25
    xmax[free_d] .=  0.25

    d[free_d]    .= 0.25

    global low      = xmin;
    global upp      = xmax;
    #end