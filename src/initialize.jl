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
    global xmin     = zeros(n);
    global xmax     = ones(n);
    global low      = xmin;
    global upp      = xmax;
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
#end