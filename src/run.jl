using LinearSolve, LinearSolvePardiso, SparseArrays,  StaticArrays
using FerriteMeshParser
using Ferrite
using IterativeSolvers
using AlgebraicMultigrid
using IncompleteLU
##using Gmsh, FerriteGmsh, Plots, Gr # Behöver byggas om....

function load_files()
    include("mesh_reader.jl")

    include("material.jl")

    include("element_routines.jl")

    include("fem.jl")

    include("assemElem.jl")

    include("assem.jl")

    include("sensitivities.jl")
end

load_files()
filename = "mesh2.txt"
coord, enod, edof = readAscii(filename);

function solver(dh)
    imax     = 25
    TOL      = 1e-6
    residual = 0.0
    iter     = 1

    ndof     = size(coord,1)*2 
    nelm     = size(enod,1)

    #K        = zeros(ndof,ndof)
    #K        = spzeros(ndof,ndof)
    
    #sparse_pattern = zeros(ndof,ndof)
    #sparse_pattern[edof,edof] .= 1.0
    #K = sparse(sparse_pattern)
    #sparse_pattern = Nothing
    
    
    # För rätt format kan z-koordinat tas bort i notepad++ med specialsök: ", 0\n" - replace
    #grid = get_ferrite_grid("data/mesh_fine_fine.inp")
    # ----------- #
    # Read "grid" #
    # ----------- #
    #grid = get_ferrite_grid("data/mesh2.inp")
    #dh = DofHandler(grid)
    #add!(dh, :u, 2)
    #close!(dh)
    K  = create_sparsity_pattern(dh)

    #  -------- #
    # Convert   #
    #  -------- #
    # coord <-- dh.grid.nodes 
    # enod  <-- dh.grid.cells

    #  ----- #
    # Init   #
    #  ----- #
    Fᵢₙₜ     = zeros(ndof)
    Fₑₓₜ     = zeros(ndof)
    a        = zeros(ndof)
    Δa       = zeros(ndof)
    res      = zeros(ndof)
    bcdof,bcval = setBC(0.05,dh)
    pdofs       = bcdof
    fdofs       = setdiff(1:ndof,pdofs)
    # ---------- #
    # Set params # // Kanske som input till solver???
    # ---------- #
    mp       = [175 80.769230769230759]
    t        = 1.0

    bcval₀   = bcval

    for n ∈ 1 : 2
        res   = res.*0
        bcval = bcval₀
        residual = 0*residual
        iter  = 0

        fill!(Δa,0.0)
        
        println("Starting equillibrium iteration at loadstep: ",n)

        # # # # # # # # # #
        # Newton solve.  #
        # # # # # # # # # #
        while (iter < imax && residual > TOL ) || iter < 2

            iter += 1
            
            a += Δa

            assemGlobal!(K,Fᵢₙₜ,dh,mp,t,a,coord,enod)
        
            solveq!(Δa, K, -Fᵢₙₜ, bcdof, bcval)

            bcval      = 0*bcval
            res        = Fᵢₙₜ - Fₑₓₜ
            res[bcdof] = 0*res[bcdof]
            residual   = norm(res,2)
            println("Iteration: ", iter, " Residual: ", residual)
        end
    end

    Fₑₓₜ = res - Fᵢₙₜ

    return a, dh, Fₑₓₜ, Fᵢₙₜ, K
end

function fictitious_solver(d,dh0)
    imax     = 25
    TOL      = 1e-6
    residual = 0.0
    iter     = 1

    ndof     = size(coord,1)*2 
    nelm     = size(enod,1)

    # För rätt format kan z-koordinat tas bort i notepad++ med specialsök: ", 0\n" - replace
    #grid = get_ferrite_grid("data/mesh_fine_fine.inp")
    # ----------- #
    # Read "grid" #
    # ----------- #
    #grid = get_ferrite_grid("data/mesh2.inp")
    #dh = DofHandler(grid)
    #add!(dh, :u, 2)
    #close!(dh)
    Kψ  = create_sparsity_pattern(dh0)


    # Borde flyttas utanför och bevaras
   # addfaceset!(dh.grid, "Γ₁", x -> norm(x[1]) ≈ 0.5)
   # addfaceset!(dh.grid, "Γ₂", x -> norm(x[2]) ≈ 0.5)
    #addfaceset!(dh.grid, "Γ₃", x -> norm(x[2]) ≈ 1.0)

  #  ΓN = union(
   #         getfaceset(grid, "Γ₁"),
    #        getfaceset(grid, "Γ₂"),
     #       #getfaceset(grid, "Γ₃"),
      #  )

    ip = Lagrange{2, RefTetrahedron, 2}()
    qr = QuadratureRule{2, RefTetrahedron}(2)
    qr_face = QuadratureRule{1, RefTetrahedron}(2)

    cv = CellVectorValues(qr, ip)
    fv = FaceVectorValues(qr_face, ip)


    #  -------- #
    # Convert   #
    #  -------- #
    # coord <-- dh.grid.nodes 
    # enod  <-- dh.grid.cells

    #  ----- #
    # Init   #
    #  ----- #
    Fᵢₙₜ     = zeros(ndof)
    Fₑₓₜ     = zeros(ndof)
    Ψ        = zeros(ndof)
    #d        = ones(ndof)*0.5
    ΔΨ       = zeros(ndof)
    res      = zeros(ndof)
    bcdof,bcval = setBC(0.0,dh0)
    pdofs       = bcdof
    fdofs       = setdiff(1:ndof,pdofs)
    # ---------- #
    # Set params # // Kanske som input till solver???
    # ---------- #
    mp       = [1.0 1.0] 
    t        = 1.0

    bcval₀   = bcval

    for n ∈ 1 : 2
        res   = res.*0
        bcval = bcval₀
        residual = 0*residual
        iter  = 0
        λ     = 0.1 * n
        fill!(ΔΨ,0.0)
        
        println("Starting equillibrium iteration at loadstep: ",n)

        # # # # # # # # # #
        # Newton solve.  #
        # # # # # # # # # #
        while (iter < imax && residual > TOL ) || iter < 2

            iter += 1
            
            Ψ += ΔΨ

            assemGlobal!(Kψ,Fᵢₙₜ,dh0,mp,t,Ψ,coord,enod,fv,λ,d,ΓN)
        
            solveq!(ΔΨ, Kψ, -Fᵢₙₜ, bcdof, bcval*0)

            bcval      = 0*bcval
            res        = Fᵢₙₜ - Fₑₓₜ
            res[bcdof] = 0*res[bcdof]
            residual   = norm(res,2)
            Ψ[bcdof]   = bcval*0.0;
            println("Iteration: ", iter, "| Residual: ", residual, "| Load level λ: ", λ)
        end
    end
    return Ψ, dh0, Kψ, Fᵢₙₜ
end

function postprocess(a,dh)
        begin
        vtk_grid("hyperelasticity3", dh) do vtkfile
            vtk_point_data(vtkfile, dh, a)
        end
    end
end

### 
#dh0 = dh
#updateCoords!(dh,a)




