using LinearSolve, SparseArrays,  StaticArrays
using FerriteMeshParser,Ferrite, IterativeSolvers, AlgebraicMultigrid, IncompleteLU
##using Gmsh, FerriteGmsh, Plots, Gr # Behöver byggas om....



#load_files()
#filename = "mesh2.txt"
#coord, enod, edof = readAscii(filename);

function solver(dh,coord)
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
    bcdof,bcval = setBC(0.02,dh)
    pdofs       = bcdof
    fdofs       = setdiff(1:ndof,pdofs)
    # ---------- #
    # Set params # // Kanske som input till solver???
    # ---------- #
    mp       = [175 80.769230769230759]
    t        = 1.0

    bcval₀   = bcval

    for n ∈ 1 : 10
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
    #println("Norm coord ", norm(coord))
    return a, dh, Fₑₓₜ, Fᵢₙₜ, K
end


function fictitious_solver(d,dh0,coord₀)
    imax     = 25
    TOL      = 1e-10
    residual = 0.0
    iter     = 1
    global λ
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
    ΔΨ       = zeros(ndof)
    res      = zeros(ndof)
    bcdof,bcval = setBC(0,dh0)

    pdofs       = bcdof
    fdofs       = setdiff(1:ndof,pdofs)
    # ---------- #
    # Set params # // Kanske som input till solver???
    # ---------- #

    bcval₀   = bcval

    for n ∈ 1 : 10
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
            assemGlobal!(Kψ,Fᵢₙₜ,dh0,mp₀,t,Ψ,coord₀,enod,fv,λ,d,ΓN)
            solveq!(ΔΨ, Kψ, -Fᵢₙₜ, bcdof, bcval)
            bcval      = bcval.*0
            res        = Fᵢₙₜ #- Fₑₓₜ
            res[bcdof] = res[bcdof].*0
            residual   = norm(res,2)
            Ψ[bcdof]   = bcval;
            println("Iteration: ", iter, " Residual: ", residual, " λ: ", λ)
        end
        #φ = deepcopy(Ψ)
        #φ[locked_d] .= 0.0
        #φ[free_d] .= 0.1
        #postprocess(φ,dh0)
    end
    return Ψ, dh0, Kψ, Fᵢₙₜ, λ
end
