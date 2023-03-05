using LinearSolve, LinearSolvePardiso, SparseArrays,  StaticArrays
using FerriteMeshParser
using Ferrite
using IterativeSolvers
using AlgebraicMultigrid
using IncompleteLU

function load_files()
    include("mesh_reader.jl")

    include("material.jl")

    include("element_routines.jl")

    include("fem.jl")
end

load_files()
filename = "mesh_fine_fine.txt"

coord, enod, edof = readAscii(filename);



function solver()
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
    grid = get_ferrite_grid("data/mesh_fine_fine.inp")
    dh = DofHandler(grid)
    add!(dh, :u, 2)
    close!(dh)
    K = create_sparsity_pattern(dh)

    Fᵢₙₜ     = zeros(ndof)
    Fₑₓₜ     = zeros(ndof)
    a        = zeros(ndof)
    Δa       = zeros(ndof)
    res      = zeros(ndof)

    F        = zeros(3,2,2)
    ef       = zeros(4,3)
    D        = zeros(3,3,3)
    Sgp      = zeros(4,3)
    es       = zeros(3,3)
    kₑ       = zeros(12,12)
    fₑ       = zeros(12)


    bcdof,bcval = setBC(0.005,dh)
    pdofs       = bcdof
    fdofs       = setdiff(1:ndof,pdofs)
    # bcval    = zeros(2)
    # bcval    = [0.0; 0.005]
    # bcdof    = Vector{Int64}(undef,2)
    #bcdof    = [1;7] - fine
    # bcdof = [1;4]

    #bcval    = zeros(6)
    #bcval    = [0.0; 0.0; 0.05; 0.05; 0.0; 0.0]
    #bcdof    = Vector{Int64}(undef,6)
    #bcdof    = [3; 4; 7; 9; 5; 6]
    
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

        while (iter < imax && residual > TOL ) || iter < 2

            iter += 1
            assembler = start_assemble(K,Fᵢₙₜ)
            a += Δa
            #@time for ie ∈ 1 : nelm
            ie = 0
            for cell in CellIterator(dh)
                ie +=1                
                fill!(kₑ, 0.0)
                fill!(fₑ, 0.0)
                cell_dofs= celldofs(cell)
                # lägg ihop alla c2tl6-funktioner 
                c2tl6_d!( F, a[cell_dofs], coord[enod[ie][2:7],:] )

                for gp = 1:3
                    ef[:,gp] = [F[gp,1,1] F[gp,1,2] F[gp,2,1] F[gp,2,2]]
                end
                
                ## D - dneohooke
                dneohookeD!(D,ef,mp)
                ## S - neohooke
                neohookeS!(Sgp,ef,mp)

                for gp = 1:3 
                    es[1,gp] = Sgp[1,gp] # Sₑ som namn istället?
                    es[2,gp] = Sgp[2,gp]
                    es[3,gp] = Sgp[4,gp]
                end 
                c2tl6_e!(kₑ, coord[enod[ie][2:7],:], t, D,  a[cell_dofs], es)
                c2tl6_f!(fₑ, coord[enod[ie][2:7],:], t, es, a[cell_dofs])

                
                assemble!(assembler, cell_dofs, kₑ,fₑ)

                # println(kₑ)
                #@inbounds Fᵢₙₜ[cell_dofs] += fₑ            
            end

            # println("F_int ",Fᵢₙₜ[fdofs])
            # println("F2 ", F2[fdofs])
            
        
            res = Fᵢₙₜ - Fₑₓₜ
            # prob = LinearProblem(K[fdofs,fdofs], -res[fdofs] - K[fdofs,pdofs]*bcval)
            # Δa[fdofs] = solve(prob, MKLPardisoIterate()).u
            solveq!(Δa, K, -Fᵢₙₜ, bcdof, bcval)

            #IterativeSolvers.cg!(Δa[fdofs], K[fdofs,fdofs], (-res[fdofs] - K[fdofs,pdofs]*bcval); maxiter=1000)
            #Δa[pdofs] = bcval
            # println("HL ", K2[fdofs,pdofs]*bcval)
            # println("HL ", K[fdofs,pdofs]*bcval)
            # if iter == 2
            #     println("a ", a)
            #     println("Δa ", Δa)
            #     println("bcdof ", bcdof)
            #     println("pdof ", pdofs)
            #     iter = imax
            #     spy(K)
            # end
            # a += Δa
            # Fᵢₙₜ  = Fᵢₙₜ .*0

            # for ie ∈ 1:nelm
            #     fₑ       = fₑ.*0

            #     c2tl6_d!( F, a[edof[ie,:]], coord[enod[ie][2:7],:] )
            #     for gp = 1:3
            #         ef[:,gp] = [F[gp,1,1] F[gp,1,2] F[gp,2,1] F[gp,2,2]]
            #     end
                
            #     ## S - neohooke
            #     neohookeS!(Sgp,ef,mp)
            #     for gp = 1:3 
            #         es[1,gp] = Sgp[1,gp] # Sₑ som namn istället?
            #         es[2,gp] = Sgp[2,gp]
            #         es[3,gp] = Sgp[4,gp]
            #     end 

            #     c2tl6_f!(fₑ, coord[enod[ie][2:7],:], t, es, a[edof[ie,:]])
            #     Fᵢₙₜ[edof[ie,:]]           += fₑ            
            # end

            bcval      = 0*bcval
            res        = Fᵢₙₜ - Fₑₓₜ
            res[bcdof] = 0*res[bcdof]
            residual   = norm(res,2)
            println("Iteration: ", iter, " Residual: ", residual)
        end
    end
end