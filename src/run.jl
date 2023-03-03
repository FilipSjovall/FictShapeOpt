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
filename = "mesh.txt"

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
    grid = get_ferrite_grid("data/mesh.inp")
    dh = DofHandler(grid)
    add!(dh, :u, 2)
    close!(dh)
    K = create_sparsity_pattern(dh)


    K2       = spzeros(ndof,ndof)
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

    bcval    = zeros(2)
    bcval    = [0.0; 0.005]
    bcdof    = Vector{Int64}(undef,2)
    #bcdof    = [1;7] - fine
    bcdof = [1;4]

    #bcval    = zeros(6)
    #bcval    = [0.0; 0.0; 0.05; 0.05; 0.0; 0.0]
    #bcdof    = Vector{Int64}(undef,6)
    #bcdof    = [3; 4; 7; 9; 5; 6]
    
    mp       = [175 80.769230769230759]
    t        = 1.0
    bcval₀   = bcval
    for n ∈ 1 : 1#10
        res   = res.*0
        bcval = bcval₀
        residual = 0*residual
        iter  = 0

        println("Starting equillibrium iteration at loadstep: ",n)

        while (iter < imax && residual > TOL ) || iter < 2
            iter += 1
            iter  = imax
            #K     = K.*0
            assembler = start_assemble(K,Fᵢₙₜ)
            #Fᵢₙₜ  = Fᵢₙₜ .*0

            #@time for ie ∈ 1 : nelm
            ie = 0
            for cell in CellIterator(dh)
                ie +=1                
                kₑ       = kₑ.*0
                fₑ       = fₑ.*0
                cell_dofs= celldofs(cell)
                # långsam?
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

                #println("fem-sparse")
                #FEMSparse.assemble_local_matrix!(assembler, edof[ie,:], kₑ)
                K2[edof[ie,:],edof[ie,:]]  += kₑ 
                #assemble!(assembler, edof[ie,:], kₑ)
                println("enod ", enod[ie,:])
                println("edof ", edof[ie,:])
                println("celldofs ", cell_dofs)
                
                assemble!(assembler, cell_dofs, kₑ,fₑ)
                #@inbounds Fᵢₙₜ[cell_dofs] += fₑ            
            end


            nd       = size(K,1)
            pdofs    = bcdof
            fdofs    = setdiff(1:nd,pdofs)
            println("kdiff ",norm(K[fdofs,fdofs])-norm(K2[fdofs,fdofs]))
            #res = Fᵢₙₜ - Fₑₓₜ
            #@time solveq!(Δa, K, -res, bcdof, bcval)
            solveq!(Δa, K, -Fᵢₙₜ, bcdof, bcval)
       
            a += Δa
            Fᵢₙₜ  = Fᵢₙₜ .*0

            for ie ∈ 1:nelm
                fₑ       = fₑ.*0

                c2tl6_d!( F, a[edof[ie,:]], coord[enod[ie][2:7],:] )
                for gp = 1:3
                    ef[:,gp] = [F[gp,1,1] F[gp,1,2] F[gp,2,1] F[gp,2,2]]
                end
                
                ## S - neohooke
                neohookeS!(Sgp,ef,mp)
                for gp = 1:3 
                    es[1,gp] = Sgp[1,gp] # Sₑ som namn istället?
                    es[2,gp] = Sgp[2,gp]
                    es[3,gp] = Sgp[4,gp]
                end 

                c2tl6_f!(fₑ, coord[enod[ie][2:7],:], t, es, a[edof[ie,:]])
                Fᵢₙₜ[edof[ie,:]]           += fₑ            
            end

            bcval      = 0*bcval
            res        = Fᵢₙₜ - Fₑₓₜ
            res[bcdof] = 0*res[bcdof]
            residual   = norm(res,2)
            println("Iteration: ", iter, " Residual: ", residual)
        end
    end
    return K, K2
end