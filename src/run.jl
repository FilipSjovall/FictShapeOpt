include("mesh_reader.jl")

include("material.jl")

include("element_routines.jl")

include("fem.jl")
filename = "mesh_fine.txt"

coord, enod, edof = readAscii(filename);


#function solver()
    imax     = 25
    TOL      = 1e-6
    residual = 0.0
    iter     = 1

    ndof     = size(coord,1)*2 
    nelm     = size(enod,1)

    K        = zeros(ndof,ndof)
    Fᵢₙₜ     = zeros(ndof)
    Fₑₓₜ     = zeros(ndof)
    a        = zeros(ndof)
    Δa       = zeros(ndof)
    res      = zeros(ndof)

    F        = zeros(3,2,2)
    ef       = zeros(4,3)
    D        = zeros(3,3,3)
    Sgp        = zeros(4,3)
    es       = zeros(3,3)
    kₑ       = zeros(12,12)
    fₑ       = zeros(12)

    bcval    = zeros(6)
    bcval    = [0.0; 0.0; 0.05; 0.05; 0.0; 0.0]
    bcdof    = Vector{Int64}(undef,6)
    bcdof    = [3; 4; 7; 9; 5; 6]
    mp       = [175 80.769230769230759]
    t        = 1.0

    #for n ∈ 1 : 10
        res   = res.*0
        #bcval = 0.0
        println("Starting equillibrium iteration at loadstep: ",n)

        #while (iter < imax && residual > TOL ) || iter < 2
            iter += 1
            K     = K.*0
            Fᵢₙₜ  = Fᵢₙₜ .*0

            for ie ∈ 1:nelm
                kₑ       = kₑ.*0
                fₑ       = fₑ.*0
                c2tl6_d!( F, a[edof[ie,:]], coord[enod[ie][2:7],:] )
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

                c2tl6_e!(kₑ, coord[enod[ie][2:7],:], t, D,  a[edof[ie,:]], es)
                c2tl6_f!(fₑ, coord[enod[ie][2:7],:], t, es, a[edof[ie,:]])
                
                #println(kₑ)    
                #println(D[:,:,1])

                K[edof[ie,:],edof[ie,:]]   .= K[edof[ie,:],edof[ie,:]] .+ kₑ # += för mer kompakt kod? 
                Fᵢₙₜ[edof[ie,:]]           = Fᵢₙₜ[edof[ie,:]] + fₑ            
            end

            res = Fᵢₙₜ - Fₑₓₜ
            solveq!(Δa,K,-res,bcdof,bcval)
        #end
    #end
#end