using LinearSolve, LinearSolvePardiso, SparseArrays, 
      StaticArrays, 
      IterativeSolvers, AlgebraicMultigrid, IncompleteLU    

#ENV["PATH"]

include("..//mesh_reader.jl")

#grid1 = createBoxMesh("box_1",0.0,0.0,1.0,1.0,0.1)

function load_files()
  
      include("..//material.jl")
  
      include("..//fem.jl")
  
      include("assemElemLin.jl")
  
      include("assemLin.jl")
  
      include("..//sensitivities.jl")

      include("..//run.jl")

      include("..//mma.jl")

      include("initLin.jl")
end
  
load_files()

function solver(dh,coord)
      imax     = 25
      TOL      = 1e-6
      residual = 0.0
      iter     = 1
  
      ndof     = size(coord,1)*2 
      nelm     = size(enod,1)

      # ----------- #
      # Read "grid" #
      # ----------- #
      K  = create_sparsity_pattern(dh)
  
      #  -------- #
      # Convert   #
      #  -------- #
      # coord <-- dh.grid.nodes 
      # enod  <-- dh.grid.cells
  
      #  ----- #
      # Init   #
      #  ----- #
      Fᵢₙₜ        = zeros(ndof)
      Fₑₓₜ        = zeros(ndof)
      a           = zeros(ndof)
      Δa          = zeros(ndof)
      res         = zeros(ndof)
      bcdof,bcval = setBCLin(0.0,dh)
      pdofs       = bcdof
      fdofs       = setdiff(1:ndof,pdofs)
      # ---------- #
      # Set params # // Kanske som input till solver???
      # ---------- #
      mp       = [175 80.769230769230759]
      t        = 1.0
  
      bcval₀   = bcval
  
      for n ∈ 1 : 10

          τ     = [1;0].*n

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
  
              assemGlobal!(K,Fᵢₙₜ,dh,mp,t,a,coord,enod,Γt,τ)
          
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
   
ie = 0
τ = [1.0 ; 0]
for cell in CellIterator(dh)
      kₑ = zeros(6,6)
      fₑ = zeros(6)
      ie = cellid(cell)
      for face in 1:nfaces(cell)
          if (cellid(cell), face) in Γt
              face_nods = [ Ferrite.facedof_indices(ip)[face][1]; Ferrite.facedof_indices(ip)[face][2] ]
              face_dofs = [ face_nods[1]*2-1; face_nods[1]*2; face_nods[2]*2-1; face_nods[2]*2 ]
              X = coord[ enod[ie][face_nods.+1] ,: ]
              fₑ[face_dofs] += tractionLoad( X, τ )
          end    
      end
end

