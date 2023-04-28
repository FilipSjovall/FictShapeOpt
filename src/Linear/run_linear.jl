

function solver(dh,coord)
      # ------------- #
      # Init-stuff    #
      # ------------- #
      imax     = 25
      TOL      = 1e-8
      residual = 0.0
      iter     = 1
      ndof     = size(coord,1)*2 
      # ------------- #
      # Assem pattern #
      # ------------- #
      K  = create_sparsity_pattern(dh)
      # ------ #
      #  Init  #
      # ------ #
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
      t        = 1.0
      bcval₀   = bcval
      for n ∈ 1 : 10
          τ        = [0.1;0.1].*n
          res      = res.*0
          bcval    = bcval₀
          residual = 0*residual
          iter     = 0
          fill!(Δa,0.0)
          println("Starting equillibrium iteration at loadstep: ",n)
          # # # # # # # # # #
          # Newton solve.  #
          # # # # # # # # # #
          while (iter < imax && residual > TOL ) || iter < 2
              iter      += 1
              a         += Δa
              assemGlobal!(K,Fᵢₙₜ,dh,mp,t,a,coord,enod,Γt,τ)
              solveq!(Δa, K, -Fᵢₙₜ, bcdof, bcval)
              bcval      = 0*bcval
              res        = Fᵢₙₜ - Fₑₓₜ
              res[bcdof] = 0*res[bcdof]
              residual   = norm(res,2)
              println("Iteration: ", iter, " Residual: ", residual)
          end
      end
      fill!(Fₑₓₜ,0.0)
      τ        = [0.1;0.1].*n
      assemGlobal!(Fₑₓₜ,dh,t,a,coord,enod,Γt,τ)
      Fₑₓₜ[bcdof] = - Fᵢₙₜ[bcdof]
      return a, dh, Fₑₓₜ, Fᵢₙₜ, K
end
   
function fictitious_solver(d,dh0,coord₀)
      # allt överflödigt bör vid tillfälle flyttas utanför 
      # lösare till ett "init-liknande script så att huvudsaklig kod hålls ren
      imax     = 25
      TOL      = 1e-10
      residual = 0.0
      iter     = 1
      global λ
      ndof     = size(coord₀,1)*2 
      nelm     = size(enod,1)
  
      Kψ       = create_sparsity_pattern(dh0)
  
      #  ----- #
      # Init   #
      #  ----- #
      Fᵢₙₜ        = zeros(ndof)
      Fₑₓₜ        = zeros(ndof)
      Ψ           = zeros(ndof)
      ΔΨ          = zeros(ndof)
      res         = zeros(ndof)
      bcdof,bcval = setBCLin(0.0,dh0) # Ha bc som argument?
  
      # Struct - problem {dh,bcs,mp}

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
              assemGlobal!(Kψ,Fᵢₙₜ,dh0,mp₀,t,Ψ,coord₀,enod,λ,d,Γ_robin)
              solveq!(ΔΨ, Kψ, -Fᵢₙₜ, bcdof, bcval)
              bcval      = bcval.*0
              res        = Fᵢₙₜ #- Fₑₓₜ
              res[bcdof] = res[bcdof].*0
              residual   = norm(res,2)
              Ψ[bcdof]   = bcval;
              println("Iteration: ", iter, " Residual: ", residual, " λ: ", λ)
          end
      end
      return Ψ, dh0, Kψ, Fᵢₙₜ, λ
  end

