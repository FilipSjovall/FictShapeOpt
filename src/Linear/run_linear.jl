#
#
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

#
function solver_C(dh, coord)

    # ---------- #
    # Set params # // Kanske som input till solver???
    # ---------- # // definiera mp här? och kanske ε ? iofs snyggare utanför!  
    t = 1.0

    # Define material parameters
    mp = [210 0.3] # [E ν]

    # ------------- #
    # Init-stuff    #
    # ------------- #
    imax = 25
    TOL = 1e-8
    residual = 0.0
    iter = 1
    # ------------- #.0
    # ------------- #
    #K = create_sparsity_pattern(dh)

    # ------ #
    #  Init  #
    # ------ #
    global Fᵢₙₜ = zeros(dh.ndofs.x)
    global Fₑₓₜ = zeros(dh.ndofs.x)
    global a    = zeros(dh.ndofs.x)
    global Δa   = zeros(dh.ndofs.x)
    global res  = zeros(dh.ndofs.x)
    global K    = create_sparsity_pattern(dh)
    # ---------- #
    # Set BCS    # 
    # ---------- # 
    # Set bcs - should be moved outside this function
    bcdof_top, bcval_top = setBCXY(-0.01, dh, Γ_top)
    bcdof_bot, bcval_bot = setBCXY(0.0, dh, Γ_bot)
    bcdof = [bcdof_top; bcdof_bot]
    bcval = [bcval_top; bcval_bot]

    ϵᵢⱼₖ = sortperm(bcdof)
    bcdof = bcdof[ϵᵢⱼₖ]
    bcval = bcval[ϵᵢⱼₖ]

    # - For Linear solver..
    pdofs = bcdof
    fdofs = setdiff(1:dh.ndofs.x, pdofs)

    bcval₀ = bcval

    for loadstep ∈ 1 : 1
        res = res .* 0
        bcval = bcval₀
        residual = 0 * residual
        iter = 0
        fill!(Δa, 0.0)
        println("Starting equillibrium iteration at loadstep: ", loadstep)

        # # # # # # # # # #
        # Newton solve.   #
        # # # # # # # # # #
        while (iter < imax && residual > TOL) || iter < 2
            iter += 1
            a += Δa
            assemGlobal!(K, Fᵢₙₜ, dh, mp, t, a, coord, enod, ε)
            solveq!(Δa, K, -Fᵢₙₜ, bcdof, bcval)
            bcval = 0 * bcval
            res = Fᵢₙₜ - Fₑₓₜ
            res[bcdof] = 0 * res[bcdof]
            residual = norm(res, 2)
            println("Iteration: ", iter, " Residual: ", residual)

            postprocess_opt(a, dh, "contact_mesh" * string(loadstep))
            #postprocess_opt(Fᵢₙₜ, dh, "contact_mesh" * string(loadstep))
            #σx, σy = StressExtract(dh, a, mp)
            #vtk_grid("contact" * string(loadstep), dh) do vtkfile
            #    vtk_point_data(vtkfile, dh, a) # displacement field
            #    vtk_point_data(vtkfile, σx, "σx")
            #    vtk_point_data(vtkfile, σy, "σy")
            #end
        end
    end
    fill!(Fₑₓₜ, 0.0)
    Fₑₓₜ[bcdof] = -Fᵢₙₜ[bcdof]

    return a, dh, Fₑₓₜ, Fᵢₙₜ, K
    
end

function fictitious_solver_C(d, dh0, coord₀)
    # allt överflödigt bör vid tillfälle flyttas utanför 
    # lösare till ett "init-liknande script så att huvudsaklig kod hålls ren
    imax = 25
    TOL = 1e-10
    residual = 0.0
    iter = 1
    global λ
    ndof = size(coord₀, 1) * 2
    nelm = size(enod, 1)



    #  ----- #
    # Init   #
    #  ----- #
    global Kψ = create_sparsity_pattern(dh)
    global Ψ = zeros(dh.ndofs.x)
    global Fᵢₙₜ = zeros(dh.ndofs.x)
    global Fₑₓₜ = zeros(dh.ndofs.x)
    global Ψ = zeros(dh.ndofs.x)
    global ΔΨ = zeros(dh.ndofs.x)
    global res = zeros(dh.ndofs.x)
    res = zeros(ndof)
    bcdof_top, bcval_top = setBCXY(0.0, dh, Γ_top)
    bcdof_bot, bcval_bot = setBCXY(0.0, dh, Γ_bot)
    bcdof = [bcdof_top; bcdof_bot]
    bcval = [bcval_top; bcval_bot]

    ϵᵢⱼₖ = sortperm(bcdof)
    bcdof = bcdof[ϵᵢⱼₖ]
    bcval = bcval[ϵᵢⱼₖ]

    # Struct - problem {dh,bcs,mp}

    pdofs = bcdof
    fdofs = setdiff(1:ndof, pdofs)
    # ---------- #
    # Set params # // Kanske som input till solver???
    # ---------- #

    bcval₀ = bcval

    for n ∈ 1 : 10
        res = res .* 0
        bcval = bcval₀
        residual = 0 * residual
        iter = 0
        λ = 0.1 * n
        fill!(ΔΨ, 0.0)

        println("Starting equillibrium iteration at loadstep: ", n)

        # # # # # # # # # #
        # Newton solve.  #
        # # # # # # # # # #
        while (iter < imax && residual > TOL) || iter < 2
            iter += 1
            Ψ += ΔΨ
            assemGlobal!(Kψ, Fᵢₙₜ, dh0, mp₀, t, Ψ, coord₀, enod, λ, d, Γ_robin)
            solveq!(ΔΨ, Kψ, -Fᵢₙₜ, bcdof, bcval)
            bcval = bcval .* 0
            res = Fᵢₙₜ #- Fₑₓₜ
            res[bcdof] = res[bcdof] .* 0
            residual = norm(res, 2)
            Ψ[bcdof] = bcval
            println("Iteration: ", iter, " Residual: ", residual, " λ: ", λ)
        end
    end
    return Ψ, dh0, Kψ, Fᵢₙₜ, λ
end
