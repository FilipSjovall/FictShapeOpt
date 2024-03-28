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
      for loadstep ∈ 1 : 10
          τ        = [0.1;0.1].*n
          res      = res.*0
          bcval    = bcval₀
          residual = 0*residual
          iter     = 0
          fill!(Δa,0.0)
          println("Starting equilibrium iteration at loadstep: ",n)
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
function solver_C(dh, coord, Δ, nloadsteps)

    # ---------- #
    # Set params # // Kanske som input till solver???
    # ---------- # // definiera mp här? och kanske ε ? iofs snyggare utanför!
    t = 1.0

    # ------------- #
    # Init-stuff    #
    # ------------- #
    imax = 200
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
    global rc   = zeros(dh.ndofs.x)
    global Fₑₓₜ = zeros(dh.ndofs.x)
    global a    = zeros(dh.ndofs.x)
    global Δa   = zeros(dh.ndofs.x)
    global res  = zeros(dh.ndofs.x)
    global K    = create_sparsity_pattern(dh)
    # ---------- #
    # Set BCS    #
    # ---------- #
    # Set bcs - should be moved outside this function
    #bcdof_top, bcval_top = setBCXY_both(Δ / nloadsteps, dh, Γ_top)
    #bcdof_bot, bcval_bot = setBCXY_both(0.0, dh, Γ_bot)
    bcdof_top, bcval_top   = setBCXY(Δ/nloadsteps, dh, n_top)
    bcdof_bot, bcval_bot   = setBCXY(0.0, dh, n_bot)
    #bcdof_left, bcval_left = setBCX(0.0, dh, n_left)
    bcdofs        = [bcdof_top; bcdof_bot]
    bcvals        = [bcval_top; bcval_bot]
    ϵᵢⱼₖ         = sortperm(bcdofs)
    global bcdofs = bcdofs[ϵᵢⱼₖ]
    global bcvals = bcvals[ϵᵢⱼₖ]
    # - For Linear solver..
    global pdofs = bcdofs
    global fdofs = setdiff(1:dh.ndofs.x, pdofs)
    bcval₀ = bcvals
    global β = 1.0
    loadstep = 0
    while loadstep < nloadsteps
        loadstep += 1
        #global ε = ε * 1.1
        res = res .* 0
        bcvals = bcval₀
        residual = 0 * residual
        iter = 0
        fill!(Δa, 0.0)
        print("\n", "Starting equilibrium iteration at loadstep: ", loadstep, "\n")
        a_old = a
        # # # # # # # # # #
        # Newton solve.   #
        # # # # # # # # # #
            while  residual > TOL || iter < 2
                iter += 1
                if iter % 20 == 0 || norm(res) > 1e3
                    a = a_old
                    bcvals = bcval₀
                    if β > 1/8
                        global β = β * 0.5
                        Δ_remaining = (Δ*nloadsteps - β * Δ - loadstep * Δ)/nloadsteps
                        remaining_steps = nloadsteps - loadstep
                        nloadsteps = loadstep + 2remaining_steps + (1 / β - 1)
                        bcvals = bcvals ./2 #
                        bcval₀= bcvals
                    end
                    fill!(Δa, 0.0)
                    println("Penalty paremeter and updated: $ε, and step length $β ")
                end

                #a += β * Δa
                a += Δa
                assemGlobal!(K, Fᵢₙₜ, dh, mp, t, a, coord, enod, ε)
                solveq!(Δa,  K, -Fᵢₙₜ, bcdofs, bcvals)
                bcvals = 0 * bcvals
                res = Fᵢₙₜ - Fₑₓₜ
                res[bcdofs] = 0 * res[bcdofs]
                residual = norm(res, 2)
                @printf "Iteration: %i | Residual: %.4e | Δ: %.4f \n" iter residual a[bcdofs[2]]
            end
            X_c,tract = plotTraction()
            p5 = plot(X_c, tract, label="λ" , marker=4, lc=:tomato, mc=:tomato, grid=false, legend=:outerleft, ylims = (0, 1.2*maximum(tract)))
            display(p5)
            σx, σy = StressExtract(dh, a, mp)
            vtk_grid("results/contact" * string(loadstep), dh) do vtkfile
                #vtk_grid("contact" * string(iter), dh) do vtkfile
                vtk_point_data(vtkfile, dh, a) # displacement field
                vtk_point_data(vtkfile, σx, "σx")
                vtk_point_data(vtkfile, σy, "σy")
            end
            fill!(Fₑₓₜ, 0.0)
            Fₑₓₜ[bcdofs] = -Fᵢₙₜ[bcdofs]
    end
    τ_c = ExtractContactTraction(a, ε, coord)
    return a, dh, Fₑₓₜ, Fᵢₙₜ, K, τ_c
end
#
# Fictitious equillibrium for shape optimization of problem with contact
function fictitious_solver_C(d, dh0, coord₀)
    # allt överflödigt bör vid tillfälle flyttas utanför
    # lösare till ett "init-liknande script så att huvudsaklig kod hålls ren
    imax = 100
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
    bcdof_top, bcval_top = setBCXY(0.0, dh, n_top)
    bcdof_bot, bcval_bot = setBCXY(0.0, dh, n_bot)
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

    for loadstep ∈ 1 : 10
        res = res .* 0
        bcval = bcval₀
        residual = 0 * residual
        iter = 0
        λ = 0.1 * loadstep
        fill!(ΔΨ, 0.0)

        println("Starting equilibrium iteration at loadstep: ", loadstep)


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
            postprocess_opt(Ψ, dh, "results/fict_def" * string(loadstep))
        end

    end
    return Ψ, dh0, Kψ, Fᵢₙₜ, λ
end
#
# Solver for hertz contact
function solver_C2(dh, coord)

    # ---------- #
    # Set params # // Kanske som input till solver???
    # ---------- # // definiera mp här? och kanske ε ? iofs snyggare utanför!
    t = 1.0

    # Define material parameters
    mp = [175 80.769230769230759]


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
    global rc = zeros(dh.ndofs.x)
    global Fₑₓₜ = zeros(dh.ndofs.x)
    global a = zeros(dh.ndofs.x)
    global Δa = zeros(dh.ndofs.x)
    global res = zeros(dh.ndofs.x)
    global K = create_sparsity_pattern(dh)
    # ---------- #
    # Set BCS    #
    # ---------- #
    # Set bcs - should be moved outside this function
    bcdof_top, bcval_top = setBCXY(-0.05, dh, Γ_top)
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

    for loadstep ∈ 1 : 10
        τ   = [0.0; 0.0001] * (loadstep-1)
        res = res .* 0
        bcval = bcval₀
        residual = 0 * residual
        iter = 0
        fill!(Δa, 0.0)
        println("Starting equilibrium iteration at loadstep: ", loadstep)

        # # # # # # # # # #
        # Newton solve.   #
        # # # # # # # # # #
        while (iter < imax && residual > TOL) || iter < 2
            iter += 1
            a += Δa
            assemGlobal!(K, Fᵢₙₜ, rc, dh, mp, t, a, coord, enod, ε, Γ_top, τ)
            solveq!(Δa, K, -Fᵢₙₜ, bcdof, bcval)
            bcval = 0 * bcval
            res = Fᵢₙₜ - Fₑₓₜ
            res[bcdof] = 0 * res[bcdof]
            residual = norm(res, 2)
            println("Iteration: ", iter, " Residual: ", residual)
        end
        σx, σy = StressExtract(dh, a, mp)
        vtk_grid("hertz" * string(loadstep), dh) do vtkfile
            vtk_point_data(vtkfile, dh, a) # displacement field
            vtk_point_data(vtkfile, σx, "σx")
            vtk_point_data(vtkfile, σy, "σy")
        end
    end
    fill!(Fₑₓₜ, 0.0)
    Fₑₓₜ[bcdof] = -Fᵢₙₜ[bcdof]
    τ_c = ExtractContactTraction(a, ε, coord)
    return a, dh, Fₑₓₜ, Fᵢₙₜ, K, τ_c

end
#
# Fictitious equillibrium for shape optimization with consistent with contact
function fictitious_solver_with_contact(d, dh0, coord₀, nloadsteps)
    # allt överflödigt bör vid tillfälle flyttas utanför
    # lösare till ett "init-liknande script så att huvudsaklig kod hålls ren
    imax = 100
    TOL = 1e-10
    residual = 0.0
    iter = 1
    global λ = 0
    ndof = size(coord₀, 1) * 2
    nelm = size(enod, 1)


    t = 1.0
    #  ----- #
    # Init   #
    #  ----- #
    global Kψ = create_sparsity_pattern(dh0)
    global Ψ = zeros(dh0.ndofs.x)
    global FΨ = zeros(dh0.ndofs.x)
    global Fₑₓₜ = zeros(dh0.ndofs.x)
    global Ψ = zeros(dh0.ndofs.x)
    global ΔΨ = zeros(dh0.ndofs.x)
    global res = zeros(dh0.ndofs.x)

    #bcdof_top_o2, _ = setBCXY_both(0.0, dh, Γ_top)
    #bcdof_bot_o2, _ = setBCXY_both(0.0, dh, Γ_bot)
    bcdof_top_o2, _  = setBCXY(0.0, dh, n_top)
    bcdof_bot_o2, _  = setBCXY(0.0, dh, n_bot)
    #bcdof_left_o2, _ = setBCX(0.0, dh, n_left)
    bcdof_o2         = [bcdof_top_o2; bcdof_bot_o2]
    ϵᵢⱼₖ            = sortperm(bcdof_o2)
    global bcdof_o2  = bcdof_o2[ϵᵢⱼₖ]
    global bcval_o2  = bcdof_o2 .* 0.0

    # Struct - problem {dh0,bcs,mp}

    global pdofs = bcdof
    global fdofs = setdiff(1:ndof, pdofs)

    # ---------- #
    # Set params # // Kanske som input till solver???
    # ---------- #

    bcval₀_o2 = bcval_o2
    Δλ = (1.0 / nloadsteps)
    #Δλ₀ = Δλ

    #for loadstep ∈ 1 : nloadsteps
    ##
    loadstep = 0
    while loadstep < nloadsteps
        loadstep +=1
        #if Δλ >  0.1 * 1/8
        #    global μ = μ * 1.1
        #end
    ##
        res = res .* 0
        bcval_o2 = bcval₀_o2
        residual = 0 * residual
        iter = 0
        global λ += Δλ #* loadstep
        fill!(ΔΨ, 0.0)
        print("\n","Starting equilibrium iteration at loadstep: ", loadstep, "\n")
        Ψ_old = Ψ
        # # # # # # # # # #
        # Newton solve.  #
        # # # # # # # # # #
        while  residual > TOL || iter < 2
            iter += 1
            if iter % 10 == 0 || norm(res) > 1e2 #&& Δλ > 1/16
                Ψ = Ψ_old
                #if Δλ > 0.1 * 1/8
                    global λ -= Δλ #* loadstep
                    Δλ        = Δλ/2
                    global λ += Δλ  #* loadstep
                    remaining_steps = nloadsteps - loadstep
                    #nloadsteps = loadstep + 2remaining_steps +  Δλ₀ / Δλ  - 1
                    nloadsteps = loadstep + round((1 - λ ) / Δλ)
                #else
                    # global μ    = μ * 1.1#0.9
                #end
                fill!(ΔΨ, 0.0)
                println("Step length updated: $Δλ, penalty parameter: $μ")
            end

            Ψ    += ΔΨ
            assemGlobal!(Kψ, FΨ, dh0, mp₀, t, Ψ, coord₀, enod, λ, d, Γ_robin, μ)
            solveq!(ΔΨ, Kψ, -FΨ, bcdof_o2, bcval_o2)
            bcval_o2      = bcval_o2 .* 0
            res           = FΨ #- Fₑₓₜ
            res[bcdof_o2] = res[bcdof_o2] .* 0
            residual      = norm(res, 2)
            Ψ[bcdof_o2]   = bcval_o2
            if loadstep < 40 && iter < 20
                postprocess_opt(Ψ, dh0, "results/fictitious_flat" * string(loadstep))


            end
            if iter < 20
                postprocess_opt(res, dh0, "results/fictres_flat" * string(iter))
                postprocess_opt(Ψ + ΔΨ, dh0, "results/fictitious_flat" * string(iter))
            end
            @printf "Iteration: %i | Residual: %.4e | λ: %.4f \n" iter residual λ
        end
    end
    return Ψ, dh0, Kψ, FΨ, λ
end
#
#
function fictitious_solver_with_contact_hook(d, dh0, coord₀, nloadsteps)
    # allt överflödigt bör vid tillfälle flyttas utanför
    # lösare till ett "init-liknande script så att huvudsaklig kod hålls ren
    TOL      = 1e-10
    residual = 0.0
    iter     = 1
    global λ = 0
    ndof     = size(coord₀, 1) * 2
    nelm     = size(enod, 1)
    t        = 1.0

    #  ----- #
    # Init   #
    #  ----- #
    global Kψ  = create_sparsity_pattern(dh0)
    global Ψ   = zeros(dh0.ndofs.x)
    global FΨ  = zeros(dh0.ndofs.x)
    global Ψ   = zeros(dh0.ndofs.x)
    global ΔΨ  = zeros(dh0.ndofs.x)
    global res = zeros(dh0.ndofs.x)

    global bcdof_o2 = bcdofs_opt
    global bcval_o2 = bcdofs_opt .* 0.0
    global pdofs    = bcdofs_opt
    global fdofs    = setdiff(1:ndof, pdofs)

    bcval₀_o2 = bcval_opt
    Δλ = (1.0 / nloadsteps)
    loadstep  = 0

    while loadstep < nloadsteps
        loadstep += 1
        res       = res .* 0
        bcval_opt = bcval₀_o2
        residual  = 0 * residual
        iter      = 0
        global λ += Δλ #* loadstep
        fill!(ΔΨ, 0.0)
        print("\n", "Starting equilibrium iteration at loadstep: ", loadstep, "\n")
        Ψ_old = Ψ

        # # # # # # # # # #
        # Newton solve.  #
        # # # # # # # # # #
        while residual > TOL || iter < 2
            iter += 1
            if iter % 20 == 0 || norm(res) > 1e2 #&& Δλ > 1/16
                Ψ = Ψ_old
                if Δλ > 0.1 * 1 / 64
                    global λ -= Δλ #* loadstep
                    Δλ = Δλ / 2
                    global λ += Δλ  #* loadstep
                    remaining_steps = nloadsteps - loadstep
                    nloadsteps      = loadstep + round((1 - λ) / Δλ)
                end
                fill!(ΔΨ, 0.0)
                println("Step length updated: $Δλ, penalty parameter: $μ")
            end

            Ψ += ΔΨ
            assemGlobal!(Kψ, FΨ, dh0, mp₀, t, Ψ, coord₀, enod, λ, d, Γ_robin, μ)
            solveq!(ΔΨ, Kψ, -FΨ, bcdofs_opt, bcval_opt)
            #
            bcval_opt       = bcval_opt .* 0
            res             = FΨ #- Fₑₓₜ
            res[bcdofs_opt] = res[bcdofs_opt] .* 0
            residual        = norm(res, 2)
            Ψ[bcdofs_opt]  .= 0.0
            if loadstep < 40 && iter < 20
                postprocess_opt(Ψ, dh0, "results/fictitious_t2" * string(loadstep))
            end
            if iter < 20
                postprocess_opt(res, dh0, "results/fictres_t2" * string(iter))
                postprocess_opt(Ψ, dh0, "results/fictitious_iter_t2" * string(iter))
            end
            @printf "Iteration: %i | Residual: %.4e | λ: %.4f \n" iter residual λ
            #Wψ = energy(dh0, Ψ, mp₀)
            #σψx, σψy = StressExtract(dh0, Ψ, mp₀)
            #@save "filter forces fat" FΨ Wψ σψx σψy
        end
        if loadstep < 11
            Ψ_hist[:, loadstep] = Ψ
            d_hist[:, loadstep] = d
        end
    end
    return Ψ, dh0, Kψ, FΨ, λ, Ψ_hist, d_hist
end
#
#
function solver_C_hook(dh, coord, Δ, nloadsteps)

    # ---------- #
    # Set params # // Kanske som input till solver???
    # ---------- # // definiera mp här? och kanske ε ? iofs snyggare utanför!
    t = 1.0

    # Define material parameters
    mp = [175 80.769230769230759]

    # ------------- #
    # Init-stuff    #
    # ------------- #
    imax = 200
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
    global rc    = zeros(dh.ndofs.x)
    global Fₑₓₜ  = zeros(dh.ndofs.x)
    global a     = zeros(dh.ndofs.x)
    global Δa    = zeros(dh.ndofs.x)
    global res   = zeros(dh.ndofs.x)
    global K     = create_sparsity_pattern(dh)

    # ------------------- #
    # Boundary conditions #
    # ------------------- #
   #=
    bcdof_left, bcval_left    = setBCXY_both(0.0, dh, n_left)
    bcdof_right, bcval_right  = setBCXY_both(Δ/nloadsteps, dh, n_right)
    bcdofs                     = [bcdof_left; bcdof_right]
    bcvals                     = [bcval_left; bcval_right]
   =#
   ## #=
     bcdof_left, bcval_left     = setBCXY_X(  0.0, dh, n_left)
    #bcdof_left, bcval_left     = setBCXY_X( -Δ / nloadsteps, dh, n_left)
    bcdof_right, bcval_right   = setBCXY_X(  Δ / nloadsteps, dh, n_right)
    bcdof_bot, bcval_bot       = setBCY(0.0, dh, n_bot)
    bcdof_top, bcval_top       = setBCY(0.0, dh, n_top)

    bcdof_bot, bcval_bot       = Vector{Int64}(), Vector{Float64}()
    bcdof_top, bcval_top       = Vector{Int64}(), Vector{Float64}()

    bcdofs                     = [bcdof_left; bcdof_right; bcdof_bot; bcdof_top]
    bcvals                     = [bcval_left; bcval_right; bcval_bot; bcval_top]
    ## =#
    ϵᵢⱼₖ                      = sortperm(bcdofs)
    global bcdofs              = bcdofs[ϵᵢⱼₖ]
    global bcvals              = bcvals[ϵᵢⱼₖ]

    # - For Linear solver..
    global pdofs = bcdofs
    global fdofs = setdiff(1:dh.ndofs.x, pdofs)

    bcval₀ = bcvals
    global β = 1.0
    #for loadstep ∈ 1 : nloadsteps
    ##
    loadstep = 0
    while loadstep < nloadsteps
        loadstep += 1
        #global ε = ε * 1.1
        ##
        res = res .* 0
        bcvals = bcval₀
        residual = 0 * residual
        iter = 0
        fill!(Δa, 0.0)
        print("\n", "Starting equilibrium iteration at loadstep: ", loadstep, "\n")
        #global ε = ε₀
        a_old = a
        # # # # # # # # # #
        # Newton solve.   #
        # # # # # # # # # #

        while residual > TOL || iter < 2
            iter += 1
            if iter % 20 == 0 || norm(res) > 1e3 && β > 1 / 8
                a = a_old
                bcvals = bcval₀
                global β = β * 0.5
                Δ_remaining = (Δ * nloadsteps - β * Δ - loadstep * Δ) / nloadsteps
                remaining_steps = nloadsteps - loadstep
                nloadsteps = loadstep + 2remaining_steps + (1 / β - 1)
                bcvals = bcvals ./ 2 #
                bcval₀ = bcvals
                fill!(Δa, 0.0)
                println("Penalty paremeter and updated: $ε, and step length $β ")
            end

            a += Δa
            assemGlobal!(K, Fᵢₙₜ, dh, mp, t, a, coord, enod, ε)
            solveq!(Δa, K, -Fᵢₙₜ, bcdofs, bcvals)
            bcvals = 0 * bcvals
            res = Fᵢₙₜ - Fₑₓₜ
            res[bcdofs] = 0 * res[bcdofs]
            residual = norm(res, 2)
            @printf "Iteration: %i | Residual: %.4e | Δ: %.4f \n" iter residual a[bcdof_right[1]]
        end
        if loadstep < 40 && iter < 20
            σx, σy = StressExtract(dh, a, mp)
            vtk_grid("results/🍌-contact" * string(loadstep), dh) do vtkfile
                vtk_point_data(vtkfile, dh, a) # displacement field
                vtk_point_data(vtkfile, σx, "σx")
                vtk_point_data(vtkfile, σy, "σy")
            end
        end
        Fₑₓₜ[bcdofs] = -Fᵢₙₜ[bcdofs]
        if loadstep < 11
            a_hist[:,loadstep] = a
        end
    end
    τ_c = ExtractContactTraction(a, ε, coord)
    return a, dh, Fₑₓₜ, Fᵢₙₜ, K, τ_c, a_hist
end

function fictitious_solver_with_contact_half(d, dh0, coord₀, nloadsteps)
    # allt överflödigt bör vid tillfälle flyttas utanför
    # lösare till ett "init-liknande script så att huvudsaklig kod hålls ren
    imax = 100
    TOL = 1e-10
    residual = 0.0
    iter = 1
    global λ = 0
    ndof = size(coord₀, 1) * 2
    nelm = size(enod, 1)


    t = 1.0
    #  ----- #
    # Init   #
    #  ----- #
    #global Kψ = create_sparsity_pattern(dh0)
    global FΨ = zeros(dh0.ndofs.x)
    global Fₑₓₜ = zeros(dh0.ndofs.x)
    global Ψ = zeros(dh0.ndofs.x)
    global ΔΨ = zeros(dh0.ndofs.x)
    global res = zeros(dh0.ndofs.x)

    #bcdof_top_o2, _ = setBCXY_both(0.0, dh, Γ_top)
    #bcdof_bot_o2, _ = setBCXY_both(0.0, dh, Γ_bot)
    bcdof_top_o2, _  = setBCY(0.0, dh, n_top)
    bcdof_bot_o2, _  = setBCY(0.0, dh, n_bot)
    bcdof_left_o2, _ = setBCX(0.0, dh, n_left)
    bcdof_o2 = [bcdof_top_o2; bcdof_bot_o2; bcdof_left_o2]
    ϵᵢⱼₖ = sortperm(bcdof_o2)
    global bcdof_o2 = bcdof_o2[ϵᵢⱼₖ]
    global bcval_o2 = bcdof_o2 .* 0.0

    # ---------- #
    # Set params # // Kanske som input till solver???
    # ---------- #

    bcval₀_o2 = bcval_o2
    n₀ = nloadsteps
    Δλ = (1.0 / nloadsteps)
    #Δλ₀ = Δλ

    #for loadstep ∈ 1 : nloadsteps
    ##
    loadstep = 0
    while loadstep < nloadsteps
        loadstep += 1
        #if Δλ >  0.1 * 1/8
        #    global μ = μ * 1.1
        #end
        ##
        res = res .* 0
        bcval_o2 = bcval₀_o2
        residual = 0 * residual
        iter = 0
        global λ += Δλ #* loadstep
        fill!(ΔΨ, 0.0)
        print("\n", "Starting equilibrium iteration at loadstep: ", loadstep, "\n")
        Ψ_old = Ψ
        # # # # # # # # # #
        # Newton solve.  #
        # # # # # # # # # #
        while residual > TOL || iter < 2
            iter += 1
            if iter % 10 == 0 || norm(res) > 1e2 && Δλ > ((1.0 / n₀) * 1/16)
                Ψ = Ψ_old
                #if Δλ > 0.1 * 1/8
                global λ -= Δλ #* loadstep
                Δλ = Δλ / 2
                global λ += Δλ  #* loadstep
                remaining_steps = nloadsteps - loadstep
                #nloadsteps = loadstep + 2remaining_steps +  Δλ₀ / Δλ  - 1
                nloadsteps = loadstep + round((1 - λ) / Δλ)
                #else
                # global μ    = μ * 1.1#0.9
                #end
                fill!(ΔΨ, 0.0)
                println("Step length updated: $Δλ, penalty parameter: $μ")
            end

            Ψ += ΔΨ
            assemGlobal!(Kψ, FΨ, dh0, mp₀, t, Ψ, coord₀, enod, λ, d, Γ_robin, μ)
            solveq!(ΔΨ, Kψ, -FΨ, bcdof_o2, bcval_o2)
            bcval_o2 = bcval_o2 .* 0
            res = FΨ #- Fₑₓₜ
            res[bcdof_o2] = res[bcdof_o2] .* 0
            residual = norm(res, 2)
            Ψ[bcdof_o2] = bcval_o2
            if loadstep < 40
                postprocess_opt(Ψ, dh0, "results/fictitious" * string(loadstep))
                #postprocess_opt(Ψ, dh0, "results/fictitious" * string(iter))
            end
            @printf "Iteration: %i | Residual: %.4e | λ: %.4f \n" iter residual λ
        end
    end
    return Ψ, dh0, Kψ, FΨ, λ
end

function solver_C_half(dh, coord, Δ, nloadsteps)

    # ---------- #
    # Set params # // Kanske som input till solver???
    # ---------- # // definiera mp här? och kanske ε ? iofs snyggare utanför!
    t = 1.0

    # ------------- #
    # Init-stuff    #
    # ------------- #
    imax = 200
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
    global rc = zeros(dh.ndofs.x)
    global Fₑₓₜ = zeros(dh.ndofs.x)
    global a = zeros(dh.ndofs.x)
    global Δa = zeros(dh.ndofs.x)
    global res = zeros(dh.ndofs.x)
    #global K = create_sparsity_pattern(dh)
    # ---------- #
    # Set BCS    #
    # ---------- #
    # Set bcs - should be moved outside this function
    #bcdof_top, bcval_top = setBCXY_both(Δ / nloadsteps, dh, Γ_top)
    #bcdof_bot, bcval_bot = setBCXY_both(0.0, dh, Γ_bot)
    bcdof_top, bcval_top   = setBCY(Δ / nloadsteps, dh, n_top)
    bcdof_bot, bcval_bot   = setBCY(0.0, dh, n_bot)
    bcdof_left, bcval_left = setBCX(0.0, dh, n_left)
    bcdofs = [bcdof_top; bcdof_bot; bcdof_left]
    bcvals = [bcval_top; bcval_bot; bcval_left]
    ϵᵢⱼₖ = sortperm(bcdofs)
    global bcdofs = bcdofs[ϵᵢⱼₖ]
    global bcvals = bcvals[ϵᵢⱼₖ]
    # - For Linear solver..
    global pdofs = bcdofs
    global fdofs = setdiff(1:dh.ndofs.x, pdofs)
    bcval₀ = bcvals
    global β = 1.0
    loadstep = 0
    while loadstep < nloadsteps
        loadstep += 1
        #global ε = ε * 1.1
        res = res .* 0
        bcvals = bcval₀
        residual = 0 * residual
        iter = 0
        fill!(Δa, 0.0)
        print("\n", "Starting equilibrium iteration at loadstep: ", loadstep, "\n")
        a_old = a
        # # # # # # # # # #
        # Newton solve.   #
        # # # # # # # # # #
        while residual > TOL || iter < 2
            iter += 1
            if iter % 20 == 0 || norm(res) > 1e3
                a = a_old
                bcvals = bcval₀
                if β > 1 / 8
                    global β = β * 0.5
                    Δ_remaining = (Δ * nloadsteps - β * Δ - loadstep * Δ) / nloadsteps
                    remaining_steps = nloadsteps - loadstep
                    nloadsteps = loadstep + 2remaining_steps + (1 / β - 1)
                    bcvals = bcvals ./ 2 #
                    bcval₀ = bcvals
                end
                fill!(Δa, 0.0)
                println("Penalty paremeter and updated: $ε, and step length $β ")
            end

            #a += β * Δa
            a += Δa
            assemGlobal!(K, Fᵢₙₜ, rc, dh, mp, t, a, coord, enod, ε)
            solveq!(Δa, K, -Fᵢₙₜ, bcdofs, bcvals)
            bcvals = 0 * bcvals
            res = Fᵢₙₜ - Fₑₓₜ
            res[bcdofs] = 0 * res[bcdofs]
            residual = norm(res, 2)
            @printf "Iteration: %i | Residual: %.4e | Δ: %.4f \n" iter residual a[bcdofs[2]]
            σx, σy = StressExtract(dh, a, mp)
            vtk_grid("contact" * string(iter), dh) do vtkfile
                vtk_point_data(vtkfile, dh, a) # displacement field
                vtk_point_data(vtkfile, σx, "σx")
                vtk_point_data(vtkfile, σy, "σy")
            end
        end
        if loadstep == 10
            # Plot traction , can be moved to function...
            τ_c = ExtractContactTraction(a, ε, coord)
            traction = ExtractContactTraction(a, ε, coord)
            X_c = []
            tract = []
            for (key, val) ∈ traction
                append!(X_c, coord[key, 1])
                append!(tract, val)
            end
            ϵᵢⱼₖ = sortperm(X_c)
            tract = tract[ϵᵢⱼₖ]
            X_c = X_c[ϵᵢⱼₖ]
            p = plot(X_c, tract, legend=false, marker=4, lc=:tomato, mc=:tomato)
            display(p)
        end
        σx, σy = StressExtract(dh, a, mp)
        vtk_grid("results/contact" * string(loadstep), dh) do vtkfile
            #vtk_grid("contact" * string(iter), dh) do vtkfile
            vtk_point_data(vtkfile, dh, a) # displacement field
            vtk_point_data(vtkfile, σx, "σx")
            vtk_point_data(vtkfile, σy, "σy")
        end
        fill!(Fₑₓₜ, 0.0)
        Fₑₓₜ[bcdofs] = -Fᵢₙₜ[bcdofs]
    end
    τ_c = ExtractContactTraction(a, ε, coord)
    return a, dh, Fₑₓₜ, Fᵢₙₜ, K, τ_c
end

##
function fictitious_solver_with_contact_hook_half(d, dh0, coord₀, nloadsteps)
    # allt överflödigt bör vid tillfälle flyttas utanför
    # lösare till ett "init-liknande script så att huvudsaklig kod hålls ren
    TOL = 1e-10
    residual = 0.0
    iter = 1
    global λ = 0
    ndof = size(coord₀, 1) * 2
    nelm = size(enod, 1)
    t = 1.0

    #  ----- #
    # Init   #
    #  ----- #
    global Kψ = create_sparsity_pattern(dh0)
    global Ψ = zeros(dh0.ndofs.x)
    global FΨ = zeros(dh0.ndofs.x)
    global Ψ = zeros(dh0.ndofs.x)
    global ΔΨ = zeros(dh0.ndofs.x)
    global res = zeros(dh0.ndofs.x)

    global bcdof_o2 = bcdofs_opt
    global bcval_o2 = bcdofs_opt .* 0.0
    global pdofs = bcdofs_opt
    global fdofs = setdiff(1:ndof, pdofs)

    bcval₀_o2 = bcval_opt
    Δλ = (1.0 / nloadsteps)
    loadstep = 0

    while loadstep < nloadsteps
        loadstep += 1
        res = res .* 0
        bcval_opt = bcval₀_o2
        residual = 0 * residual
        iter = 0
        global λ += Δλ #* loadstep
        fill!(ΔΨ, 0.0)
        print("\n", "Starting equilibrium iteration at loadstep: ", loadstep, "\n")
        Ψ_old = Ψ

        # # # # # # # # # #
        # Newton solve.  #
        # # # # # # # # # #
        while residual > TOL || iter < 2
            iter += 1
            if iter % 20 == 0 || norm(res) > 1e2 #&& Δλ > 1/16
                Ψ = Ψ_old
                if Δλ > 0.1 * 1 / 64
                    global λ -= Δλ #* loadstep
                    Δλ = Δλ / 2
                    global λ += Δλ  #* loadstep
                    remaining_steps = nloadsteps - loadstep
                    nloadsteps = loadstep + round((1 - λ) / Δλ)
                else
                    global μ = μ * 0.9
                end
                fill!(ΔΨ, 0.0)
                println("Step length updated: $Δλ, penalty parameter: $μ")
            end

            Ψ += ΔΨ
            assemGlobal!(Kψ, FΨ, dh0, mp₀, t, Ψ, coord₀, enod, λ, d, Γ_robin, μ)
            solveq!(ΔΨ, Kψ, -FΨ, bcdofs_opt, bcval_opt)
            # assemGlobal!(Kψ, FΨ, dh0, mp₀, t, Ψ, coord₀, enod, λ, d, Γ_robin, μ)

            bcval_opt = bcval_opt .* 0
            res = FΨ #- Fₑₓₜ
            res[bcdofs_opt] = res[bcdofs_opt] .* 0
            residual = norm(res, 2)
            Ψ[bcdofs_opt] .= 0.0
            # if loadstep < 40 && iter < 20
            #     postprocess_opt(Ψ, dh0, "results/fictitious" * string(loadstep))
            # end
            if iter < 20
                postprocess_opt(res, dh0, "results/fictres" * string(iter))
                postprocess_opt(Ψ, dh0, "results/fictitious_iter" * string(iter))
            end
            @printf "Iteration: %i | Residual: %.4e | λ: %.4f \n" iter residual λ
        end
    end
    return Ψ, dh0, Kψ, FΨ, λ
end
#
#
function solver_C_hook_half(dh, coord, Δ, nloadsteps)

    # ---------- #
    # Set params # // Kanske som input till solver???
    # ---------- # // definiera mp här? och kanske ε ? iofs snyggare utanför!
    t = 1.0

    # ------------- #
    # Init-stuff    #
    # ------------- #
    imax = 200
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
    global rc = zeros(dh.ndofs.x)
    global Fₑₓₜ = zeros(dh.ndofs.x)
    global a = zeros(dh.ndofs.x)
    global Δa = zeros(dh.ndofs.x)
    global res = zeros(dh.ndofs.x)
    global K = create_sparsity_pattern(dh)

    # ------------------- #
    # Boundary conditions #
    # ------------------- #
    # L structure
    bcdof_left, bcval_left = setBCXY_X(0.0, dh, n_left)
    bcdof_top, bcval_top   = setBCXY_both(0.0, dh, n_top)
    #bcdof_top, bcval_top   = Vector{Int64}(), Vector{Float64}()
    # Cylinder
    bcdof_cyl, bcval_cyl   = setBCXY_X(Δ/nloadsteps, dh, n_cyl)
    # Collect
    bcdofs = [bcdof_left; bcdof_top; bcdof_cyl]
    bcvals = [bcval_left; bcval_top; bcval_cyl]
    #
    ϵᵢⱼₖ = sortperm(bcdofs)
    global bcdofs = bcdofs[ϵᵢⱼₖ]
    global bcvals = bcvals[ϵᵢⱼₖ]

    # - For Linear solver..
    global pdofs = bcdofs
    global fdofs = setdiff(1:dh.ndofs.x, pdofs)

    bcval₀ = bcvals
    global β = 1.0
    #for loadstep ∈ 1 : nloadsteps
    ##
    loadstep = 0
    while loadstep < nloadsteps
        loadstep += 1
        #global ε = ε * 1.1
        ##
        res = res .* 0
        bcvals = bcval₀
        residual = 0 * residual
        iter = 0
        fill!(Δa, 0.0)
        print("\n", "Starting equilibrium iteration at loadstep: ", loadstep, "\n")
        #global ε = ε₀
        a_old = a
        # # # # # # # # # #
        # Newton solve.   #
        # # # # # # # # # #

        #@show β
        while residual > TOL || iter < 2
            iter += 1
            if iter % 20 == 0 || norm(res) > 1e3 && β > 1 / 8
                a = a_old
                bcvals = bcval₀
                global β = β * 0.5
                #Δ_remaining = (Δ * nloadsteps - β * Δ - loadstep * Δ) / nloadsteps
                remaining_steps = nloadsteps - loadstep
                nloadsteps = loadstep + 2remaining_steps + (1 / β - 1)
                bcvals = bcvals ./ 2 #
                bcval₀ = bcvals
                fill!(Δa, 0.0)
                println("Penalty paremeter and updated: $ε, and step length $β ")
            end

            #a += β * Δa
            a += Δa
            assemGlobal!(K, Fᵢₙₜ, rc, dh, mp, t, a, coord, enod, ε)
            solveq!(Δa, K, -Fᵢₙₜ, bcdofs, bcvals)
            bcvals = 0 * bcvals
            res = Fᵢₙₜ - Fₑₓₜ
            res[bcdofs] = 0 * res[bcdofs]
            residual = norm(res, 2)
            @printf "Iteration: %i | Residual: %.4e | Δ: %.4f \n" iter residual a[bcdof_cyl[1]]
        end
        if loadstep < 40 && iter < 20
            σx, σy = StressExtract(dh, a, mp)
            vtk_grid("results/contact" * string(loadstep), dh) do vtkfile
                vtk_point_data(vtkfile, dh, a) # displacement field
                vtk_point_data(vtkfile, σx, "σx")
                vtk_point_data(vtkfile, σy, "σy")
            end
        end
        Fₑₓₜ[bcdofs] = -Fᵢₙₜ[bcdofs]
    end
    τ_c = ExtractContactTraction(a, ε, coord)
    return a, dh, Fₑₓₜ, Fᵢₙₜ, K, τ_c
end

function solver_C_U(dh, coord, Δ, nloadsteps)

    # ---------- #
    # Set params # // Kanske som input till solver???
    # ---------- # // definiera mp här? och kanske ε ? iofs snyggare utanför!
    t = 1.0

    # Define material parameters
    #mp = [175 80.769230769230759]

    # ------------- #
    # Init-stuff    #
    # ------------- #
    imax = 200
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
    global rc = zeros(dh.ndofs.x)
    global Fₑₓₜ = zeros(dh.ndofs.x)
    global a = zeros(dh.ndofs.x)
    global Δa = zeros(dh.ndofs.x)
    global res = zeros(dh.ndofs.x)
    global K = create_sparsity_pattern(dh)

    # ------------------- #
    # Boundary conditions #
    # ------------------- #
    # L structure
    bcdof_bot, bcval_bot = setBCXY_both(0.0, dh, n_bot)
    bcdof_top, bcval_top = setBCXY(-Δ/nloadsteps, dh, n_top)
    # Collect
    bcdofs = [bcdof_bot; bcdof_top]
    bcvals = [bcval_bot; bcval_top]
    #
    ϵᵢⱼₖ = sortperm(bcdofs)
    global bcdofs = bcdofs[ϵᵢⱼₖ]
    global bcvals = bcvals[ϵᵢⱼₖ]

    # - For Linear solver..
    global pdofs = bcdofs
    global fdofs = setdiff(1:dh.ndofs.x, pdofs)

    bcval₀ = bcvals
    global β = 1.0
    #for loadstep ∈ 1 : nloadsteps
    ##
    loadstep = 0
    while loadstep < nloadsteps
        loadstep += 1
        #global ε = ε * 1.1
        ##
        res = res .* 0
        bcvals = bcval₀
        residual = 0 * residual
        iter = 0
        fill!(Δa, 0.0)
        print("\n", "Starting equilibrium iteration at loadstep: ", loadstep, "\n")
        #global ε = ε₀
        a_old = a
        # # # # # # # # # #
        # Newton solve.   #
        # # # # # # # # # #

        #@show β
        while residual > TOL || iter < 2
            iter += 1
            if iter % 20 == 0 || norm(res) > 1e3 && β > 1 / 8
                a = a_old
                bcvals = bcval₀
                global β = β * 0.5
                Δ_remaining = (Δ * nloadsteps - β * Δ - loadstep * Δ) / nloadsteps
                remaining_steps = nloadsteps - loadstep
                nloadsteps = loadstep + 2remaining_steps + (1 / β - 1)
                bcvals = bcvals ./ 2 #
                bcval₀ = bcvals
                fill!(Δa, 0.0)
                println("Penalty paremeter and updated: $ε, and step length $β ")
            end

            #a += β * Δa
            a += Δa
            assemGlobal!(K, Fᵢₙₜ, rc, dh, mp, t, a, coord, enod, ε)
            solveq!(Δa, K, -Fᵢₙₜ, bcdofs, bcvals)
            bcvals = 0 * bcvals
            res = Fᵢₙₜ - Fₑₓₜ
            res[bcdofs] = 0 * res[bcdofs]
            residual = norm(res, 2)
            @printf "Iteration: %i | Residual: %.4e | Δ: %.4f \n" iter residual a[bcdof_top[1]]
        end
        if loadstep < 40 && iter < 20
            σx, σy = StressExtract(dh, a, mp)
            vtk_grid("results/contact" * string(loadstep), dh) do vtkfile
                vtk_point_data(vtkfile, dh, a) # displacement field
                vtk_point_data(vtkfile, σx, "σx")
                vtk_point_data(vtkfile, σy, "σy")
            end
        end
        Fₑₓₜ[bcdofs] = -Fᵢₙₜ[bcdofs]
    end
    τ_c = ExtractContactTraction(a, ε, coord)
    return a, dh, Fₑₓₜ, Fᵢₙₜ, K, τ_c
end


function fictitious_solver_hook(d, dh0, coord₀, nloadsteps)
    # allt överflödigt bör vid tillfälle flyttas utanför
    # lösare till ett "init-liknande script så att huvudsaklig kod hålls ren
    TOL = 1e-10
    residual = 0.0
    iter = 1
    global λ = 0
    ndof = size(coord₀, 1) * 2
    nelm = size(enod, 1)
    t = 1.0

    #  ----- #
    # Init   #
    #  ----- #
    global Kψ  = create_sparsity_pattern(dh0)
    global Ψ   = zeros(dh0.ndofs.x)
    global FΨ  = zeros(dh0.ndofs.x)
    global Ψ   = zeros(dh0.ndofs.x)
    global ΔΨ  = zeros(dh0.ndofs.x)
    global res = zeros(dh0.ndofs.x)

    global bcdof_o2 = bcdofs_opt
    global bcval_o2 = bcdofs_opt .* 0.0
    global pdofs = bcdofs_opt
    global fdofs = setdiff(1:ndof, pdofs)

    bcval₀_o2 = bcval_opt
    Δλ = (1.0 / nloadsteps)
    loadstep = 0

    while loadstep < nloadsteps
        loadstep += 1
        res = res .* 0
        bcval_opt = bcval₀_o2
        residual = 0 * residual
        iter = 0
        global λ += Δλ #* loadstep
        fill!(ΔΨ, 0.0)
        print("\n", "Starting equilibrium iteration at loadstep: ", loadstep, "\n")
        while residual > TOL || iter < 2
            iter += 1
            Ψ += ΔΨ
            assemGlobal!(Kψ, FΨ, dh0, mp₀, t, Ψ, coord₀, enod, λ, d, Γ_robin)
            solveq!(ΔΨ, Kψ, -FΨ, bcdofs_opt, bcval_opt)

            bcval_opt = bcval_opt .* 0
            res = FΨ #- Fₑₓₜ
            res[bcdofs_opt] = res[bcdofs_opt] .* 0
            residual = norm(res, 2)
            Ψ[bcdofs_opt] .= 0.0
            if iter < 20
                postprocess_opt(res, dh0, "results/fictres" * string(iter))
                postprocess_opt(Ψ, dh0, "results/fictitious_iter" * string(iter))
            end
            @printf "Iteration: %i | Residual: %.4e | λ: %.4f \n" iter residual λ
        end
    end
    return Ψ, dh0, Kψ, FΨ, λ
end
#
function solver_Lab(dh, coord, Δ, nloadsteps)
    # ---------- #
    # Set params #
    # ---------- #
    t = 1.0
    # ------------- #
    # Init-stuff    #
    # ------------- #
    imax     = 200
    TOL      = 1e-8
    residual = 0.0
    iter     = 1
    # --------------- #
    #  Init matrices  #
    # --------------- #
    global Fᵢₙₜ = zeros(dh.ndofs.x)
    global rc = zeros(dh.ndofs.x)
    global Fₑₓₜ = zeros(dh.ndofs.x)
    global a = zeros(dh.ndofs.x)
    global Δa = zeros(dh.ndofs.x)
    global res = zeros(dh.ndofs.x)
    global K = create_sparsity_pattern(dh)
    # ------------------- #
    # Boundary conditions #
    # ------------------- #
    bcdof_bot, bcval_bot = setBCY(0.0, dh, n_bot)
    bcdof_top, bcval_top = setBCY(Δ / nloadsteps, dh, n_top)
    bcdof_right, bcval_right = setBCX(0.0, dh, n_sym)

    bcdofs = [bcdof_bot; bcdof_top; bcdof_right]
    bcvals = [bcval_bot; bcval_top; bcval_right]
    ϵᵢⱼₖ  = sortperm(bcdofs)
    global bcdofs = bcdofs[ϵᵢⱼₖ]
    global bcvals = bcvals[ϵᵢⱼₖ]
    # - - - - - - - - - #
    # For Linear solver #
    # - - - - - - - - - #
    global pdofs = bcdofs
    global fdofs = setdiff(1:dh.ndofs.x, pdofs)
    bcval₀   = bcvals

    loadstep = 0
    global β = 1.
    while loadstep < nloadsteps
        loadstep += 1
        #τ         = [0.0; 1e1]* loadstep/nloadsteps
        # global ε = ε * 1.2
        res = res .* 0
        bcvals = bcval₀
        residual = 0 * residual
        iter = 0
        fill!(Δa, 0.0)
        print("\n", "Starting equilibrium iteration at loadstep: ", loadstep, "\n")
        a_old = a
        # # # # # # # # # #
        # Newton solve.   #
        # # # # # # # # # #
        while residual > TOL || iter < 2
            iter += 1
            if iter % 10 == 0 || norm(res) > 1e3
                    a = a_old
                    bcvals = bcval₀
                    if β > 1/8
                        global β = β * 0.5
                        Δ_remaining = (Δ*nloadsteps - β * Δ - loadstep * Δ)/nloadsteps
                        remaining_steps = nloadsteps - loadstep
                        nloadsteps = loadstep + 2remaining_steps + (1 / β - 1)
                        bcvals = bcvals ./2 #
                        bcval₀= bcvals
                    end
                    fill!(Δa, 0.0)
                    println("Step length $β ")
            end
            a += Δa
            assemGlobal!(K, Fᵢₙₜ, dh, t, a, coord, enod, ε, mp₁, mp₂)
            #@show Fᵢₙₜ[contact_dofs]
            #assemGlobal!(K, Fᵢₙₜ, dh, t, a, coord, enod, ε, mp₁, mp₂, τ)
            solveq!(Δa, K, -Fᵢₙₜ, bcdofs, bcvals)
            bcvals = 0 * bcvals
            res = Fᵢₙₜ - Fₑₓₜ
            res[bcdofs] = 0 * res[bcdofs]
            residual = norm(res, 2)
            @printf "Iteration: %i | Residual: %.4e | Δ: %.4f \n" iter residual a[bcdof_top[1]]
        #
        end
        if loadstep < 40 && iter < 20
            σx, σy = StressExtract(dh, a, mp₁) # måste ändra så att vi kör med mp₁ & mp₂
            vtk_grid("results/🍌-contact" * string(loadstep), dh) do vtkfile
                vtk_point_data(vtkfile, dh, a )
                vtk_point_data(vtkfile, σx, "σx")
                vtk_point_data(vtkfile, σy, "σy")
            end
        end
        Fₑₓₜ[bcdofs] = -Fᵢₙₜ[bcdofs]
    end
    # X_c,tract = plotTraction()
    # if length(tract) > 0
    #     p5 = plot(X_c, tract, label="λ" , marker=4, lc=:tomato, mc=:tomato, grid=false, legend=:outerleft, ylims = (0, 1.2*maximum(tract)) )
    #     display(p5)
    # end
    return a, dh, Fₑₓₜ, Fᵢₙₜ, K
end

function fictitious_solver_with_contact_lab(d, dh0, coord₀, nloadsteps)
    TOL = 1e-10
    residual = 0.0
    iter = 1
    global λ = 0
    t = 1.0
    #  ----- #
    # Init   #
    #  ----- #
    global Kψ  = create_sparsity_pattern(dh0)
    global Ψ   = zeros(dh0.ndofs.x)
    global FΨ  = zeros(dh0.ndofs.x)
    global Ψ   = zeros(dh0.ndofs.x)
    global ΔΨ  = zeros(dh0.ndofs.x)
    global res = zeros(dh0.ndofs.x)
    global bcdof_o2 = bcdofs_opt
    global bcval_o2 = bcdofs_opt .* 0.0
    global pdofs = bcdofs_opt
    global fdofs = setdiff(1:dh0.ndofs.x, pdofs)
    bcval₀_o2 = bcval_opt
    Δλ = (1.0 / nloadsteps)
    loadstep = 0
    while loadstep < nloadsteps
        loadstep += 1
        res = res .* 0
        bcval_opt = bcval₀_o2
        residual = 0 * residual
        iter = 0
        global λ += Δλ #* loadstep
        fill!(ΔΨ, 0.0)
        print("\n", "Starting equilibrium iteration at loadstep: ", loadstep, "\n")
        Ψ_old = Ψ
        # # # # # # # # # #
        # Newton solve.  #
        # # # # # # # # # #
        while residual > TOL || iter < 2
            iter += 1
            if iter % 20 == 0 || norm(res) > 1e2 #&& Δλ > 1/16
                Ψ = Ψ_old
                if Δλ > 0.1 * 1 / 64
                    global λ -= Δλ #* loadstep
                    Δλ = Δλ / 2
                    global λ += Δλ  #* loadstep
                    remaining_steps = nloadsteps - loadstep
                    nloadsteps = loadstep + round((1 - λ) / Δλ)
                end
                fill!(ΔΨ, 0.0)
                println("Step length updated: $Δλ, penalty parameter: $μ")
            end
            Ψ += ΔΨ
            assemGlobal!(Kψ, FΨ, dh0, mp₀, t, Ψ, coord₀, enod, λ, d, Γ_robin, μ)
            solveq!(ΔΨ, Kψ, -FΨ, bcdofs_opt, bcval_opt)
            #
            bcval_opt = bcval_opt .* 0
            res       = FΨ #- Fₑₓₜ
            res[bcdofs_opt] = res[bcdofs_opt] .* 0
            residual        = norm(res, 2)
            Ψ[bcdofs_opt]  .= 0.0
            @printf "Iteration: %i | Residual: %.4e | λ: %.4f \n" iter residual λ
            if loadstep < 40 && iter < 20
                postprocess_opt(Ψ, dh0, "results/fictitious_t2" * string(loadstep))
            end
        end
    end
    return Ψ, dh0, Kψ, FΨ, λ
end

function solver_arc(dh, coord, Δ, nloadsteps)

    # ---------- #
    # Set params #
    # ---------- #
    t = 1.0

    # Define material parameters
    mp = [175 80.769230769230759]

    # ------------- #
    # Init-stuff    #
    # ------------- #
    imax = 200
    TOL = 1e-8
    residual = 0.0
    iter = 1

    # ------ #
    #  Init  #
    # ------ #
    global Fᵢₙₜ = zeros(dh.ndofs.x)
    global rc = zeros(dh.ndofs.x)
    global Fₑₓₜ = zeros(dh.ndofs.x)
    global a = zeros(dh.ndofs.x)
    global Δa = zeros(dh.ndofs.x)
    global res = zeros(dh.ndofs.x)
    global K = create_sparsity_pattern(dh)

    # ------------------- #
    # Boundary conditions #
    # ------------------- #
    bcdof_left, bcval_left = setBCXY_X(0.0, dh, n_left)
    bcdof_right, bcval_right = setBCXY_Y(Δ / nloadsteps, dh, n_right)

    bcdof_bot, bcval_bot = Vector{Int64}(), Vector{Float64}()
    bcdof_top, bcval_top = Vector{Int64}(), Vector{Float64}()

    bcdofs = [bcdof_left; bcdof_right; bcdof_bot; bcdof_top]
    bcvals = [bcval_left; bcval_right; bcval_bot; bcval_top]
    ## =#
    ϵᵢⱼₖ = sortperm(bcdofs)
    global bcdofs = bcdofs[ϵᵢⱼₖ]
    global bcvals = bcvals[ϵᵢⱼₖ]

    # - For Linear solver..
    global pdofs = bcdofs
    global fdofs = setdiff(1:dh.ndofs.x, pdofs)

    bcval₀ = bcvals
    global β = 1.0
    #for loadstep ∈ 1 : nloadsteps
    ##
    loadstep = 0
    while loadstep < nloadsteps
        loadstep += 1
        #global ε = ε * 1.1
        ##
        res = res .* 0
        bcvals = bcval₀
        residual = 0 * residual
        iter = 0
        fill!(Δa, 0.0)
        print("\n", "Starting equilibrium iteration at loadstep: ", loadstep, "\n")
        # # # # # # # # # #
        # Newton solve.   #
        # # # # # # # # # #

        while residual > TOL || iter < 2
            iter += 1
            a += Δa
            assemGlobal!(K, Fᵢₙₜ, dh, mp, t, a, coord, enod)
            solveq!(Δa, K, -Fᵢₙₜ, bcdofs, bcvals)
            bcvals = 0 * bcvals
            res = Fᵢₙₜ - Fₑₓₜ
            res[bcdofs] = 0 * res[bcdofs]
            residual = norm(res, 2)
            @printf "Iteration: %i | Residual: %.4e | Δ: %.4f \n" iter residual a[bcdof_right[1]]
        end
        if loadstep < 40 && iter < 20
            σx, σy = StressExtract(dh, a, mp)
            vtk_grid("results/🍌-contact" * string(loadstep), dh) do vtkfile
                vtk_point_data(vtkfile, dh, a) # displacement field
                vtk_point_data(vtkfile, σx, "σx")
                vtk_point_data(vtkfile, σy, "σy")
            end
        end
        Fₑₓₜ[bcdofs] = -Fᵢₙₜ[bcdofs]
        if loadstep < 11
            a_hist[:, loadstep] = a
        end
    end
    return a, dh, Fₑₓₜ, Fᵢₙₜ, K, a_hist
end

function solver_hook(dh, coord, Δ, nloadsteps)

    # ---------- #
    # Set params #
    # ---------- #
    t = 1.0

    # Define material parameters
    mp = [175 80.769230769230759]

    # ------------- #
    # Init-stuff    #
    # ------------- #
    imax = 200
    TOL = 1e-8
    residual = 0.0
    iter = 1

    # ------ #
    #  Init  #
    # ------ #
    global Fᵢₙₜ = zeros(dh.ndofs.x)
    global rc = zeros(dh.ndofs.x)
    global Fₑₓₜ = zeros(dh.ndofs.x)
    global a = zeros(dh.ndofs.x)
    global Δa = zeros(dh.ndofs.x)
    global res = zeros(dh.ndofs.x)
    global K = create_sparsity_pattern(dh)

    # ------------------- #
    # Boundary conditions #
    # ------------------- #
    bcdof_left, bcval_left = setBCXY_X(0.0, dh, n_left)
    #bcdof_left, bcval_left     = setBCXY_X( -Δ / nloadsteps, dh, n_left)
    bcdof_right, bcval_right = setBCXY_X(Δ / nloadsteps, dh, n_right)
    #bcdof_bot, bcval_bot = setBCY(0.0, dh, n_bot)
    #bcdof_top, bcval_top = setBCY(0.0, dh, n_top)

    bcdof_bot, bcval_bot = Vector{Int64}(), Vector{Float64}()
    bcdof_top, bcval_top = Vector{Int64}(), Vector{Float64}()

    bcdofs = [bcdof_left; bcdof_right; bcdof_bot; bcdof_top]
    bcvals = [bcval_left; bcval_right; bcval_bot; bcval_top]
    ## =#
    ϵᵢⱼₖ = sortperm(bcdofs)
    global bcdofs = bcdofs[ϵᵢⱼₖ]
    global bcvals = bcvals[ϵᵢⱼₖ]

    # - For Linear solver..
    global pdofs = bcdofs
    global fdofs = setdiff(1:dh.ndofs.x, pdofs)

    bcval₀ = bcvals
    global β = 1.0
    #for loadstep ∈ 1 : nloadsteps
    ##
    loadstep = 0
    while loadstep < nloadsteps
        loadstep += 1
        #global ε = ε * 1.1
        ##
        res = res .* 0
        bcvals = bcval₀
        residual = 0 * residual
        iter = 0
        fill!(Δa, 0.0)
        print("\n", "Starting equilibrium iteration at loadstep: ", loadstep, "\n")
        # # # # # # # # # #
        # Newton solve.   #
        # # # # # # # # # #

        while residual > TOL || iter < 2
            iter += 1
            a += Δa
            assemGlobal!(K, Fᵢₙₜ, dh, mp, t, a, coord, enod)
            solveq!(Δa, K, -Fᵢₙₜ, bcdofs, bcvals)
            bcvals = 0 * bcvals
            res = Fᵢₙₜ - Fₑₓₜ
            res[bcdofs] = 0 * res[bcdofs]
            residual = norm(res, 2)
            @printf "Iteration: %i | Residual: %.4e | Δ: %.4f \n" iter residual a[bcdof_right[1]]
        end
        if loadstep < 40 && iter < 20
            σx, σy = StressExtract(dh, a, mp)
            vtk_grid("results/🍌-contact" * string(loadstep), dh) do vtkfile
                vtk_point_data(vtkfile, dh, a) # displacement field
                vtk_point_data(vtkfile, σx, "σx")
                vtk_point_data(vtkfile, σy, "σy")
            end
        end
        Fₑₓₜ[bcdofs] = -Fᵢₙₜ[bcdofs]
        if loadstep < 11
            a_hist[:, loadstep] = a
        end
    end
    return a, dh, Fₑₓₜ, Fᵢₙₜ, K, a_hist
end
