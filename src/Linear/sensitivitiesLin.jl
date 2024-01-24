# Sensitivity of fictious material with respect to design parameter 'd'
function drψ(dr_dd,dh,Ψ,λ,d,Γ_robin,coord₀)
    assembler = start_assemble(dr_dd)
    ie = 0
    for cell in CellIterator(dh)
        ie += 1
        cell_dofs= celldofs(cell)
        dfe = zeros(6,6) ## 12x1 eller 12x12 ??
        for face in 1:nfaces(cell)
            if (cellid(cell), face) in Γ_robin
                face_nods = [ Ferrite.facedof_indices(ip)[face][1]; Ferrite.facedof_indices(ip)[face][2] ]
                face_dofs = [ face_nods[1]*2-1; face_nods[1]*2; face_nods[2]*2-1; face_nods[2]*2 ]
                X         = coord₀[ enod[ie][face_nods.+1] ,: ]
                dfe[face_dofs,face_dofs],_ = Robin(X,Ψ[cell_dofs[face_dofs]],d[cell_dofs[face_dofs]],λ)
            end
        end
        # Dampened
        #assemble!(assembler, cell_dofs, -0.1dfe)
        # Regular
        assemble!(assembler, cell_dofs, -dfe)
    end
    return dr_dd
end

# Shape sensitivity of equillibrium residual
function drᵤ_dx(dr,dh,mp,t,a,coord,enod,τ, Γt)
    assembler = start_assemble(dr)
    ie = 0
    drₑ = zeros(6,6)
    for cell in CellIterator(dh)
        ie += 1
        cell_dofs= celldofs(cell)
        drₑ      = assem_dr(coord[enod[ie][2:end],:],a[cell_dofs],mp,t)
        dre      = zeros(6,6)
        for face in 1:nfaces(cell)
            if (cellid(cell), face) in Γt
                face_nods                 = [ Ferrite.facedof_indices(ip)[face][1]; Ferrite.facedof_indices(ip)[face][2] ]
                face_dofs                 = [ face_nods[1]*2-1; face_nods[1]*2; face_nods[2]*2-1; face_nods[2]*2 ]
                X                         = coord[ enod[ie][face_nods.+1] ,: ]
                dre[face_dofs,face_dofs]  = dTractionLoad(X,τ)
            end
        end
        assemble!(assembler, cell_dofs, drₑ+dre)
    end
    return dr
end

# Shape sensitivity of equillibrium residual
function drᵤ_dx(dr, dh, mp, t, a, coord, enod)
    assembler = start_assemble(dr)
    ie = 0
    drₑ = zeros(6, 6)
    for cell in CellIterator(dh)
        ie += 1
        cell_dofs = celldofs(cell)
        drₑ = assem_dr(coord[enod[ie][2:end], :], a[cell_dofs], mp, t)
        dre = zeros(6, 6)
        assemble!(assembler, cell_dofs, drₑ + dre)
    end
    return dr
end


function drᵤ_dx_c(∂rᵤ_∂x, dh, mp, t, a, coord, enod, ε)
    assembler = start_assemble(∂rᵤ_∂x)
    ie = 0
    drₑ = zeros(6, 6)
    for cell in CellIterator(dh)
        ie += 1
        cell_dofs = celldofs(cell)
        drₑ = assem_dr(coord[enod[ie][2:end], :], a[cell_dofs], mp, t)
        dre = zeros(6, 6)
        assemble!(assembler, cell_dofs, drₑ + dre)
    end

    X_ordered = getXinDofOrder(dh, X, coord)
    drc       = ForwardDiff.jacobian(x -> contact_residual_ordered(x, a, ε), X_ordered)
    # dr[contact_dofs,contact_dofs] += drc[contact_dofs,contact_dofs]
    ∂rᵤ_∂x   -= drc

    return ∂rᵤ_∂x
end
# Shape sensitivity of external forces
function dFext_dx(dF,dh,mp,t,a,coord,enod,τ, Γt)
    assembler = start_assemble(dF)
    ie = 0
    for cell in CellIterator(dh)
        ie += 1
        cell_dofs= celldofs(cell)
        dre      = zeros(6,6)
        for face in 1:nfaces(cell)
            if (cellid(cell), face) in Γt
                face_nods                 = [ Ferrite.facedof_indices(ip)[face][1]; Ferrite.facedof_indices(ip)[face][2] ]
                face_dofs                 = [ face_nods[1]*2-1; face_nods[1]*2; face_nods[2]*2-1; face_nods[2]*2 ]
                X                         = coord[ enod[ie][face_nods.+1] ,: ]
                dre[face_dofs,face_dofs]  = dTractionLoad(X,τ)
            end
        end
        assemble!(assembler, cell_dofs, -dre)
    end
    return dF
end

# Shape sensitivity of equillibrium residual including contact
function drᵤ_dx_c(∂rᵤ_∂x, dh, mp, t, a, coord, enod, ε)
    assembler = start_assemble(∂rᵤ_∂x)
    ie = 0
    drₑ = zeros(6, 6)
    for cell in CellIterator(dh)
        ie += 1
        cell_dofs = celldofs(cell)
        drₑ = assem_dr(coord[enod[ie][2:end], :], a[cell_dofs], mp, t)
        dre = zeros(6, 6)
        assemble!(assembler, cell_dofs, drₑ + dre)
    end

    X_ordered = getXinDofOrder(dh, X, coord)
    drc       = ForwardDiff.jacobian(x -> contact_residual_ordered(x, a, ε), X_ordered)
    # dr[contact_dofs,contact_dofs] += drc[contact_dofs,contact_dofs]
    ∂rᵤ_∂x   -= drc

    return ∂rᵤ_∂x
end

# Shape sensitivity of equillibrium residual including contact for two materials
function drᵤ_dx_c(∂rᵤ_∂x, dh, t, a, coord, enod, ε, mp₁, mp₂)
    assembler = start_assemble(∂rᵤ_∂x)
    ie = 0
    drₑ = zeros(6, 6)
    for cell in CellIterator(dh)
        ie += 1
        cell_dofs = celldofs(cell)
        if ie ∈ dh.grid.cellsets["top mesh"]
            mp = mp₁
        else
            mp = mp₂
        end
        drₑ = assem_dr(coord[enod[ie][2:end], :], a[cell_dofs], mp, t)
        dre = zeros(6, 6)
        assemble!(assembler, cell_dofs, drₑ + dre)
    end

    X_ordered = getXinDofOrder(dh, X, coord)
    drc = ForwardDiff.jacobian(x -> contact_residual_ordered(x, a, ε), X_ordered)
    # dr[contact_dofs,contact_dofs] += drc[contact_dofs,contact_dofs]
    ∂rᵤ_∂x -= drc

    return ∂rᵤ_∂x
end

# Shape sensitivity of fictitious residual including contact
function drΨ_dx_c(∂rΨ_∂x, dh0, mp₀, t, Ψ, coord₀, enod, λ, d, Γ_robin, μ)
    assembler = start_assemble(∂rΨ_∂x)
    ie = 0
    for cell in CellIterator(dh0)
        ie += 1
        cell_dofs = celldofs(cell)
        kₑ = zeros(6, 6)
        fₑ = zeros(6)
        kₑ, _ = assemElem(coord₀[enod[ie][2:end], :], Ψ[cell_dofs], mp₀, t)
        ke = zeros(6, 6)
        fe = zeros(6)
        for face in 1:nfaces(cell)
            if (cellid(cell), face) in Γ_robin
                face_nods = [Ferrite.facedof_indices(ip)[face][1]; Ferrite.facedof_indices(ip)[face][2]]
                face_dofs = [face_nods[1] * 2 - 1; face_nods[1] * 2; face_nods[2] * 2 - 1; face_nods[2] * 2]
                Xf = coord₀[enod[ie][face_nods.+1], :]
                ke[face_dofs, face_dofs], _ = Robin(Xf, Ψ[cell_dofs[face_dofs]], d[cell_dofs[face_dofs]], λ)
            end
        end
        assemble!(assembler, cell_dofs, kₑ + ke)
    end
    # Contact
    X_ordered = getXfromCoord(coord₀)
    #X_ordered = getXinDofOrder(dh0, X, coord₀)
    Kc = ForwardDiff.jacobian(x -> contact_residual(x, Ψ, μ), X_ordered)
    ∂rΨ_∂x[contact_dofs, contact_dofs] -= Kc[contact_dofs, contact_dofs]

    #X_ordered = getXinDofOrder(dh0, X, coord₀)
    #drc       = ForwardDiff.jacobian(x -> contact_residual_ordered(x, Ψ, μ), X_ordered)
    #∂rΨ_∂x   -= drc
    return ∂rΨ_∂x
end

# Objective function for p-norm of contact pressure
function contact_pnorm(X::AbstractVector{T1}, a::AbstractVector{T2}, ε, p) where {T1,T2}

   # Order displacements according to nodes and not dofs
   a_ordered = getDisplacementsOrdered(dh, a)

   # Scaling
   κ = gap_scaling(X)

   # convert X to Real for compatibility with ForwardDiff
   #X_float = real.(X)  + real.(a_ordered) # a ska vara sorterad på samma sätt som X, detta måste fixas!!!!!!!!!
   X_float = real.(X + a_ordered) # a ska vara sorterad på samma sätt som X, detta måste fixas!!!!!!!!!

   # Extract the coordinate vector (nbr_nodes x 2 )
   coordu = getCoordfromX(X_float)

   # Create dictionaries that are needed for the Mortar2D package
   elements, element_types, slave_elements, slave_element_ids, master_element_ids, coords = create_contact_list(dh, Γs, Γm, coordu)

   # Compute nodal normals
   normals = Mortar2D.calculate_normals(elements, element_types, coords)

   # Compute the projected gap function
   g = gap_function(X_float)

   #println("norm(scaling): ",norm(κ))
   # Assemble D and M matrices and the slave and master dofs corresponding to the mortar segmentation
   slave_dofs, master_dofs, D, M = Mortar2D.calculate_mortar_assembly(elements, element_types, coords, slave_element_ids, master_element_ids)

   # ---------------- #
   # (∑ₐ λₐᵖ )^(1/p)  #
   # ---------------- #

   # Loop over master side dofs
    g₀ = 0.0
   for (i, A) in enumerate(slave_dofs)
        λ_A = penalty(g[i, :] ⋅ normals[slave_dofs[i]], ε)
        g₀ += (λ_A  * (1 / κ[i]))^p
   end

   # ---------------------------------- #
   # ∫ᵧ g 𝛅λ dγ = 0 for penalty methods  #
   # ---------------------------------- #

   return (g₀)^(1/p)
end

# Objective function for p-norm of contact pressure | ordered
function contact_pnorm_ordered(X::AbstractVector{T1}, a::AbstractVector{T2}, ε,p) where {T1,T2}

    # Order  X
    X_ordered = getX_from_Dof_To_Node_order(dh, X)

    #X_ordered = getXinDofOrder(dh, X, coord)

    r_c = contact_pnorm(X_ordered, a, ε,p)

    return r_c
end

# Objective function for p-norm of contact pressure | Other version
function contact_pnorm_s(X::AbstractVector{T1}, a::AbstractVector{T2}, ε, p) where {T1,T2}

    # Order displacements according to nodes and not dofs
    a_ordered = getDisplacementsOrdered(dh, a)

    # Scaling
    κ = gap_scaling(X)

    # convert X to Real for compatibility with ForwardDiff
    #X_float = real.(X)  + real.(a_ordered) # a ska vara sorterad på samma sätt som X, detta måste fixas!!!!!!!!!
    X_float = real.(X + a_ordered) # a ska vara sorterad på samma sätt som X, detta måste fixas!!!!!!!!!

    # Extract the coordinate vector (nbr_nodes x 2 )
    coordu = getCoordfromX(X_float)

    # Create dictionaries that are needed for the Mortar2D package
    elements, element_types, slave_elements, slave_element_ids, master_element_ids, coords = create_contact_list(dh, Γs, Γm, coordu)

    # Compute nodal normals
    normals = Mortar2D.calculate_normals(elements, element_types, coords)

    # Compute the projected gap function
    g = gap_function(X_float)

    #println("norm(scaling): ",norm(κ))
    # Assemble D and M matrices and the slave and master dofs corresponding to the mortar segmentation
    slave_dofs, master_dofs, D, M = Mortar2D.calculate_mortar_assembly(elements, element_types, coords, slave_element_ids, master_element_ids)

    # ---------------- #
    # (∑ₐ λₐᵖ )^(1/p)  #
    # ---------------- #

    # Loop over master side dofs
    g₀ = 0.0
    Ω  = 0.0
    p_mean = []
    for (i, A) in enumerate(slave_dofs)
        λ_A = penalty(g[i, :] ⋅ normals[slave_dofs[i]], ε)
        #g₀ += (λ_A * (1 / κ[i]))
        if λ_A != 0
            append!(p_mean, λ_A * (1 / κ[i]))
        end
    end

    #g₀ = (g₀/length(p_mean) - mean(p_mean)^2)
    g₀ = var(p_mean)

    # ---------------------------------- #
    # ∫ᵧ g 𝛅λ dγ = 0 for penalty methods #
    # ---------------------------------- #
    return g₀
end

# Objective function for p-norm of contact pressure | Other version |
function contact_pnorm_ordered_s(X::AbstractVector{T1}, a::AbstractVector{T2}, ε, p) where {T1,T2}

    # Order  X
    X_ordered = getX_from_Dof_To_Node_order(dh, X)

    #X_ordered = getXinDofOrder(dh, X, coord)

    r_c = contact_pnorm_s(X_ordered, a, ε, p)

    return r_c
end

# Compute total contact area
function contact_area(X::AbstractVector{T1}, a::AbstractVector{T2}) where {T1,T2}

    # Order displacements according to nodes and not dofs
    a_ordered = getDisplacementsOrdered(dh, a)

    # Scaling
    κ = gap_scaling(X)

    # convert X to Real for compatibility with ForwardDiff
    #X_float = real.(X)  + real.(a_ordered) # a ska vara sorterad på samma sätt som X, detta måste fixas!!!!!!!!!
    X_float = real.(X + a_ordered) # a ska vara sorterad på samma sätt som X, detta måste fixas!!!!!!!!!

    # Extract the coordinate vector (nbr_nodes x 2 )
    coordu = getCoordfromX(X_float)

    # Create dictionaries that are needed for the Mortar2D package
    elements, element_types, slave_elements, slave_element_ids, master_element_ids, coords = create_contact_list(dh, Γs, Γm, coordu)

    # Compute nodal normals
    normals = Mortar2D.calculate_normals(elements, element_types, coords)

    # Compute the projected gap function
    g = gap_function(X_float)

    #println("norm(scaling): ",norm(κ))
    # Assemble D and M matrices and the slave and master dofs corresponding to the mortar segmentation
    slave_dofs, master_dofs, D, M = Mortar2D.calculate_mortar_assembly(elements, element_types, coords, slave_element_ids, master_element_ids)

    # ---------------- #
    # (∑ₐ λₐᵖ )^(1/p)  #
    # ---------------- #

    # Loop over master side dofs
    Ω = 0.0
    for (i, A) in enumerate(slave_dofs)
        Ω += (penalty(g[i, :] ⋅ normals[slave_dofs[i]], 1.0)) / κ[i]
    end

    return Ω
end

# Compute total contact area | input-dofs ordered to suit AD
function contact_area_ordered(X::AbstractVector{T1}, a::AbstractVector{T2}) where {T1,T2}

    # Order  X
    X_ordered = getX_from_Dof_To_Node_order(dh, X)

    #X_ordered = getXinDofOrder(dh, X, coord)

    r_c = contact_pnorm_s(X_ordered, a, ε, p)

    return r_c
end
