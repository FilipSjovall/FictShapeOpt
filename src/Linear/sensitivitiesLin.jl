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
        assemble!(assembler, cell_dofs, -dfe)
    end
    return dr_dd
end

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

function drᵤ_dx_c(dr, dh, mp, t, a, coord, enod, ε)
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

    X_ordered = getXinDofOrder(dh, X, coord)
    drc       = ForwardDiff.jacobian(x -> contact_residual_ordered(x, a, ε), X_ordered)

    # dr[contact_dofs,contact_dofs] += drc[contact_dofs,contact_dofs]
    dr -= drc

    return dr
end


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

   # Initialize the nodal gap vector.
   gₙ = zeros(eltype(X_float), length(slave_dofs))

   # Loop to compute weigted gap at each node
   for i ∈ eachindex(gₙ)
      gₙ[i] = g[i, :] ⋅ normals[slave_dofs[i]]
   end

   # Initialize r_c
   #τ_c = zeros(eltype(X_float), size(X)) # sparse...?
   #τ_c = zeros(eltype(X_float), length(contact_dofs))
   g₀ = 0.0

   # ---------- #
   # ∫ᵧ 𝛅g λ dγ  #
   # ---------- #

   # Loop over master side dofs
   #for C in master_dofs
   for (i, A) in enumerate(slave_dofs)
    λ_A = penalty(g[i, :] ⋅ normals[slave_dofs[i]], ε)
      g₀ += (λ_A  * (1 / κ[i]))^p
   end

   # ---------------------------------- #
   # ∫ᵧ g 𝛅λ dγ = 0 for penalty methods  #
   # ---------------------------------- #

   return (g₀)^(1/p)
end

function contact_pnorm_ordered(X::AbstractVector{T1}, a::AbstractVector{T2}, ε,p) where {T1,T2}

    # Order  X
    X_ordered = getX_from_Dof_To_Node_order(dh, X)

    r_c = contact_pnorm(X_ordered, a, ε,p)

    return r_c
end
