function create_contact_list(dh,Γs,Γm, coord_dual)
    i                  = 0
    elements           = Dict{Int64, Vector{Int64}}()
    coords             = Dict{Int64, Vector{Real}}()
    slave_elements     = Dict{Int64, Vector{Int64}}()
    element_types      = Dict{Int64, Symbol}()
    slave_element_ids  = Vector{Int64}()
    master_element_ids = Vector{Int64}()
    for face in  Γs
       i        += 1
       face_el   = face[1]
       face_nods = Ferrite.faces(dh.grid.cells[face_el])[1]
       push!(elements,face_el=>[face_nods[1],face_nods[2]])
       push!(coords,face_nods[1]=>coord_dual[face_nods[1],:])
       push!(coords,face_nods[2]=>coord_dual[face_nods[2],:])
       push!(slave_element_ids,face_el)
       push!(slave_elements,face_el=>[face_nods[1],face_nods[2]])
       push!(element_types, face_el => :Seg2)
    end

    for face in  Γm
       i        += 1
       face_el   = face[1]
       face_nods = Ferrite.faces(dh.grid.cells[face_el])[1]
       push!(elements,face_el=>[face_nods[1],face_nods[2]])
       push!(coords,face_nods[1]=>coord_dual[face_nods[1],:])
       push!(coords,face_nods[2]=>coord_dual[face_nods[2],:])
       push!(master_element_ids,face_el)
       push!(element_types, face_el => :Seg2)
    end
    return elements,element_types, slave_elements, slave_element_ids, master_element_ids, coords
end

function penalty(g, ε)
    if g < 0.0
        p = 0.0
    else
        p = -ε * g # funkar utan minus
    end
   return p
end

function gap_function(X::AbstractVector{T}) where {T}
   # convert X to Real for compatibility with ForwardDiff
   X_float = real.(X)

   # Extract the coordinate vector (nbr_nodes x 2 )
   #coord  = getCoordVector(X_float)
   coordu = getCoordfromX(X_float)

   # Create dictionaries that are needed for the Mortar2D package
   elements, element_types, slave_elements, slave_element_ids, master_element_ids, coords = create_contact_list(dh, Γs, Γm, coordu)


   # Assemble D and M matrices and the slave and master dofs corresponding to the mortar segmentation
   slave_dofs, master_dofs, D, M = Mortar2D.calculate_mortar_assembly(elements, element_types, coords, slave_element_ids, master_element_ids)

   # Loop over slave dofs to compute the nodal gap vector.
   g0 = zeros(eltype(X_float), length(slave_dofs), 2)

   # Loops are fast with the LLVM compiler
   for (j, A) in (enumerate(slave_dofs))
      slave = [0; 0]
      for B in slave_dofs
         slave += D[A, B] * coords[B]
      end
      master = [0; 0]
      #for C in master_dofs
      for C in intersect(master_dofs,1:size(M,2))
         master += M[A, C] * coords[C]
      end
      # To compute the projected gap vector we multiply g[j,:] with the normal at node j
      g0[j, :] = slave - master
   end
   return g0
end

function gap_scaling(X::AbstractVector{T}) where {T}
   # convert X to Real for compatibility with ForwardDiff
   X_float = real.(X)

   # Extract the coordinate vector (nbr_nodes x 2 )
   coord = getCoordfromX(X_float)

   # Create dictionaries that are needed for the Mortar2D package
   elements, element_types, slave_elements, slave_element_ids, master_element_ids, coords = create_contact_list(dh, Γs, Γm, coord)

   # Assemble D and M matrices and the slave and master dofs corresponding to the mortar segmentation
   slave_nods, master_dofs, D, M = Mortar2D.calculate_mortar_assembly(elements, element_types, coords, slave_element_ids, master_element_ids)

   #  # Define scaling
   κ = ones(eltype(X_float), length(slave_nods))

   for (i, a) in enumerate(slave_nods)
      for (j, d) in enumerate(slave_nods)
         κ[i] += D[a, d]
      end
   end
   return κ
end

function contact_residual(X::AbstractVector{T1}, a::AbstractVector{T2}, ε::Number) where {T1,T2}

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

   # Assemble D and M matrices and the slave and master dofs corresponding to the mortar segmentation
   slave_dofs, master_dofs, D, M = Mortar2D.calculate_mortar_assembly(elements, element_types, coords, slave_element_ids, master_element_ids)



   # Compute the projected gap function
   g = zeros(eltype(X_float), length(slave_dofs), 2)

   # Loops are fast with the LLVM compiler
   for (j, A) in (enumerate(slave_dofs))
      slave = [0; 0]
      for B in intersect(slave_dofs, 1:size(D, 2))
         slave += D[A, B] * coords[B]
      end
      master = [0; 0]
      #for C in master_dofs
      for C in intersect(master_dofs, 1:size(M, 2))
         master += M[A, C] * coords[C]
      end
      # To compute the projected gap vector we multiply g[j,:] with the normal at node j
      g[j, :] = slave - master
   end


   # Initialize r_c
   r_c = zeros(eltype(X_float), size(X)) # sparse...?

   # ---------- #
   # ∫ᵧ 𝛅g λ dγ  #
   # ---------- #
   for (i, A) in enumerate(slave_dofs)
      λ_A = penalty(g[i, :] ⋅ normals[slave_dofs[i]], ε)
      for B in intersect(slave_dofs, 1:size(D, 2))
         B_dofs = register[B, :]  # Extract nodal degrees of freedom
         r_c[B_dofs] += D[A, B] * λ_A * normals[A] * (1 / κ[i])
      end
      for C in intersect(master_dofs, 1:size(M, 2))
         C_dofs = register[C, :] # Extract nodal degrees of freedom
         r_c[C_dofs] += -M[A, C] * λ_A * normals[A] * (1 / κ[i])
      end
   end

   # ---------------------------------- #
   # ∫ᵧ g 𝛅λ dγ = 0 for penalty methods  #
   # ---------------------------------- #
   return r_c
end

function contact_residual_simple(rc,a::AbstractVector{T}) where T
   rc = contact_residual(X_ordered, a, ε)
   nothing
end

function contact_residual_ordered(X::AbstractVector{T1}, a::AbstractVector{T2}, ε) where {T1,T2}

   # Order  X
   X_ordered = getX_from_Dof_To_Node_order(dh, X)

   r_c       = contact_residual(X_ordered, a, ε)

   return r_c
end

function contact_traction(X::AbstractVector{T1}, a::AbstractVector{T2}, ε) where {T1,T2}

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
   τ_c = zeros(eltype(X_float), size(coordu,1)) # sparse...?

   # ---------- #
   # ∫ᵧ 𝛅g λ dγ  #
   # ---------- #

   # Loop over master side dofs
   #for C in master_dofs
   for (i, A) in enumerate(slave_dofs)
      λ_A = penalty(g[i, :] ⋅ normals[slave_dofs[i]], ε)
      τ_c[A] = λ_A  * (1 / κ[i]) #  ∫ N^s N^s λ n dγ
        println("Traction | ", A, " ", λ_A, " normals | ", normals[slave_dofs[i]], " gap | ", gₙ[i])
   end
   @show coords
   # ---------------------------------- #
   # ∫ᵧ g 𝛅λ dγ = 0 for penalty methods  #
   # ---------------------------------- #

   return τ_c
end

function contact_residual_reduced(X::AbstractVector{T1}, a_c::AbstractVector{T2}, a_f::AbstractVector{T3}, ε::Number) where {T1,T2,T3}

    a_total  = similar(X)

    a_total[contact_dofs]  = a_c

    a_total[freec_dofs]    = a_f

    # Order displacements according to nodes and not dofs

    a_ordered = getDisplacementsOrdered(dh, a_total)

    # Scaling
    κ = gap_scaling(X)

    # convert X to Real for compatibility with ForwardDiff
    #X_float = real.(X)  + real.(a_ordered) # a ska vara sorterad på samma sätt som X, detta måste fixas!!!!!!!!!
    X_float = real.(X + a_ordered)

    # Extract the coordinate vector (nbr_nodes x 2 )
    coordu = getCoordfromX(X_float)

    # Create dictionaries that are needed for the Mortar2D package
    elements, element_types, slave_elements, slave_element_ids, master_element_ids, coords = create_contact_list(dh, Γs, Γm, coordu)

    # Compute nodal normals
    normals = Mortar2D.calculate_normals(elements, element_types, coords)

    # Assemble D and M matrices and the slave and master dofs corresponding to the mortar segmentation
    slave_dofs, master_dofs, D, M = Mortar2D.calculate_mortar_assembly(elements, element_types, coords, slave_element_ids, master_element_ids)

    # Compute the projected gap function
    g = zeros(eltype(X_float), length(slave_dofs), 2)

    # Loops are fast with the LLVM compiler
    for (j, A) in (enumerate(slave_dofs))
        slave = [0; 0]
        for B in slave_dofs
            slave += D[A, B] * coords[B]
        end
        master = [0; 0]
        #for C in master_dofs
        for C in intersect(master_dofs, 1:size(M, 2))
            master += M[A, C] * coords[C]
        end
        # To compute the projected gap vector we multiply g[j,:] with the normal at node j
        g[j, :] = slave - master
    end

    # Initialize r_c
    r_c = zeros(eltype(X_float), length(contact_dofs) ) # sparse...?

    # ---------- #
    # ∫ᵧ 𝛅g λ dγ  #
    # ---------- #
    for (i, A) in enumerate(slave_dofs)
        λ_A = penalty(g[i, :] ⋅ normals[slave_dofs[i]], ε)
        for (j,B) in enumerate(slave_dofs)
            # Extract nodal degrees of freedom

            #B_dofs = [2j-1, 2j]
            #B_dofs = [2order_nod[j]-1, 2order_nod[j]]
            nod = order[B]
            B_dofs = [2nod - 1, 2nod]
            r_c[B_dofs] += D[A, B] * λ_A * normals[A] * (1 / κ[i])
        end
        for (j,C) in enumerate(intersect(master_dofs, 1:size(M, 2)))
            # Extract nodal degrees of freedom

            #C_dofs = [2length(slave_dofs), 2length(slave_dofs)] + [2j-1, 2j]
            #C_dofs = [length(slave_dofs), length(slave_dofs)] + [2order_nod[j] - 1, 2order_nod[j]]
            nod = order[C]
            C_dofs = [2nod - 1, 2nod]
            r_c[C_dofs] += -M[A, C] * λ_A * normals[A] * (1 / κ[i])
        end
    end

    # ---------------------------------- #
    # ∫ᵧ g 𝛅λ dγ = 0 for penalty methods  #
    # ---------------------------------- #
    return r_c
end
