function create_contact_list(dh, Γs, Γm, coord_dual)
    i = 0
    elements           = Dict{Int64,Vector{Int64}}()
    coords             = Dict{Int64,Vector{Real}}()
    slave_elements     = Dict{Int64,Vector{Int64}}()
    element_types      = Dict{Int64,Symbol}()
    slave_element_ids  = Vector{Int64}()
    master_element_ids = Vector{Int64}()
    # Slave surface
    for face in Γs
        i += 1
        face_el = face[1]
        #face_nods = reverse(Ferrite.faces(dh.grid.cells[face_el])[face[2]]) # " face = (El, nod) "
        face_nods = Ferrite.faces(dh.grid.cells[face_el])[face[2]] # " face = (El, nod) "
        push!(elements,face_el => [face_nods[1],face_nods[2]])
        push!(coords, face_nods[1] => coord_dual[face_nods[1], :])
        push!(coords, face_nods[2] => coord_dual[face_nods[2], :])
        push!(slave_element_ids, face_el)
        push!(slave_elements, face_el => [face_nods[1], face_nods[2]])
        push!(element_types, face_el => :Seg2)
    end
    # Master surface
    for face in Γm
        i += 1
        face_el = face[1]
        face_nods = Ferrite.faces(dh.grid.cells[face_el])[face[2]]
        push!(elements, face_el => [face_nods[1], face_nods[2]])
        push!(coords, face_nods[1] => coord_dual[face_nods[1], :])
        push!(coords, face_nods[2] => coord_dual[face_nods[2], :])
        push!(master_element_ids, face_el)
        push!(element_types, face_el => :Seg2)
    end
    return elements, element_types, slave_elements, slave_element_ids, master_element_ids, coords
end

function penalty(g, ε)
    if g < 0.0
        p = 0.0 #0.1* ε * g
    else
        p = ε * g #
    end
    return p
end

# softplus(x) = log(exp(x) + 1)
# gelu(x) = 0.5x * (1 + tanh(√(2 / π) * (x + 0.044715x^3)))
# mish(x) = x * tanh(softplus(x))

function penalty_filter(g, ε)
    if g < 0.0
        # p = 0.
        if g ≥ -0.0001
           #p = 0.1 * ε * g
           p = 0.0
        else
           p = 0.0
        end
    else
        p = ε * g #
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
    #for (j, A) in (enumerate(slave_dofs))
    # for (j, A) in (enumerate(intersect(slave_dofs, 1:min(size(D, 2), size(M, 1)))))
    # #for (j, A) in (enumerate(intersect(slave_dofs, 1:min(size(D, 2)))))
    #     slave = [0.; 0.]
    #     for B in slave_dofs
    #         slave += D[A, B] * coords[B]
    #     end
    #     master = [0.; 0.]
    #     #for C in master_dofs
    #     for C in intersect(master_dofs, 1:size(M, 2))
    #         master += M[A, C] * coords[C]
    #     end
    #     # To compute the projected gap vector we multiply g[j,:] with the normal at node j
    #     g0[j, :] = slave - master
    # end
    for (j, A) in (enumerate(intersect(slave_dofs, 1:min(size(D, 1)))))
        slave = [0; 0]
        for (jj,B) in (enumerate(intersect(slave_dofs, 1:min(size(D, 1))))) # slave_dofs
            slave += D[A, B] * coords[B]
        end
        master = [0; 0]
        #for C in master_dofs
        for C in intersect(master_dofs, 1:size(M, 2))
            master += M[A, C] * coords[C]
        end
        # To compute the projected gap vector we multiply g[j,:] with the normal at node j
        g0[j, :] = slave - master
    end
    return g0
end

function gap_scaling(X::AbstractVector{T}) where {T}
    #
    # Vi måste räkna ut M,D igen här för det odeformerade tillstånded, möjligtvis kan detta sparas undan...
    #

    # convert X to Real for compatibility with ForwardDiff
    X_float = real.(X)

    # Extract the coordinate vector (nbr_nodes x 2 )
    coord = getCoordfromX(X_float)

    # Create dictionaries that are needed for the Mortar2D package
    elements, element_types, slave_elements, slave_element_ids, master_element_ids, coords = create_contact_list(dh, Γs, Γm, coord)

    # Assemble D and M matrices and the slave and master dofs corresponding to the mortar segmentation
    slave_nods, master_dofs, D, M = Mortar2D.calculate_mortar_assembly(elements, element_types, coords, slave_element_ids, master_element_ids)

    #  # Define scaling
    κ = zeros(eltype(X_float), length(slave_nods))
    for (i, a) in (enumerate(intersect(slave_nods, 1:min(size(D, 2)))))     #
        for (j, d) in (enumerate(intersect(slave_nods, 1:min(size(D, 2))))) #
            κ[i] += D[a, d]
        end
    end
    #@show κ
    #κ .= 1
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
    #for (j, A) in (enumerate(slave_dofs))
    for (j, A) in (enumerate(intersect(slave_dofs, 1:min(size(D, 2)))))
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
    #for (i, A) in enumerate(slave_dofs)
    for (i, A) in (enumerate(intersect(slave_dofs, 1:min(size(D, 2), size(M, 1)))))
        λ_A = penalty(g[i, :] ⋅ normals[slave_dofs[i]], ε)
        if (λ_A != 0 && κ[i] != 0)
            for B in intersect(slave_dofs, 1:size(D, 2))
                B_dofs = register[B, :]  # Extract nodal degrees of freedom
                r_c[B_dofs] += D[A, B] * λ_A * normals[A] * (1 / κ[i])
            end
            for C in intersect(master_dofs, 1:size(M, 2))
                C_dofs = register[C, :] # Extract nodal degrees of freedom
                r_c[C_dofs] += -M[A, C] * λ_A * normals[A] * (1 / κ[i])
            end
        end
    end

    # ---------------------------------- #
    # ∫ᵧ g 𝛅λ dγ = 0 for penalty methods #
    # ----------------------------------#
    return r_c
end

function contact_residual_simple(rc, a::AbstractVector{T}) where {T}
    rc = contact_residual(X_ordered, a, ε)
    nothing
end

function contact_residual_ordered(X::AbstractVector{T1}, a::AbstractVector{T2}, ε) where {T1,T2}

    # Order  X
    X_ordered = getX_from_Dof_To_Node_order(dh, X)

    r_c = contact_residual(X_ordered, a, ε)

    return r_c
end

function contact_traction(X::AbstractVector{T1}, a::AbstractVector{T2}, ε) where {T1,T2}

    # Order displacements according to nodes and not dofs
    a_ordered = getDisplacementsOrdered(dh, a)

    # Scaling
    κ = gap_scaling(X)
    #κ .= 0.005
    # convert X to Real for compatibility with ForwardDiff
    X_float = real.(X + a_ordered) # a ska vara sorterad på samma sätt som X, stämmer väl nu?

    # Extract the coordinate vector (nbr_nodes x 2 )
    coordu = getCoordfromX(X_float)

    # Create dictionaries that are needed for the Mortar2D package
    elements, element_types, slave_elements, slave_element_ids, master_element_ids, coords = create_contact_list(dh, Γs, Γm, coordu)

    # Compute nodal normals
    normals = Mortar2D.calculate_normals(elements, element_types, coords)

    # Compute the projected gap function
    g = gap_function(X_float)

    # Assemble D and M matrices and the slave and master dofs corresponding to the mortar segmentation
    slave_dofs, master_dofs, D, M = Mortar2D.calculate_mortar_assembly(elements, element_types, coords, slave_element_ids, master_element_ids)

    # Initialize the nodal gap vector.
    gₙ = zeros(eltype(X_float), length(slave_dofs))

    # Loop to compute weigted gap at each node
    for i ∈ eachindex(gₙ)
        gₙ[i] = g[i, :] ⋅ normals[slave_dofs[i]]
    end

    # Initialize r_c
    τ_c = Dict{Int64,Real}()

    # ---------- #
    # ∫ᵧ 𝛅g λ dγ #
    # ---------- #

    # Loop over master side dofs
    for (i, A) in enumerate(slave_dofs)
        λ_A = penalty(g[i, :] ⋅ normals[A] , ε)
        #@show λ_A 1/κ[i] g[i, :] ⋅ normals[A]
        if λ_A != 0 && κ[i] != 0
            τ = λ_A * normals[A]/ κ[i]
            #push!(τ_c, A => λ_A )
            push!(τ_c, A => λ_A / κ[i])
            #push!(τ_c, A => τ[2])
        else
            push!(τ_c, A => 0.0)
        end
        #
        #push!(τ_c, A => g[i,:] ⋅ [0.0 1.0] / κ[i])
    end
    # ---------------------------------- #
    # ∫ᵧ g 𝛅λ dγ = 0 for penalty methods #
    # ---------------------------------- #
    return τ_c
end

function contact_traction_vector(X::AbstractVector{T1}, a::AbstractVector{T2}, ε) where {T1,T2}

    # Order displacements according to nodes and not dofs
    a_ordered = getDisplacementsOrdered(dh, a)

    # Scaling
    κ = gap_scaling(X)
    #κ .= 0.005
    # convert X to Real for compatibility with ForwardDiff
    X_float = real.(X + a_ordered) # a ska vara sorterad på samma sätt som X, stämmer väl nu?

    # Extract the coordinate vector (nbr_nodes x 2 )
    coordu = getCoordfromX(X_float)

    # Create dictionaries that are needed for the Mortar2D package
    elements, element_types, slave_elements, slave_element_ids, master_element_ids, coords = create_contact_list(dh, Γs, Γm, coordu)

    # Compute nodal normals
    normals = Mortar2D.calculate_normals(elements, element_types, coords)

    # Compute the projected gap function
    g = gap_function(X_float)

    # Assemble D and M matrices and the slave and master dofs corresponding to the mortar segmentation
    slave_dofs, master_dofs, D, M = Mortar2D.calculate_mortar_assembly(elements, element_types, coords, slave_element_ids, master_element_ids)

    # Initialize the nodal gap vector.
    gₙ = zeros(eltype(X_float), length(slave_dofs))

    # Loop to compute weigted gap at each node
    for i ∈ eachindex(gₙ)
        gₙ[i] = g[i, :] ⋅ normals[slave_dofs[i]]
    end

    # Initialize r_c
    #τ_c = Dict{Int64,Real}()
    τ_c = Dict{Int64, Vector{Float64} }()

    # ---------- #
    # ∫ᵧ 𝛅g λ dγ #
    # ---------- #
    # Loop over master side dofs
    for (i, A) in enumerate(slave_dofs)
        λ_A = penalty(g[i, :] ⋅ normals[A] , ε)
        if (λ_A != 0 && κ[i] != 0)
            τ = λ_A * normals[A]/ κ[i]
            push!(τ_c, A => τ)
        else
            τ = 0.0 * normals[A]
            push!(τ_c, A => τ)
        end
    end
    # ---------------------------------- #
    # ∫ᵧ g 𝛅λ dγ = 0 for penalty methods #
    # ---------------------------------- #
    return τ_c
end


function contact_residual_reduced(X::AbstractVector{T1}, a_c::AbstractVector{T2}, a_f::AbstractVector{T3}, ε::Number) where {T1,T2,T3}
    a_total = similar(X)
    a_total.= 0

    # Contact dofs for AD
    a_total[contact_dofs] = a_c

    # Non-contact dofs
    a_total[freec_dofs] = a_f

    # Order displacements according to nodes and not dofs
    a_ordered = getDisplacementsOrdered(dh, a_total)

    # Scaling
    κ = gap_scaling(X)

    # convert X to Real for compatibility with ForwardDiff
    X_float = real.(X + a_ordered)

    # Extract the coordinate vector (nbr_nodes x 2 )
    coordu = getCoordfromX(X_float)

    # Create dictionaries that are needed for the Mortar2D package
    # byt mot update_coords ... ? kanske kan loopa över nycklarna i Contact.coords med...
    elements, element_types, slave_elements, slave_element_ids, master_element_ids, coords = create_contact_list(dh, Γs, Γm, coordu)

    # Compute nodal normals
    normals = Mortar2D.calculate_normals(elements, element_types, coords)

    # Assemble D and M matrices and the slave and master dofs corresponding to the mortar segmentation
    slave_dofs, master_dofs, D, M = Mortar2D.calculate_mortar_assembly(elements, element_types, coords, slave_element_ids, master_element_ids)

    # Compute the projected gap function
    g = zeros(eltype(X_float), length(slave_dofs), 2)
    #
    # byt mot funktion?
    for (j, A) in (enumerate(intersect(slave_dofs, 1:min(size(D, 1)))))
        slave = [0; 0]
        for (jj,B) in (enumerate(intersect(slave_dofs, 1:min(size(D, 1))))) # slave_dofs
            slave += D[A, B] * coords[B]
        end
        master = [0; 0]
        for C in intersect(master_dofs, 1:size(M, 2))
            master += M[A, C] * coords[C]
        end
        # To compute the projected gap vector we multiply g[j,:] with the normal at node j
        g[j, :] = slave - master
    end
    r_c = zeros(eltype(X_float), length(contact_dofs)) # sparse...?
    # ----------- #
    # ∫ᵧ λ 𝛅g dγ  #
    # ----------- #
    # skall detta vara en funktion?
    for (i, A) in (enumerate(intersect(slave_dofs, 1:min(size(D, 1)))))
        λ_A = penalty(g[i, :] ⋅ normals[slave_dofs[i]], ε)
        if ( λ_A != 0.0 && κ[i] != 0 )
            for (j, B) in (enumerate(intersect(slave_dofs, 1:size(D, 2))))
                # Extract nodal degrees of freedom
                nod          = order[B]
                B_dofs       = [2nod - 1, 2nod]
                r_c[B_dofs] += D[A, B] * λ_A * normals[A] * ( 1 / κ[i] )
            end
            for (jj, C) in enumerate(intersect(master_dofs, 1:size(M, 2)))
                # Extract nodal degrees of freedom
                nod          = order[C]
                C_dofs       = [2nod - 1, 2nod]
                r_c[C_dofs] += -M[A, C] * λ_A * normals[A] * ( 1 / κ[i] )
            end
        end
    end
    return r_c
end


function contact_residual_reduced_filter(X::AbstractVector{T1}, a_c::AbstractVector{T2}, a_f::AbstractVector{T3}, ε::Number) where {T1,T2,T3}

    a_total = similar(X)

    a_total[contact_dofs] = a_c

    a_total[freec_dofs] = a_f

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
    #for (j, A) in (enumerate(slave_dofs))
    for (j, A) in (enumerate(intersect(slave_dofs, 1:min(size(D, 2)))))
        slave = [0; 0]
        for B in intersect(slave_dofs, 1:min(size(D, 2)))
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
    #@show g
    #@show normals
    #@show D
    #@show M
    # Initialize r_c
    r_c = zeros(eltype(X_float), length(contact_dofs)) # sparse...?

    # ---------- #
    # ∫ᵧ 𝛅g λ dγ  #
    # ---------- #
    #for (i, A) in enumerate(slave_dofs)
    for (i, A) in (enumerate(intersect(slave_dofs, 1:min(size(D, 2)))))
        λ_A = penalty_filter(g[i, :] ⋅ normals[slave_dofs[i]], ε)
        if ( λ_A != 0 && κ[i] != 0 )
            for (j, B) in (enumerate(intersect(slave_dofs, 1:size(D, 2))))
                # Extract nodal degrees of freedom
                nod = order[B]
                B_dofs = [2nod - 1, 2nod]
                r_c[B_dofs] += D[A, B] * λ_A * normals[A] * (1 / κ[i])
            end
            for (j, C) in enumerate(intersect(master_dofs, 1:size(M, 2)))
                # Extract nodal degrees of freedom
                nod = order[C]
                C_dofs = [2nod - 1, 2nod]
                r_c[C_dofs] += -M[A, C] * λ_A * normals[A] * (1 / κ[i])
            end
        end
    end
    # ---------------------------------- #
    # ∫ᵧ g 𝛅λ dγ = 0 for penalty methods  #
    # ---------------------------------- #
    return r_c
end
