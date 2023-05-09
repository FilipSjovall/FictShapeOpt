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

