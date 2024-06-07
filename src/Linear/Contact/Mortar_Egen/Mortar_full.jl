

struct ContactParameters{e,coord,sel,sids,mids}
    elements :: e
    coords   :: coord
    slave_elements :: sel
    slave_element_ids :: sids
    master_element_ids :: mids
end


struct ContactParameters
    elements::Dict{Int, Vector{Int}}
    coords::Vector{Coord}
    slave_elements::Dict{Int, Vector{Int}}
    slave_element_ids::Vector{Int}
    master_element_ids::Vector{Int}
end

function create_contact_list(dh<:DofHandler, Γs, Γm, coord_dual)
    elements = Dict{Int, Vector{Int}}()
    coords = Dict{Int, Vector{Real}}()
    slave_elements = Dict{Int, Vector{Int}}()
    slave_element_ids = Int[]
    master_element_ids = Int[]

    for (i, face) in enumerate(Γs)
        face_el = face[1]
        face_nods = Ferrite.faces(dh.grid.cells[face_el])[face[2]]
        elements[face_el] = [face_nods[1], face_nods[2]]
        coords[face_nods[1]] = coord_dual[face_nods[1], :]
        coords[face_nods[2]] = coord_dual[face_nods[2], :]
        push!(slave_element_ids, face_el)
        slave_elements[face_el] = [face_nods[1], face_nods[2]]
    end

    for (i, face) in enumerate(Γm)
        face_el = face[1]
        face_nods = Ferrite.faces(dh.grid.cells[face_el])[face[2]]
        elements[face_el] = [face_nods[1], face_nods[2]]
        coords[face_nods[1]] = coord_dual[face_nods[1], :]
        coords[face_nods[2]] = coord_dual[face_nods[2], :]
        push!(master_element_ids, face_el)
    end

    return ContactParameters(elements, coords, slave_elements, slave_element_ids, master_element_ids)
end

function update_coords!(params::ContactParameters, new_coords::Vector{Coord})
    params.coords = new_coords
end


function compute_face_coords(dh<:DofHandler, Γs, Γm, coord_dual)
    for (i, face) in enumerate(Γs)
        face_el = face[1]
        face_nods = Ferrite.faces(dh.grid.cells[face_el])[face[2]]
        coords_new[face_nods[1]] = coord_dual[face_nods[1], :]
        coords_new[face_nods[2]] = coord_dual[face_nods[2], :]
    end
    for (i, face) in enumerate(Γm)
        face_el = face[1]
        face_nods = Ferrite.faces(dh.grid.cells[face_el])[face[2]]
        coords_new[face_nods[1]] = coord_dual[face_nods[1], :]
        coords_new[face_nods[2]] = coord_dual[face_nods[2], :]
    end
    return coords_new
end


Δn = ForwardDiff.gradient(x->calculate_normals(elements,x),X) # funkar?

Δg = ForwardDiff.gradient(x->gap_function(elements,x),X) # funkar?

ΔM =
ΔD =
