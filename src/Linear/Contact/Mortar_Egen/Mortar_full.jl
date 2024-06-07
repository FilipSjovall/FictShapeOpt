

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


function compute_residual_and_tangent(ContactParameters,mfl)
    normals = calculate_normals(elements, X)

    for sid in ContactParameters.slave_element_ids
        sdofs   = contact_dofs[sid]
        xs1,xs2 = ContactParameters.coord[sid]
        ns1,ns2 = normals[sid]
        for mid in ContactParameters.master_element_ids
            mdofs   = contact_dofs[mid]
            xm1,xm2 = ContactParameters.coord[mid]
            xi1a = master_toslave(...)
            xi1b = master_toslave(...)
            xi1 = clamp.([xi1a; xi1b], -1.0, 1.0)
            l = 1/2*abs(xi1[2]-xi1[1])
            isapprox(l, 0.0) && continue # no contribution in this master element
            for ip in integration_points
                w   = weights[ip]
                ξ_s = (1-η)/2*ξ₁ + (1+η)/2*ξ₂
                N1  = [(1-ξ_s)/2 (1+ξ_s)/2]
                N2  = [(1-ξ_m)/2 (1+ξ_m)/2] # ??
                x_s = N1*[xs1 xs2]
                x_m = N2*[xm1 xm2]
                n_s = N1*[ns1 ns2]
                g   = n_s ⋅ (x_s - x_m)
                λ   = ε*max(g⋅n,0)
                rc[slave_dofs]  += N1'*N1*w*J*λ
                rc[master_dofs] -= N1'*N2*w*J*λ
            end
        end
    end
end
