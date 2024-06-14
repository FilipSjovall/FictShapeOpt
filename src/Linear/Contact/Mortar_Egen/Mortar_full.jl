

struct ContactParameters{e,coord,sel,sids,mids}
    elements :: e
    coords   :: coord
    slave_elements :: sel
    slave_element_ids :: sids
    master_element_ids :: mids
end


#struct ContactParameters
#    elements::Dict{Int, Vector{Int}}
#    coords::Vector{Coord}
#    slave_elements::Dict{Int, Vector{Int}}
#    slave_element_ids::Vector{Int}
#    master_element_ids::Vector{Int}
#end

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


# ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! #
# ContactParameters bara konstanter ! ! ! ! ! ! ! ! # Borde kunna skippa en stor del av kontakten med den här förenklingen + att det borde gå snabbare
# ContactParameters.coord = X (referenskoordinater) # Vidare kanske en closest_nodes eller liknande kan konstrueras för varje slave_element id eller något..
# ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! # Parallelisering om vi kan färlägga segmenteringen på något sätt.. kanske inte värt?
                                                    # Borde i varje fall fortsätta med reduktion av frihetsgrader


function compute_residual_and_tangent(ContactParameters,a)
    normals = calculate_normals(elements, X)

    for sid in ContactParameters.slave_element_ids
        sdofs   = contact_dofs[sid]
        xs1,xs2 = ContactParameters.coord[sid] + a[sdofs]
        ns1,ns2 = normals[sid]
        for mid in ContactParameters.master_element_ids
            mdofs   = contact_dofs[mid]
            xm1,xm2 = ContactParameters.coord[mid] + a[mdofs]
            distance([xs1,xs2],[xm1,xm2]) > tol && continue # if elements far apart -> skip to next
            xi1a = master_to_slave(xm1, xs1, xs2, ns1, ns2) # ?
            xi1b = master_to_slave(xm2, xs1, xs2, ns1, ns2) # ?
            xi1 = clamp.([xi1a; xi1b], -1.0, 1.0)
            l = 1/2*abs(xi1[2]-xi1[1])
            J = 1/2*norm(xs1-xs2) # ∂ξ / ∂η
            isapprox(l, 0.0) && continue # no contribution in this master element
            for ip in integration_points
            # for (w,η) in integration_points
                w   = weights[ip]
                η   = points[ip]
                ξ_s = (1-η)/2*ξ₁ + (1+η)/2*ξ₂
                j   = ... # ∂x¹e / ∂ξ
                N1  = [(1-ξ_s)/2 (1+ξ_s)/2]
                x_s = N1*[xs1 xs2]
                n_s = N1*[ns1 ns2]
                ξ_m = slave_to_master(xs, ns, xm1, xm2)
                N2  = [(1-ξ_m)/2 (1+ξ_m)/2] # ??
                x_m = N2*[xm1 xm2]
                # τ_s = [-n[2] n[1]]
                g   = n_s ⋅ (x_s - x_m) # saknas en formfunktion?
                λ   = ε*max(g⋅n,0)
                rc[sdofs] += N1'*N1*w*J*λ*n
                rc[mdofs] -= N1'*N2*w*J*λ*n
            end
        end
    end
end


function distance(xs,xm)

end
