function mortar_matrix(slave_element_id, elements, element_types,
                                   coords, normals, segmentation)
    # Initialization + calculate jacobian
    De = zeros(2, 2)
    Me = Dict{Int, Matrix{Real}}()
    scon = elements[slave_element_id]
    xs1 = coords[scon[1]]
    xs2 = coords[scon[2]]
    ns1 = normals[scon[1]]
    ns2 = normals[scon[2]]
    J = norm(xs2-xs1)/2.0
    # 1. calculate De
    for (mid, (xi1, xi2)) in segmentation[slave_element_id]
        s = abs(xi2-xi1)/2.0
        for (w, ip) in seg_integration_points[3]
            xi_s = (1-ip)/2*xi1 + (1+ip)/2*xi2
            N1 = [(1-xi_s)/2 (1+xi_s)/2]
            De += w*N1'*N1*s*J
        end
    end
    # 2. calculate Me
    for (mid, (xi1, xi2)) in segmentation[slave_element_id]
        Me[mid] = zeros(2, 2)
        mcon = elements[mid]
        xm1 = coords[mcon[1]]
        xm2 = coords[mcon[2]]
        s = abs(xi2-xi1)/2.0
        for (w, xi) in seg_integration_points[3]
            xi_s = (1-xi)/2*xi1 + (1+xi)/2*xi2
            N1 = [(1-xi_s)/2 (1+xi_s)/2]
            n_s = N1[1]*ns1 + N1[2]*ns2
            x_g = N1[1]*xs1 + N1[2]*xs2
            xi_m = project_from_slave_to_master(Val{:Seg2}, x_g, n_s, xm1, xm2)
            N2 = [(1-xi_m)/2 (1+xi_m)/2]
            Me[mid] += w*N1'*N2*s*J
        end
    end
    return De, Me


end
