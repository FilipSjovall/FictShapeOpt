const MortarSegmentation = Dict{Int, Vector{Tuple{Int, Vector{Real}}}}

function compute_segments(slave_element_ids, master_element_ids, elements,
                            element_types, coords, normals)
    S = MortarSegmentation()
    for sid in slave_element_ids
        haskey(S, sid) || (S[sid] = [])
        scon = elements[sid]
        xs1 = coords[scon[1]]
        xs2 = coords[scon[2]]
        ns1 = normals[scon[1]]
        ns2 = normals[scon[2]]
        for mid in master_element_ids
            mcon = elements[mid]
            xm1 = coords[mcon[1]]
            xm2 = coords[mcon[2]]
            # first project from slave to master, to find out
            # are we completely outside of domain
            # test
            norm( (xs1 + xs2)/2  - (xm1+xm2)/2 ) > 0.1 && continue
            # end of test
            ξ2a = project_from_slave_to_master(Val{:Seg2}, xs1, ns1, xm1, xm2)
            ξ2b = project_from_slave_to_master(Val{:Seg2}, xs2, ns2, xm1, xm2)
            ξ2a >  1.0 && ξ2b >  1.0 && continue
            ξ2a < -1.0 && ξ2b < -1.0 && continue
            ξ1a = project_from_master_to_slave(Val{:Seg2}, xm1, xs1, xs2, ns1, ns2)
            ξ1b = project_from_master_to_slave(Val{:Seg2}, xm2, xs1, xs2, ns1, ns2)
			ξ1a ≈ -1000 && ξ1b ≈ -1000 && continue # kommer aldrig hända?
            ξ1 = clamp.([ξ1a; ξ1b], -1.0, 1.0)
            isapprox(abs(ξ1[2]-ξ1[1]), 0.0) && continue
            push!(S[sid], (mid, ξ1))
        end
    end
    return S
end
