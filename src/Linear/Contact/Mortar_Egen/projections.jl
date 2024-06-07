function slave_to_master(xs, ns, xm1, xm2)
    # ------
    nom = ns[1]*(xm1[2] + xm2[2] - 2*xs[2]) - ns[2]*(xm1[1] + xm2[1] - 2*xs[1])
    denom = ns[1]*(xm1[2] - xm2[2]) + ns[2]*(xm2[1] - xm1[1])
    return nom/denom
end

function master_to_slave(xm, xs1, xs2, ns1, ns2)
    #---------------------------------------------------------
    # Borde bytas till analytisk lösn
    #---------------------------------------------------------
    N1s(ξ)  = (1.0 - ξ) / 2
    N2s(ξ)  = (1.0 + ξ) / 2
    Ns(ξ)   = [N1s(ξ) N2s(ξ)]

    dN1s(ξ) = -1.0 / 2.0
    dN2s(ξ) =  1.0 / 2.0
    dNs(ξ)  = [dN1s(ξ) dN2s(ξ)]

    xs(ξ)   = Ns(ξ)  * [vcat(xs1,0.0) vcat(xs2,0.0)]'
    dxs(ξ)  = dNs(ξ) * [vcat(xs1,0.0) vcat(xs2,0.0)]'

    ns(ξ)   = Ns(ξ)  * [ vcat(ns1, 0.0) vcat(ns2, 0.0) ]'
    dns(ξ)  = dNs(ξ) * [ vcat(ns1, 0.0) vcat(ns2, 0.0) ]'

    res(ξ)  = cross(vec(xs(ξ))-vec(vcat(xm,0.0)),vec(ns(ξ)))[3]
    dres(ξ) = cross(vec(dxs(ξ)), vec(ns(ξ)))[3] + cross(vec(xs(ξ)) - vec(vcat(xm, 0.0)), vec(dns(ξ)))[3]


    ξᵢ    = 0.0
    ξᵢ₊₁ = 0.0
    dξ    = 0.0
    for i ∈ 1:20
        dξ    = -res(ξᵢ) ./ dres(ξᵢ)
        ξᵢ₊₁  = ξᵢ + dξ
        if norm(ξᵢ₊₁ - ξᵢ) < 1e-6
            return ξᵢ₊₁
        end
        ξᵢ = ξᵢ₊₁
    end
    return -1000.
end
