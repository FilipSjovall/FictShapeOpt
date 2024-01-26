# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/Mortar2D.jl/blob/master/LICENSE

"""
    project_from_slave_to_master(Val{:Seg2}, xs, ns, xm1, xm2)

Find the projection of a slave node `xs`, having normal vector `ns`, onto master
elements with nodes (`xm1`, `xm2`). Returns master element dimensionless parameter
xi, that is,

    xm = 1/2*(1-xi)*xm1 + 1/2*(1+xi)*xm2

# Example

```jldoctest
xm1 = [7.0, 2.0]
xm2 = [4.0, -2.0]
xs = [0.0, 0.0]
ns = [3.0/5.0, -4.0/5.0]
xi2 = project_from_slave_to_master(Val{:Seg2}, xs, ns, xm1, xm2)
round(xi2, digits=6)

# output

1.833333
```

"""
function project_from_slave_to_master(::Type{Val{:Seg2}}, xs, ns, xm1, xm2)
    # ------
    # Egen lösning på  om ytor är för långt ifrån varandra
    #if norm(xs1 - xm) > 0.2 || norm(xs2 - xm) > 0.2
    #    return -2
    #end
    # ------
    nom = ns[1]*(xm1[2] + xm2[2] - 2*xs[2]) - ns[2]*(xm1[1] + xm2[1] - 2*xs[1])
    denom = ns[1]*(xm1[2] - xm2[2]) + ns[2]*(xm2[1] - xm1[1])
    return nom/denom
end

"""
    project_from_master_to_slave(Val{:Seg2}, xm, xs1, xs2, ns1, ns2)

Find the projection of a master node `xm`, to the slave surface with nodes
(`xs1`, `xs2`), in direction of slave surface normal defined by (`ns1`, `ns2`).
Returns slave element dimensionless parameter, that is, to find coordinates
in slave side:

    xs = 1/2*(1-xi)*xs1 + 1/2*(1+xi)*xs2

# Example

```jldoctest
xm = [4.0, -2.0]
xs1 = [0.0, 0.0]
xs2 = [4.0, 3.0]
ns1 = [3.0/5.0, -4.0/5.0]
ns2 = sqrt(2)/2*[1, -1]
xi1 = project_from_master_to_slave(Val{:Seg2}, xm, xs1, xs2, ns1, ns2)
round(xi1, digits=6)

# output

-0.281575
```

"""
function project_from_master_to_slave(::Type{Val{:Seg2}}, xm, xs1, xs2, ns1, ns2)
    #---------------------------------------------------------
        #=
        # Special case when normal is constant. Then we have
        # unique solution and linear equation to solve.
        if isapprox(ns1, ns2,rtol=1e-6) # ska denna vara 10^-6 istället?
            nom = -2*ns1[1]*xm[2] + ns1[1]*xs1[2] + ns1[1]*xs2[2] + 2*ns1[2]*xm[1] - ns1[2]*xs1[1] - ns1[2]*xs2[1]
            denom = ns1[1]*xs1[2] - ns1[1]*xs2[2] - ns1[2]*xs1[1] + ns1[2]*xs2[1]
            return nom/denom
        end

        nom_b = 2*ns1[1]*xm[2] - 2*ns1[1]*xs1[2] - 2*ns1[2]*xm[1] + 2*ns1[2]*xs1[1] - 2*ns2[1]*xm[2] + 2*ns2[1]*xs2[2] + 2*ns2[2]*xm[1] - 2*ns2[2]*xs2[1]
        denom_b = ns1[1]*xs1[2] - ns1[1]*xs2[2] - ns1[2]*xs1[1] + ns1[2]*xs2[1] - ns2[1]*xs1[2] + ns2[1]*xs2[2] + ns2[2]*xs1[1] - ns2[2]*xs2[1]
        nom_c = -2*ns1[1]*xm[2] + ns1[1]*xs1[2] + ns1[1]*xs2[2] + 2*ns1[2]*xm[1] - ns1[2]*xs1[1] - ns1[2]*xs2[1] - 2*ns2[1]*xm[2] + ns2[1]*xs1[2] + ns2[1]*xs2[2] + 2*ns2[2]*xm[1] - ns2[2]*xs1[1] - ns2[2]*xs2[1]
        denom_c = ns1[1]*xs1[2] - ns1[1]*xs2[2] - ns1[2]*xs1[1] + ns1[2]*xs2[1] - ns2[1]*xs1[2] + ns2[1]*xs2[2] + ns2[2]*xs1[1] - ns2[2]*xs2[1]

        a = 1.0
        b = nom_b / denom_b
        c = nom_c / denom_c
        d = b^2 - 4*a*c

        if d < 0
            @warn("Mortar2D.calculate_projections(): negative discriminant $d")
            @warn("xm = $xm, xs1 = $xs1, xs2 = $xs2, ns1 = $ns1, ns2 = $ns2")
        end
        if d < 0
            error("ERROR: DomainError - negative discriminant")
        else
            sols = [-b + sqrt(d), -b - sqrt(d)]/(2.0*a)
        end
        return sols[argmin(abs.(sols))]
        =#
    #----------------------------
    ## #=


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
        dξ   = -res(ξᵢ) ./ dres(ξᵢ)
        dξ   = clamp.(dξ, -0.5, 0.5)
        ξᵢ₊₁ = clamp.(ξᵢ + dξ, -1.0, 1.0)
        if norm(ξᵢ₊₁ - ξᵢ) < 1e-6 #&& norm(res(ξᵢ₊₁)) < 1e-6 && norm(dξ)<1e-6
            #@show ξᵢ norm(ξᵢ₊₁ - ξᵢ)
            return ξᵢ₊₁
        end
        ξᵢ = ξᵢ₊₁
    end

    return -10.
    ## =#
end
