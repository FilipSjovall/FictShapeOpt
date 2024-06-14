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

    # Δx && Δy

    # Förtydliga n1 // n2

    a =   Δx[1] * ( n1[2] - n2[2] ) + Δx[2] * ( n2[2] - n1[2] )
        - Δy[1] * ( n1[1] - n2[2] ) - Δy[2] * ( n2[1] - n1[1] )
    b =   Δx[2] * n2[2] - Δx[2] * n2[2] + Δy[1] * n1[1] - Δy[2] * n2[1]
    c =   Δx[1] * ( n1[2] + n2[2] ) + Δx[2] * ( n2[2] + n1[2] )
        - Δy[1] * ( n1[1] + n2[2] ) - Δy[2] * ( n2[1] + n1[1] )
    if (b^2 -4 * a * c) ≈ 0
        return ξₚ = -b/2a
    end
    ξp₁ = ( - b + √(b^2 -4 * a * c) )/2a
    ξp₂ = ( - b - √(b^2 -4 * a * c) )/2a

end
