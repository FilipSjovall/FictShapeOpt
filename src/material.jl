# Material routines
#
# Current implementation includes only 2D plane strain for Neo-Hookean strain energy
#
using LinearAlgebra

function neohooke1!(S,eff,mp)

    Kₘ = mp[1]
    Gₘ = mp[2]

    D₁ = 2/Kₘ

    Fₚ  = zeros(3,3)

    Fₚ[1,:] = [eff[1] eff[2] 0.0]
    Fₚ[2,:] = [eff[3] eff[4] 0.0]
    Fₚ[3,:] = [0.0 0.0 1.0]

    C = transpose(Fₚ)*Fₚ

    C⁻¹ = inv(C)
    J   = det(Fₚ)

    S   = Kₘ/2.0 * (J^2 - 1.0)*C⁻¹ + Gₘ*J^(-2/3)*( I - tr(C)/3.0*C⁻¹ )
    return S
end

function neohookeS!(Sₑ,eff,mp)
   ngp = size(eff,2)
   for gp ∈ 1 : ngp 
        neohooke1!(Sₑ[:,gp],eff[:,gp],mp)
   end
   return Sₑ
end

function dneohooke1(eff,mp)
    D  = zeros(3,3) 
    Kₘ = mp[1]
    Gₘ = mp[2]

    D₁ = 2/Kₘ

    F  = zeros(3,3)

    F[1,:] = [eff[1] eff[2] 0.0]
    F[2,:] = [eff[3] eff[4] 0.0]
    F[3,:] = [0.0 0.0 1.0]

    C = transpose(F)*F

    C⁻¹ = inv(C)
    J   = det(F)

    a₁ = Kₘ * J * 2.0 + 2.0/9.0 * Gₘ * J^(-2/3) * tr(C)
    a₂ = 2.0 * Gₘ / 3.0 * J^(-2.0/3.0)
    a₃ = Gₘ/3.0 * J^(-2/3)*tr(C) - Kₘ/2.0 * (J^2 - 1.0)

    ix = Matrix{Int32}(undef,6,4)
    ix[1,:] = [1 1 1 1]
    ix[2,:] = [1 1 2 2]
    ix[3,:] = [1 1 1 2]
    ix[4,:] = [2 2 2 2]
    ix[5,:] = [2 2 1 2]
    ix[6,:] = [1 2 1 2]

        Ds = zeros(6)
        for el ∈ 1 : 6
            i  = ix[el,1]
            jj = ix[el,2]
            kk = ix[el,3]
            l  = ix[el,4]
            Ds[el] = a₁ * C⁻¹[i,jj] * C⁻¹[kk,l] - 
                     a₂ * ( I[i,jj] * C⁻¹[kk,l] + C⁻¹[i,jj] * I[kk,l] ) +
                     a₃ * ( C⁻¹[i,kk] * C⁻¹[jj,l] + C⁻¹[i,l] + C⁻¹[jj,kk] ) 
        end

        

        D[1,:] = [Ds[1] Ds[2] Ds[3]]
        D[2,:] = [Ds[2] Ds[4] Ds[5]]
        D[3,:] = [Ds[3] Ds[5] Ds[6]]
        return D
end

function dneohookeD!(Dgp,eff,mp)
    ngp = size(eff,2)
    for gp ∈ 1 : ngp 
       Dgp[:,:,gp] = dneohooke1(eff[:,gp], mp)
    end
    return Dgp
end
