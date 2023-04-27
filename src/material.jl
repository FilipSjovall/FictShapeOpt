# Material routines
#
# Current implementation includes only 2D plane strain for Neo-Hookean strain energy
#
using LinearAlgebra

function neohooke1(eff,mp)

    Kₘ = mp[1]
    Gₘ = mp[2]

    #D₁ = 2/Kₘ
    C   = zeros(3,3)
    Fₚ  = zeros(3,3)

    Fₚ[1,:] = [eff[1] eff[2] 0.0]
    Fₚ[2,:] = [eff[3] eff[4] 0.0]
    Fₚ[3,:] = [0.0 0.0 1.0]

    C = transpose(Fₚ)*Fₚ

    C⁻¹ = inv(C)
    J   = det(Fₚ)

    Stress   = Kₘ/2.0 * (J^2 - 1.0).*C⁻¹ + Gₘ*(J^(-2/3)).*( I(3) - tr(C)/3.0.*C⁻¹ )
    S        = [Stress[1,1]; Stress[2,2]; Stress[3,3]; Stress[1,2]]
    return S
end

function neohookeS!(Sₑ,eff,mp)
   ngp = size(eff,2)
   for gp ∈ 1 : ngp
        Sₑ[:,gp] = neohooke1(eff[:,gp],mp)
   end
end

function dneohooke1(eff,mp)
    D  = zeros(3,3) 
    Kₘ = mp[1]
    Gₘ = mp[2]

    #D₁ = 2/Kₘ

    F  = zeros(3,3)

    F[1,:] = [eff[1] eff[2] 0.0]
    F[2,:] = [eff[3] eff[4] 0.0]
    F[3,:] = [0.0 0.0 1.0]

 
    C = transpose(F)*F

    #C_inv = inv(C)
    C⁻¹ = inv(C)
    J   = det(F)

    a₁ = Kₘ * J ^ 2.0 + 2.0/9.0 * Gₘ * J^(-2/3) * tr(C)
    a₂ = 2.0 * Gₘ / 3.0 * J^(-2.0/3.0)
    a₃ = Gₘ/3.0 * J^(-2/3)*tr(C) - Kₘ/2.0 * (J^2 - 1.0)

    #D1111 = a1*C_inv[1,1]*C_inv[1,1] - a2*(C_inv[1,1]+C_inv[1,1]) + a3*(C_inv[1,1]*C_inv[1,1]+C_inv[1,1]*C_inv[1,1]);
    #D1122 = a1*C_inv[1,1]*C_inv[2,2] - a2*(C_inv[2,2]+C_inv[1,1]) + a3*(C_inv[1,2]*C_inv[1,2]+C_inv[1,2]*C_inv[1,2]);
    #D1112 = a1*C_inv[1,1]*C_inv[1,2] - a2* C_inv[1,2]             + a3*(C_inv[1,1]*C_inv[1,2]+C_inv[1,2]*C_inv[1,1]);
    #D2222 = a1*C_inv[2,2]*C_inv[2,2] - a2*(C_inv[2,2]+C_inv[2,2]) + a3*(C_inv[2,2]*C_inv[2,2]+C_inv[2,2]*C_inv[2,2]);
    #D1222 = a1*C_inv[1,2]*C_inv[2,2] - a2* C_inv[1,2]             + a3*(C_inv[1,2]*C_inv[2,2]+C_inv[1,2]*C_inv[2,2]);
    #D1212 = a1*C_inv[1,2]*C_inv[1,2]                              + a3*(C_inv[1,1]*C_inv[2,2]+C_inv[1,2]*C_inv[2,1]);
    #D     = [D1111 D1122 D1112 ;
    #         D1122 D2222 D1222 ;
    #         D1112 D1222 D1212];
    
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
                 a₃ * ( C⁻¹[i,kk] * C⁻¹[jj,l] + C⁻¹[i,l] * C⁻¹[jj,kk] ) 
    end
#
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
end


function StVenantS(eff,mp)
    D  = zeros(3,3) 
    Kₘ = mp[1]
    Gₘ = mp[2]

    #D₁ = 2/Kₘ
    E = 1e0
    ν = 0.3

    my= E/(2*(1+ν)); 
    k = E/ (3*(1-2*ν));

    F  = zeros(3,3)

    F[1,:] = [eff[1] eff[2] 0.0]
    F[2,:] = [eff[3] eff[4] 0.0]
    F[3,:] = [0.0 0.0 1.0]

    display(F)
 
    C = transpose(F)*F

    Eg = 1/2 * ( C - I(3) )
    I1 = tr(Eg)
    #S = Kₘ*tr(E)*I(3) + 2*Gₘ*E
    S  = 2*my*(Eg + ν/(1-2*ν) * I1 * I(3));

    return [S[1,1] S[2,2] S[3,3] S[1,2]]
end

function dStVenant(eff,mp)
    D  = zeros(3,3) 
    K = mp[1]
    G = mp[2]

    #E = 9*K*G/(3*K+G)
    #v = (3*K-2*G)/(3*K+G)  * (1/2)

    E = 1e0 
    v = 0.3

    F  = zeros(3,3)

    F[1,:] = [eff[1] eff[2] 0.0]
    F[2,:] = [eff[3] eff[4] 0.0]
    F[3,:] = [0.0 0.0 1.0]

    D = E/((1+v)*(1-2*v)) .* [1-v v 0;v 1-v 0; 0 0 0.5*(1-2*v)];

end