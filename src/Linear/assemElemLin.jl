using StaticArrays
using LinearAlgebra
using Tensors

## Gör om allt till Tensorer?

ngp = 3

ξ      = [2.0/3.0 1.0/6.0 1.0/6.0; 1.0/6.0 2.0/3.0 1.0/6.0; 1.0/6.0 1.0/6.0 2.0/3.0 ]
w      = [1.0/3.0 1.0/3.0 1.0/3.0]
index  = [1 2 3; 4 5 6; 7 8 9] 

P₀     = [0.0 1.0 0.0; 0.0 0.0 1.0]


N        = Matrix{Float64}(undef,3,3)
#N = SMatrix{3,6,Float64}

N[:,1]   = ξ[1,:]
N[:,2]   = ξ[2,:]
N[:,3]   = ξ[3,:]

#M        = reinterpret(Vec{3,Float64},vec(N))
dNᵣ      = Matrix{Float64}(undef,3,9)
dNᵣ[:,1] = [1.0 0.0 0.0]
dNᵣ[:,2] = [0.0 1.0 0.0]
dNᵣ[:,3] = [0.0 0.0 1.0]
#
dNᵣ[:,4] = [1.0 0.0 0.0]
dNᵣ[:,5] = [0.0 1.0 0.0]
dNᵣ[:,6] = [0.0 0.0 1.0]
#
dNᵣ[:,7] = [1.0 0.0 0.0]
dNᵣ[:,8] = [0.0 1.0 0.0]
dNᵣ[:,9] = [0.0 0.0 1.0]
#dNr      = reinterpret(Vec{3,Float64},vec(dNᵣ))

Jᵀ      = zeros(3,3)
Jᵀ[:,1] = [1.0 1.0 1.0]

Stress  = zeros(2,2)
Bₗ₀     = zeros(3,6)
H₀      = zeros(4,6)
A       = zeros(3,4)
B₀      = zeros(3,6)
dBₗ₀    = zeros(3,6)
dH₀     = zeros(4,6)
dA      = zeros(3,4)
dB₀     = zeros(3,6)
R       = zeros(4,4)
eff     = zeros(2,2)
S       = zeros(3)
∂S_∂x   = zeros(3)
dre     = zeros(6,6)


function assemGP(coord,ed,gp,mp,t)
    # Jacobian N x = Jᵀ N ξ
    Jᵀ[:,2:3]               = transpose(dNᵣ[:,index[gp,:]]) * coord # ??
    J⁻                      = inv(Jᵀ)
    detJ                    = det(Jᵀ)
    dNₓ                     = P₀ * J⁻ * transpose(dNᵣ[:,index[gp,:]])

    # Gradient matrices " ∇N " 
    H₀[1,1:2:5]  = dNₓ[1,:]
    H₀[2,1:2:5]  = dNₓ[2,:]
    H₀[3,2:2:6]  = dNₓ[1,:]
    H₀[4,2:2:6]  = dNₓ[2,:]

    Bₗ₀[1,1:2:5] = dNₓ[1,:]
    Bₗ₀[2,2:2:6] = dNₓ[2,:]  
    Bₗ₀[3,1:2:5] = dNₓ[2,:]
    Bₗ₀[3,2:2:6] = dNₓ[1,:]

    # ∇u = ∇ₓN u
    A_temp       = H₀*ed
    A[1,:]       = [A_temp[1] 0.0 A_temp[3] 0.0]
    A[2,:]       = [0.0 A_temp[2] 0.0 A_temp[4]]
    A[3,:]       = [A_temp[2] A_temp[1] A_temp[4] A_temp[3]]

    # ∇N + ∂N∂x ⋅ ∇u 
    B₀           = Bₗ₀ + A*H₀

    # Deformation gradient F = I + ∇u
    eff[1,1]     = A_temp[1] + 1.0
    eff[1,2]     = A_temp[2]
    eff[2,1]     = A_temp[3]
    eff[2,2]     = A_temp[4] + 1.0

    # Vec(F)
    ef = [eff[1,1] eff[1,2] eff[2,1] eff[2,2]]

    # Stress and material tangent: S = 0.5 ∂W / ∂C, D = 0.25 ∂²W / ∂C²    
    es = neohooke1(ef,mp)
    D  = dneohooke1(ef,mp)

    # Reformulate as matrix
    S             = [es[1]; es[2]; es[4]]
    Stress[1,:]   = [S[1] S[3]]
    Stress[2,:]   = [S[3] S[2]]
    R[1:2,1:2]    = Stress
    R[3:4,3:4]    = Stress

    # Internal force vector ∫ δE : S dΩ or ∫ ∇ₓδuᵢ : P dΩ
    fₑ            = transpose(B₀)*S*detJ*t*w[gp]/2
    # ∫ δEᵀ D ΔE + ∇δuᵀ S ∇Δu dΩ 
    kₑ            = ( transpose(B₀)*D*B₀ + transpose(H₀)*R*H₀ )*detJ*t*w[gp]/2 
    return kₑ, fₑ
end

function assemElem(coord,ed,mp,t)
    kₑ = zeros(6,6)
    fₑ = zeros(6)
    for gp ∈ 1 : ngp
        ke, fe= assemGP(coord,ed,gp,mp,t)
        kₑ += ke
        fₑ += fe
    end
    return kₑ, fₑ
end

function assemS(coord, ed, mp, t, gp)
    # Jacobian N x = Jᵀ N ξ
    Jᵀ[:, 2:3] = transpose(dNᵣ[:, index[gp, :]]) * coord # ??
    J⁻ = inv(Jᵀ)
    dNₓ = P₀ * J⁻ * transpose(dNᵣ[:, index[gp, :]])

    # Gradient matrices " ∇N " 
    H₀[1, 1:2:5] = dNₓ[1, :]
    H₀[2, 1:2:5] = dNₓ[2, :]
    H₀[3, 2:2:6] = dNₓ[1, :]
    H₀[4, 2:2:6] = dNₓ[2, :]

    Bₗ₀[1, 1:2:5] = dNₓ[1, :]
    Bₗ₀[2, 2:2:6] = dNₓ[2, :]
    Bₗ₀[3, 1:2:5] = dNₓ[2, :]
    Bₗ₀[3, 2:2:6] = dNₓ[1, :]

    # ∇u = ∇ₓN u
    A_temp = H₀ * ed
    A[1, :] = [A_temp[1] 0.0 A_temp[3] 0.0]
    A[2, :] = [0.0 A_temp[2] 0.0 A_temp[4]]
    A[3, :] = [A_temp[2] A_temp[1] A_temp[4] A_temp[3]]

    # Deformation gradient F = I + ∇u
    eff[1, 1] = A_temp[1] + 1.0
    eff[1, 2] = A_temp[2]
    eff[2, 1] = A_temp[3]
    eff[2, 2] = A_temp[4] + 1.0

    # Vec(F)
    ef = [eff[1, 1] eff[1, 2] eff[2, 1] eff[2, 2]]
    # Stress and material tangent: S = 0.5 ∂W / ∂C, D = 0.25 ∂²W / ∂C²    
    es = neohooke1(ef, mp)
    F  = zeros(3,3)
    Sₑ = zeros(3,3)
    F[1, :] = [eff[1] eff[2] 0.0]
    F[2, :] = [eff[3] eff[4] 0.0]
    F[3, :] = [0.0 0.0 1.0]
    Sₑ[1,:] = [es[1] es[4] 0.0  ]
    Sₑ[2,:] = [es[4] es[2] 0.0  ]
    Sₑ[3,:] = [0.0   0.0   es[3]]
    J = det(F)
    σ_temp = 1 / J * F * Sₑ * F'
    return σ_temp
end

function assem_dr(coord,ed,mp,t)
    drₑ = zeros(6,6)
    for gp ∈ 1 : ngp
        dre  = dr_GP(coord,ed,gp,mp,t)
        drₑ += dre
    end
    return drₑ
end

function init_∂X()
    dX = zeros(6,3,2)
    for j in 1:3
        for i in 1:2
            dof = (j-1)*2 + i
            dX[dof,j,i] = 1.0
        end
    end
    return dX
end

function dr_GP(coord,ed,gp,mp,t)
    @inbounds Jᵀ[:,2:3]     = transpose(dNᵣ[:,index[gp,:]]) * coord
    J⁻                      = inv(Jᵀ)
    detJ                    = det(Jᵀ)
    dNₓ                     = P₀ * J⁻ * transpose(dNᵣ[:,index[gp,:]])

    # Gradient matrices  
    @inbounds H₀[1,1:2:5]  = dNₓ[1,:]
    @inbounds H₀[2,1:2:5]  = dNₓ[2,:]
    @inbounds H₀[3,2:2:6]  = dNₓ[1,:]
    @inbounds H₀[4,2:2:6]  = dNₓ[2,:]

    @inbounds Bₗ₀[1,1:2:5] = dNₓ[1,:]
    @inbounds Bₗ₀[2,2:2:6] = dNₓ[2,:]  
    @inbounds Bₗ₀[3,1:2:5] = dNₓ[2,:]
    @inbounds Bₗ₀[3,2:2:6] = dNₓ[1,:]
    

    A_temp = H₀*ed
    @inbounds A[1,:]        = [A_temp[1] 0.0 A_temp[3] 0.0]
    @inbounds A[2,:]        = [0.0 A_temp[2] 0.0 A_temp[4]]
    @inbounds A[3,:]        = [A_temp[2] A_temp[1] A_temp[4] A_temp[3]]

    B₀                      = Bₗ₀ + A*H₀
    dX = init_∂X();
    for dof in 1:6
        ## Sensitivities
        ∂G_∂x                   = - dNₓ * dX[dof,:,:] * dNₓ
        ∂J_∂x                   =   detJ  * tr(dNₓ*dX[dof,:,:]) 

        dH₀[1,1:2:5]  = ∂G_∂x[1,:]
        dH₀[2,1:2:5]  = ∂G_∂x[2,:]
        dH₀[3,2:2:6]  = ∂G_∂x[1,:]
        dH₀[4,2:2:6]  = ∂G_∂x[2,:]
        
        dBₗ₀[1,1:2:5] = ∂G_∂x[1,:]
        dBₗ₀[2,2:2:6] = ∂G_∂x[2,:]  
        dBₗ₀[3,1:2:5] = ∂G_∂x[2,:]
        dBₗ₀[3,2:2:6] = ∂G_∂x[1,:]


        dA_temp        = dH₀ * ed

        dA[1,:]        = [dA_temp[1] 0.0 dA_temp[3] 0.0]
        dA[2,:]        = [0.0 dA_temp[2] 0.0 dA_temp[4]]
        dA[3,:]        = [dA_temp[2] dA_temp[1] dA_temp[4] dA_temp[3]]
        
        dB₀            = dBₗ₀ + dA * H₀ + A * dH₀ 

        ## Material response
        #
        # Deformation gradient
        eff[1,1] = A_temp[1] + 1.0
        eff[1,2] = A_temp[2]
        eff[2,1] = A_temp[3] 
        eff[2,2] = A_temp[4] + 1.0
 
        ef = [eff[1,1] eff[1,2] eff[2,1] eff[2,2]]

        # Stress and material tangent
        es = neohooke1(ef,mp)
        D  = dneohooke1(ef,mp)

        S                       = [es[1]; es[2]; es[4]]
        ∂S_∂x                   = D * (dBₗ₀ + 0.5 * ( dA * H₀ + A * dH₀ ) ) * ed

        Stress[1,:]   = [S[1] S[3]]
        Stress[2,:]   = [S[3] S[2]]
        R[1:2,1:2]    = Stress
        R[3:4,3:4]    = Stress


        dre[:,dof]              = ( transpose(dB₀)*S + transpose(B₀)*∂S_∂x ) *detJ*t*w[gp]/2 + transpose(B₀)*S*∂J_∂x*t*w[gp]/2
    end
    return dre
end

function Robin(coorde,Ψe,de,λ)
    L           = norm(coorde[1,:] - coorde[2,:])
    Kc          = zeros(4,4)

    # Byt mot 4x4
    N1N1        =  L/3
    N1N2        =  L/6
    N2N2        =  L/3 
    Kc[1,1:2:3] = [N1N1 N1N2]
    Kc[2,2:2:4] = [N1N1 N1N2]
    Kc[3,1:2:3] = [N1N2 N2N2]
    Kc[4,2:2:4] = [N1N2 N2N2]
    #Kc = L/6 *  [2 0 1 0; 0 2 0 1; 1 0 2 0; 0 1 0 2] # - snabbast implementering? inga allokeringar vid sidan...
    return  Kc, Kc * (Ψe - λ*de) # [0;∫N1;0;∫N2;0;∫N3]*0.5
end

function tractionLoad(coorde,τ)
    L = sqrt((coorde[1,1] - coorde[2,1])^2 + (coorde[1,2] - coorde[2,2])^2) 
    ∫Nᵀ = (L/2) * [1.0 0.0 1.0 0.0; 0.0 1.0 0.0 1.0]  
    return -∫Nᵀ' * τ 
    #return -[τ[1];τ[1];τ[1];τ[1]]
end

function dTractionLoad(coorde,τ)
    L  = sqrt((coorde[1,1] - coorde[2,1])^2 + (coorde[1,2] - coorde[2,2])^2) 
    Δx = ( coorde[1,1] - coorde[2,1] ) / (2L) 
    Δy = ( coorde[1,2] - coorde[2,2] ) / (2L) 

    return -[ Δx*τ[1]  Δx*τ[2]  Δx*τ[1]  Δx*τ[2];
              Δy*τ[1]  Δy*τ[2]  Δy*τ[1]  Δy*τ[2];
             -Δx*τ[1] -Δx*τ[2] -Δx*τ[1] -Δx*τ[2];
             -Δy*τ[1] -Δy*τ[2] -Δy*τ[1] -Δy*τ[2]]'
end

function dΩ(X)
    dVol = 0.0
    for gp ∈ 1 : ngp
        Jᵀ[:,2:3]   = transpose(dNᵣ[:,index[gp,:]]) * X
        detJ        = det(Jᵀ)
        dVol       += detJ * w[gp]*t/2
    end
    return dVol
end

function ∂dΩ∂x(X)
    dVole = zeros(6)
    for gp ∈ 1 : ngp
        @inbounds Jᵀ[:,2:3]     = transpose(dNᵣ[:,index[gp,:]]) * X
        detJ                    = det(Jᵀ)
        J⁻                      = inv(Jᵀ)
        dNₓ                     = P₀ * J⁻ * transpose(dNᵣ[:,index[gp,:]])
        dX = init_∂X();
        for dof in 1:6
            dVole[dof] += detJ * tr(dNₓ*dX[dof,:,:]) * w[gp]*t/2
        end
    end
    return dVole
end

function volume_sens(dh,coord)
    dVol= zeros(length(coord))
    ie  = 0
    for cell in CellIterator(dh)
        ie     += 1
        cell_dofs = celldofs(cell)
        dVol[cell_dofs]   += ∂dΩ∂x(coord[enod[ie][2:end],:])
    end
    return dVol
end

function defgradGP(coord,ed,gp,mp,t)
    Jᵀ[:, 2:3] = transpose(dNᵣ[:, index[gp, :]]) * coord # ??
    J⁻ = inv(Jᵀ)
    detJ = det(Jᵀ)
    dNₓ = P₀ * J⁻ * transpose(dNᵣ[:, index[gp, :]])

    # Gradient matrices " ∇N " 
    H₀[1, 1:2:5] = dNₓ[1, :]
    H₀[2, 1:2:5] = dNₓ[2, :]
    H₀[3, 2:2:6] = dNₓ[1, :]
    H₀[4, 2:2:6] = dNₓ[2, :]

    Bₗ₀[1, 1:2:5] = dNₓ[1, :]
    Bₗ₀[2, 2:2:6] = dNₓ[2, :]
    Bₗ₀[3, 1:2:5] = dNₓ[2, :]
    Bₗ₀[3, 2:2:6] = dNₓ[1, :]

    # ∇u = ∇ₓN u
    A_temp = H₀ * ed
    A[1, :] = [A_temp[1] 0.0 A_temp[3] 0.0]
    A[2, :] = [0.0 A_temp[2] 0.0 A_temp[4]]
    A[3, :] = [A_temp[2] A_temp[1] A_temp[4] A_temp[3]]

    # ∇N + ∂N∂x ⋅ ∇u 
    B₀ = Bₗ₀ + A * H₀

    # Deformation gradient F = I + ∇u
    eff[1, 1] = A_temp[1] + 1.0
    eff[1, 2] = A_temp[2]
    eff[2, 1] = A_temp[3]
    eff[2, 2] = A_temp[4] + 1.0
    # Vec(F)
    ef = [eff[1, 1] eff[1, 2] eff[2, 1] eff[2, 2]]

    # Stress and material tangent: S = 0.5 ∂W / ∂C, D = 0.25 ∂²W / ∂C²    

end