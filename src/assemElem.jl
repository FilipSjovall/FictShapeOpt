using StaticArrays
using LinearAlgebra


## Gör om allt till Tensorer?

ngp = 3

ξ      = [2.0/3.0 1.0/3.0 1.0/3.0; 1.0/3.0 2.0/3.0 1.0/3.0; 1.0/3.0 1.0/3.0 2.0/3.0 ]
w      = [1.0/3.0 1.0/3.0 1.0/3.0]
index  = [1 2 3; 4 5 6; 7 8 9] 
#index  = [1 4 7; 2 5 8; 3 6 9] 
P₀     = [0.0 1.0 0.0; 0.0 0.0 1.0]
#P₀     = [0.0 0.0; 1.0 0.0; 0.0 1.0]


Np        = Matrix{Real}(undef,3,6)
#N = SMatrix{3,6,Float64}

Np[:,1]   = ξ[1,:].*(2.0.*ξ[1,:] .- 1.0)
Np[:,2]   = ξ[2,:].*(2.0.*ξ[2,:] .- 1.0)
Np[:,3]   = ξ[3,:].*(2.0.*ξ[3,:] .- 1.0)
Np[:,4]   = 4.0.*ξ[1,:].*ξ[2,:]
Np[:,5]   = 4.0.*ξ[2,:].*ξ[3,:]
Np[:,6]   = 4.0.*ξ[1,:].*ξ[3,:]
 
N = Np
#N = SMatrix{3,6,Float64}(Np)

SdNᵣ      = Matrix{Float64}(undef,6,9)
#dNᵣ      = @SMatrix zeros(6,9)
SdNᵣ[:,1] = [4.0*ξ[1,1]-1.0 0.0 0.0 4.0*ξ[1,2] 0.0 4.0*ξ[1,3]]
SdNᵣ[:,2] = [0.0 4.0*ξ[1,2]-1.0 0.0 4.0*ξ[1,1] 4.0*ξ[1,3] 0.0]
SdNᵣ[:,3] = [0.0 0.0 4.0*ξ[1,3]-1.0 0.0 4.0*ξ[1,2] 4.0*ξ[1,1]]

SdNᵣ[:,4] = [4.0*ξ[2,1]-1.0 0.0 0.0 4.0*ξ[2,2] 0.0 4.0*ξ[2,3]]
SdNᵣ[:,5] = [0.0 4.0*ξ[2,2]-1.0 0.0 4.0*ξ[2,1] 4.0*ξ[2,3] 0.0]
SdNᵣ[:,6] = [0.0 0.0 4.0*ξ[2,3]-1.0 0.0 4.0*ξ[2,2] 4.0*ξ[2,1]]

SdNᵣ[:,7] = [4.0*ξ[3,1]-1.0 0.0 0.0 4.0*ξ[3,2] 0.0 4.0*ξ[3,3]]
SdNᵣ[:,8] = [0.0 4.0*ξ[3,2]-1.0 0.0 4.0*ξ[3,1] 4.0*ξ[3,3] 0.0]
SdNᵣ[:,9] = [0.0 0.0 4.0*ξ[3,3]-1.0 0.0 4.0*ξ[3,2] 4.0*ξ[3,1]]
dNᵣ = SMatrix{6,9,Float64}(SdNᵣ)
#dNᵣ = SdNᵣ


Jᵀ      = zeros(3,3)
Jᵀ[:,1] = [1.0 1.0 1.0]

Stress  = zeros(2,2)
Bₗ₀     = zeros(3,12)
H₀      = zeros(4,12)
A       = zeros(3,4)
B₀      = zeros(3,12)

dBₗ₀    = zeros(3,12)
dH₀     = zeros(4,12)
dA      = zeros(3,4)
dB₀     = zeros(3,12)

R       = zeros(4,4)
eff     = zeros(2,2)



#S       = zeros(2,2)
S       = zeros(3)
∂S_∂x   = zeros(3)


dre      = zeros(12,12)

function assemGP(coord,ed,gp,mp,t)
    
    @inbounds Jᵀ[:,2:3]     = transpose(dNᵣ[:,index[gp,:]]) * coord
    J⁻                      = inv(Jᵀ)
    detJ                    = det(Jᵀ)
    dNₓ                     = P₀ * J⁻ * transpose(dNᵣ[:,index[gp,:]])

    # Gradient matrices  
    @inbounds H₀[1,1:2:11]  = dNₓ[1,:]
    @inbounds H₀[2,1:2:11]  = dNₓ[2,:]
    @inbounds H₀[3,2:2:12]  = dNₓ[1,:]
    @inbounds H₀[4,2:2:12]  = dNₓ[2,:]

    @inbounds Bₗ₀[1,1:2:11] = dNₓ[1,:]
    @inbounds Bₗ₀[2,2:2:12] = dNₓ[2,:]  
    @inbounds Bₗ₀[3,1:2:11] = dNₓ[2,:]
    @inbounds Bₗ₀[3,2:2:12] = dNₓ[1,:]

    A_temp = H₀*ed
    @inbounds A[1,:]        = [A_temp[1] 0.0 A_temp[3] 0.0]
    @inbounds A[2,:]        = [0.0 A_temp[2] 0.0 A_temp[4]]
    @inbounds A[3,:]        = [A_temp[2] A_temp[1] A_temp[4] A_temp[3]]

    B₀                      = Bₗ₀ + A*H₀


    
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

    @inbounds Stress[1,:]   = [S[1] S[3]]
    @inbounds Stress[2,:]   = [S[3] S[2]]
    
    @inbounds R[1:2,1:2]    = Stress
    @inbounds R[3:4,3:4]    = Stress
    
    fₑ            = transpose(B₀)*S*detJ*t*w[gp]/2
    kₑ            = (transpose(B₀)*D*B₀ + transpose(H₀)*R*H₀)*detJ*t*w[gp]/2 
    
    return kₑ, fₑ
end

function RobinIntegral(ke,ge,cell,ΓN,fv,uₑ,λ,dₑ)
    for face in 1:nfaces(cell)
        if (cellid(cell), face) in ΓN
            Ferrite.reinit!(fv, cell, face)
            for q_point in 1:getnquadpoints(fv)
                dΓ = getdetJdV(fv, q_point)
                for i in 1:12#ndofs
                    #δui = shape_value(fv, q_point, i)
                    Ni  = shape_value(fv, q_point, i)
                    u_n = function_value(fv,q_point,uₑ)
                    d_n = function_value(fv,q_point,dₑ)
                    ge[i]   += Ni ⋅ (u_n - λ * d_n) * dΓ
                    for j in 1:12
                        Nj = shape_value(fv, q_point, j)
                        #ge[i]   += Ni ⋅ Nj * (uₑ[j] - λ * dₑ[j]) * dΓ # Borde vara *t på samtliga integraler
                        ke[i,j] += Ni ⋅ Nj * dΓ
                    end
                end
            end
        end
    end
    return ke,ge
end

function assemElem(coord,ed,mp,t)
    kₑ = zeros(12,12)
    fₑ = zeros(12)
    for gp ∈ 1 : 3
        ke, fe= assemGP(coord,ed,gp,mp,t)
        kₑ += ke
        fₑ += fe
    end
    return kₑ, fₑ
end

function init_∂X()
    dX = zeros(12,6,2)
    for j in 1:6
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
    @inbounds H₀[1,1:2:11]  = dNₓ[1,:]
    @inbounds H₀[2,1:2:11]  = dNₓ[2,:]
    @inbounds H₀[3,2:2:12]  = dNₓ[1,:]
    @inbounds H₀[4,2:2:12]  = dNₓ[2,:]

    @inbounds Bₗ₀[1,1:2:11] = dNₓ[1,:]
    @inbounds Bₗ₀[2,2:2:12] = dNₓ[2,:]  
    @inbounds Bₗ₀[3,1:2:11] = dNₓ[2,:]
    @inbounds Bₗ₀[3,2:2:12] = dNₓ[1,:]

    A_temp = H₀*ed
    @inbounds A[1,:]        = [A_temp[1] 0.0 A_temp[3] 0.0]
    @inbounds A[2,:]        = [0.0 A_temp[2] 0.0 A_temp[4]]
    @inbounds A[3,:]        = [A_temp[2] A_temp[1] A_temp[4] A_temp[3]]

    B₀                      = Bₗ₀ + A*H₀
    for dof in 1:12
        ## Sensitivities
        ∂G_∂x                   = - dNₓ * dX[dof,:,:] * dNₓ
        ∂J_∂x                   =   detJ  * tr(dNₓ*dX[dof,:,:]) 

        dH₀[1,1:2:11]  = ∂G_∂x[1,:]
        dH₀[2,1:2:11]  = ∂G_∂x[2,:]
        dH₀[3,2:2:12]  = ∂G_∂x[1,:]
        dH₀[4,2:2:12]  = ∂G_∂x[2,:]
        
        dBₗ₀[1,1:2:11] = ∂G_∂x[1,:]
        dBₗ₀[2,2:2:12] = ∂G_∂x[2,:]  
        dBₗ₀[3,1:2:11] = ∂G_∂x[2,:]
        dBₗ₀[3,2:2:12] = ∂G_∂x[1,:]

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

        @inbounds Stress[1,:]   = [S[1] S[3]]
        @inbounds Stress[2,:]   = [S[3] S[2]]
        
        @inbounds R[1:2,1:2]    = Stress
        @inbounds R[3:4,3:4]    = Stress

        #fₑ            = transpose(B₀)*S*detJ*t*w[gp]/2
   #     println("Shapes del 1: ")
   #     println(size(dB₀),size(S))
   #     println(size(transpose(dB₀)*S))

        dre[:,dof]                      = ( transpose(dB₀)*S + transpose(B₀)*∂S_∂x) *detJ*t*w[gp]/2 + transpose(B₀)*S*∂J_∂x*t*w[gp]/2
    end
    return dre
end
