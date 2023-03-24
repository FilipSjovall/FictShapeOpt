using StaticArrays
using LinearAlgebra
ngp = 3

ξ      = [2.0/3.0 1.0/3.0 1.0/3.0; 1.0/3.0 2.0/3.0 1.0/3.0; 1.0/3.0 1.0/3.0 2.0/3.0 ]
w      = [1.0/3.0 1.0/3.0 1.0/3.0]
index  = [1 2 3; 4 5 6; 7 8 9] 
#index  = [1 4 7; 2 5 8; 3 6 9] 
P₀     = [0.0 1.0 0.0; 0.0 0.0 1.0]
#P₀     = [0.0 0.0; 1.0 0.0; 0.0 1.0]


Np        = Matrix{Float64}(undef,3,6)
#N = SMatrix{3,6,Float64}

Np[:,1]   = ξ[1,:].*(2.0.*ξ[1,:] .- 1.0)
Np[:,2]   = ξ[2,:].*(2.0.*ξ[2,:] .- 1.0)
Np[:,3]   = ξ[3,:].*(2.0.*ξ[3,:] .- 1.0)
Np[:,4]   = 4.0.*ξ[1,:].*ξ[2,:]
Np[:,5]   = 4.0.*ξ[2,:].*ξ[3,:]
Np[:,6]   = 4.0.*ξ[1,:].*ξ[3,:]
 
N = SMatrix{3,6,Float64}(Np)

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

Jᵀ      = zeros(3,3)
Jᵀ[:,1] = [1.0 1.0 1.0]
Stress  = zeros(2,2)
Bₗ₀      = zeros(3,12)
H₀      = zeros(4,12)
A       = zeros(3,4)
B₀      = zeros(3,12)
R       = zeros(4,4)
eff     = zeros(2,2)



S = zeros(2,2)


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
                t = 1 * getnormal(fv, q_point)
                dΓ = getdetJdV(fv, q_point)
                for i in 1:12#ndofs
                    #δui = shape_value(fv, q_point, i)
                    Ni = shape_value(fv, q_point, i)
                    #ge[i] -= (δui ⋅ t) * dΓ
                    for j in 1:12
                        Nj = shape_value(fv, q_point, j)
                        ge[i]   += Ni ⋅ Nj * (uₑ[i] - λ * dₑ[i]) * dΓ
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
