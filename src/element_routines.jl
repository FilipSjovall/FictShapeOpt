# ###################
# Element routines for a 6-node triangular total Lagrangian element
# ###################
#
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
Bₗ₀           = zeros(3,12)
H₀            = zeros(4,12)
A             = zeros(3,4)
B₀            = zeros(3,12)
R             = zeros(4,4)
#kₑ       = Matrix{Float64}(undef,12,12)

#kₑ       = zeros(12,12)

S = zeros(2,2)

function c2tl6_e!(kₑ,coord,t,D,ed,es)

    for gp ∈ 1:3
        Jᵀ[:,2:3]     = transpose(dNᵣ[:,index[gp,:]]) * coord
        J⁻            = inv(Jᵀ)
        detJ          = det(Jᵀ)
        dNₓ           = P₀ * J⁻ * transpose(dNᵣ[:,index[gp,:]])

        
        @inbounds Bₗ₀[1,1:2:11] = dNₓ[1,:]
        @inbounds Bₗ₀[2,2:2:12] = dNₓ[2,:]  
        @inbounds Bₗ₀[3,1:2:11] = dNₓ[2,:]
        @inbounds Bₗ₀[3,2:2:12] = dNₓ[1,:]

        
        @inbounds H₀[1,1:2:11]  = dNₓ[1,:]
        @inbounds H₀[2,1:2:11]  = dNₓ[2,:]
        @inbounds H₀[3,2:2:12]  = dNₓ[1,:]
        @inbounds H₀[4,2:2:12]  = dNₓ[2,:]

        A_temp = H₀*ed
        
        @inbounds A[1,:]        = [A_temp[1] 0.0 A_temp[3] 0.0]
        @inbounds A[2,:]        = [0.0 A_temp[2] 0.0 A_temp[4]]
        @inbounds A[3,:]        = [A_temp[2] A_temp[1] A_temp[4] A_temp[3]]

        
        @inbounds B₀            = Bₗ₀ + A*H₀
        @inbounds Stress[1,:]        = [es[1,gp] es[3,gp]]
        @inbounds Stress[2,:]        = [es[3,gp] es[2,gp]]

        
        @inbounds R[1:2,1:2]    = Stress
        @inbounds R[3:4,3:4]    = Stress
        @inbounds kₑ[:,:]      += (transpose(B₀)*D[:,:,gp]*B₀ + transpose(H₀)*R*H₀)*detJ*t*w[gp]/2 

    end
end


function c2tl6_f!(fₑ,coord,t,es,ed)
    @simd for gp ∈ 1:3
        Jᵀ[:,2:3]     = transpose(dNᵣ[:,index[gp,:]]) * coord
        J⁻            = inv(Jᵀ)
        detJ          = det(Jᵀ)
        dNₓ           = P₀ * J⁻ * transpose(dNᵣ[:,index[gp,:]])

        
        @inbounds Bₗ₀[1,1:2:11] = dNₓ[1,:]
        @inbounds Bₗ₀[2,2:2:12] = dNₓ[2,:]  
        @inbounds Bₗ₀[3,1:2:11] = dNₓ[2,:]
        @inbounds Bₗ₀[3,2:2:12] = dNₓ[1,:]

        
        @inbounds H₀[1,1:2:11]  = dNₓ[1,:]
        @inbounds H₀[2,1:2:11]  = dNₓ[2,:]
        @inbounds H₀[3,2:2:12]  = dNₓ[1,:]
        @inbounds H₀[4,2:2:12]  = dNₓ[2,:]

        A_temp = H₀*ed
        @inbounds A[1,:]        = [A_temp[1] 0.0 A_temp[3] 0.0]
        @inbounds A[2,:]        = [0.0 A_temp[2] 0.0 A_temp[4]]
        @inbounds A[3,:]        = [A_temp[2] A_temp[1] A_temp[4] A_temp[3]]

        B₀            = Bₗ₀ + A*H₀
        S             = [es[1,gp]; es[2,gp]; es[3,gp]]

        fₑ[:]         += transpose(B₀)*S*detJ*t*w[gp]/2
    end
end

function c2tl6_d!(eff,ed,coord)
    for gp ∈ 1:3
        @inbounds Jᵀ[:,2:3]     = transpose(dNᵣ[:,index[gp,:]]) * coord
        J⁻            = inv(Jᵀ)
        dNₓ           = P₀ * J⁻ * transpose(dNᵣ[:,index[gp,:]])

        
        @inbounds H₀[1,1:2:11]  = dNₓ[1,:]
        @inbounds H₀[2,1:2:11]  = dNₓ[2,:]
        @inbounds H₀[3,2:2:12]  = dNₓ[1,:]
        @inbounds H₀[4,2:2:12]  = dNₓ[2,:]

        temp        = H₀*ed
        eff[gp,1,1] = temp[1] + 1.0
        eff[gp,1,2] = temp[2]
        eff[gp,2,1] = temp[3]
        eff[gp,2,2] = temp[4] + 1.0
    end
end
