# ###################
# Element routines for a 6-node triangular total Lagrangian element
# ###################
#
ngp = 3

ξ      = [2.0/3.0 1.0/3.0 1.0/3.0; 1.0/3.0 2.0/3.0 1.0/3.0; 1.0/3.0 1.0/3.0 2.0/3.0 ]
w      = [1.0/3.0 1.0/3.0 1.0/3.0]
index  = [1 4 7; 2 5 8; 3 6 9] 
P₀     = [0.0 1.0 0.0; 0.0 0.0 1.0]
#P₀     = [0.0 0.0; 1.0 0.0; 0.0 1.0]


N        = Matrix{Float64}(undef,3,6)
N[:,1]   = ξ[1,:].*(2.0.*ξ[1,:] .- 1.0)
N[:,2]   = ξ[2,:].*(2.0.*ξ[2,:] .- 1.0)
N[:,3]   = ξ[3,:].*(2.0.*ξ[3,:] .- 1.0)
N[:,4]   = 4.0.*ξ[1,:].*ξ[2,:]
N[:,5]   = 4.0.*ξ[2,:].*ξ[3,:]
N[:,6]   = 4.0.*ξ[1,:].*ξ[3,:]
 
dNᵣ      = Matrix{Float64}(undef,6,9)
dNᵣ[:,1] = [4.0*ξ[1,1]-1.0 0.0 0.0 4.0*ξ[1,2] 0.0 4.0*ξ[1,3]]
dNᵣ[:,2] = [0.0 4.0*ξ[1,2]-1.0 0.0 4.0*ξ[1,1] 4.0*ξ[1,3] 0.0]
dNᵣ[:,3] = [0.0 0.0 4.0*ξ[1,3]-1.0 0.0 4.0*ξ[1,2] 4.0*ξ[1,1]]

dNᵣ[:,4] = [4.0*ξ[2,1]-1.0 0.0 0.0 4.0*ξ[2,2] 0.0 4.0*ξ[2,3]]
dNᵣ[:,5] = [0.0 4.0*ξ[2,2]-1.0 0.0 4.0*ξ[2,1] 4.0*ξ[2,3] 0.0]
dNᵣ[:,6] = [0.0 0.0 4.0*ξ[2,3]-1.0 0.0 4.0*ξ[2,2] 4.0*ξ[2,1]]

dNᵣ[:,7] = [4.0*ξ[3,1]-1.0 0.0 0.0 4.0*ξ[3,2] 0.0 4.0*ξ[3,3]]
dNᵣ[:,8] = [0.0 4.0*ξ[3,2]-1.0 0.0 4.0*ξ[3,1] 4.0*ξ[3,3] 0.0]
dNᵣ[:,9] = [0.0 0.0 4.0*ξ[3,3]-1.0 0.0 4.0*ξ[3,2] 4.0*ξ[3,1]]

#kₑ       = Matrix{Float64}(undef,12,12)

#kₑ       = zeros(12,12)

S = zeros(2,2)

function c2tl6_e!(kₑ,coord,t,D,ed,es)
    Jᵀ      = zeros(3,3)
    Jᵀ[:,1] = [1.0 1.0 1.0]
    Stress  = zeros(2,2)
    for gp ∈ 1:3
        Jᵀ[:,2:3]     = transpose(dNᵣ[:,index[gp,:]]) * coord
        J⁻            = inv(Jᵀ)
        detJ          = det(Jᵀ)
        dNₓ           = P₀ * J⁻ * transpose(dNᵣ[:,index[gp,:]])

        Bₗ₀           = zeros(3,12)
        Bₗ₀[1,1:2:11] = dNₓ[1,:]
        Bₗ₀[2,2:2:12] = dNₓ[2,:]  
        Bₗ₀[3,1:2:11] = dNₓ[2,:]
        Bₗ₀[3,2:2:12] = dNₓ[1,:]

        H₀            = zeros(4,12)
        H₀[1,1:2:11]  = dNₓ[1,:]
        H₀[2,1:2:11]  = dNₓ[2,:]
        H₀[3,2:2:12]  = dNₓ[1,:]
        H₀[4,2:2:12]  = dNₓ[2,:]

        A_temp = H₀*ed
        A             = zeros(3,4)
        A[1,:]        = [A_temp[1] 0.0 A_temp[3] 0.0]
        A[2,:]        = [0.0 A_temp[2] 0.0 A_temp[4]]
        A[3,:]        = [A_temp[2] A_temp[1] A_temp[4] A_temp[3]]

        B₀            = zeros(3,12)
        B₀            = Bₗ₀ + A*H₀
        Stress[1,:]        = [es[1,gp] es[3,gp]]
        Stress[2,:]        = [es[3,gp] es[2,gp]]

        R             = zeros(4,4)
        R[1:2,1:2]    = Stress
        R[3:4,3:4]    = Stress
        kₑ            .= kₑ .+ (transpose(B₀)*D[:,:,gp]*B₀ + transpose(H₀)*R*H₀)*detJ*t*w[gp]/2 

    end
    
    #return kₑ
end


function c2tl6_f!(fₑ,coord,t,es,ed)
    Jᵀ      = zeros(3,3)
    Jᵀ[:,1] = [1.0 1.0 1.0]
    
    for gp ∈ 1:3
        Jᵀ[:,2:3]     = transpose(dNᵣ[:,index[gp,:]]) * coord
        J⁻            = inv(Jᵀ)
        detJ          = det(Jᵀ)
        dNₓ           = P₀ * J⁻ * transpose(dNᵣ[:,index[gp,:]])

        Bₗ₀           = zeros(3,12)
        Bₗ₀[1,1:2:11] = dNₓ[1,:]
        Bₗ₀[2,2:2:12] = dNₓ[2,:]  
        Bₗ₀[3,1:2:11] = dNₓ[2,:]
        Bₗ₀[3,2:2:12] = dNₓ[1,:]

        H₀            = zeros(4,12)
        H₀[1,1:2:11]  = dNₓ[1,:]
        H₀[2,1:2:11]  = dNₓ[2,:]
        H₀[3,2:2:12]  = dNₓ[1,:]
        H₀[4,2:2:12]  = dNₓ[2,:]

        A_temp = H₀*ed
        A             = zeros(3,4)
        A[1,:]        = [A_temp[1] 0.0 A_temp[3] 0.0]
        A[2,:]        = [0.0 A_temp[2] 0.0 A_temp[4]]
        A[3,:]        = [A_temp[2] A_temp[1] A_temp[4] A_temp[3]]

        B₀            = zeros(3,12)
        B₀            = Bₗ₀ + A*H₀
        S             = [es[1,gp]; es[2,gp]; es[3,gp]]

        fₑ            .= fₑ.+ transpose(B₀)*S*detJ*t*w[gp]/2
    end
    return fₑ
end

function c2tl6_d!(eff,ed,coord)
    Jᵀ      = zeros(3,3)
    Jᵀ[:,1] = [1.0 1.0 1.0]
    #eff = zeros(3,2,2)
    for gp ∈ 1:3
        Jᵀ[:,2:3]     = transpose(dNᵣ[:,index[gp,:]]) * coord
        J⁻            = inv(Jᵀ)
        dNₓ           = P₀ * J⁻ * transpose(dNᵣ[:,index[gp,:]])

        H₀            = zeros(4,12)
        H₀[1,1:2:11]  = dNₓ[1,:]
        H₀[2,1:2:11]  = dNₓ[2,:]
        H₀[3,2:2:12]  = dNₓ[1,:]
        H₀[4,2:2:12]  = dNₓ[2,:]

        temp            = H₀*ed
        eff[gp,1,1] = temp[1] + 1.0
        eff[gp,1,2] = temp[2]
        eff[gp,2,1] = temp[3]
        eff[gp,2,2] = temp[4] + 1.0
    end
    return eff
end
