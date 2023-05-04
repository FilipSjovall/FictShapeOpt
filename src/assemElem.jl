using StaticArrays
using LinearAlgebra


## Gör om allt till Tensorer?

ngp = 3

#ξ      = [2.0/3.0 1.0/3.0 1.0/3.0; 1.0/3.0 2.0/3.0 1.0/3.0; 1.0/3.0 1.0/3.0 2.0/3.0 ]
ξ      = [2.0/3.0 1.0/6.0 1.0/6.0; 1.0/6.0 2.0/3.0 1.0/6.0; 1.0/6.0 1.0/6.0 2.0/3.0 ]
w      = [1.0/3.0 1.0/3.0 1.0/3.0]
index  = [1 2 3; 4 5 6; 7 8 9] 
#w      = [-9/32 25/96 25/96 25/96]
#ξ      = [1/3 1/3 1-2/3; 3/5 1/5 1-4/5; 1/5 3/5 1-4/5; 1/5 1/5 1-2/5]
#index  = [1 2 3; 4 5 6; 7 8 9; 10 11 12] 

P₀     = [0.0 1.0 0.0; 0.0 0.0 1.0]
#P₀     = [0.0 0.0; 1.0 0.0; 0.0 1.0]


Np        = Matrix{Real}(undef,3,6)
Npf       = Matrix{Real}(undef,3,6)
Nf        = Array{Real}(undef,2,12,3)
#N = SMatrix{3,6,Float64}

Np[:,1]   = ξ[1,:].*(2.0.*ξ[1,:] .- 1.0)
Np[:,2]   = ξ[2,:].*(2.0.*ξ[2,:] .- 1.0)
Np[:,3]   = ξ[3,:].*(2.0.*ξ[3,:] .- 1.0)
Np[:,4]   = 4.0.*ξ[1,:].*ξ[2,:]
Np[:,5]   = 4.0.*ξ[2,:].*ξ[3,:]
Np[:,6]   = 4.0.*ξ[1,:].*ξ[3,:]



N = Np
dNᵣ      = @SMatrix zeros(6,9)
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
    
    Jᵀ[:,2:3]               = transpose(dNᵣ[:,index[gp,:]]) * coord

    J⁻                      = inv(Jᵀ)
    detJ                    = det(Jᵀ)
    dNₓ                     = P₀ * J⁻ * transpose(dNᵣ[:,index[gp,:]])

    # Gradient matrices  
    H₀[1,1:2:11]  = dNₓ[1,:]
    H₀[2,1:2:11]  = dNₓ[2,:]
    H₀[3,2:2:12]  = dNₓ[1,:]
    H₀[4,2:2:12]  = dNₓ[2,:]

    Bₗ₀[1,1:2:11] = dNₓ[1,:]
    Bₗ₀[2,2:2:12] = dNₓ[2,:]  
    Bₗ₀[3,1:2:11] = dNₓ[2,:]
    Bₗ₀[3,2:2:12] = dNₓ[1,:]

    A_temp = H₀*ed

    A[1,:]        = [A_temp[1] 0.0 A_temp[3] 0.0]
    A[2,:]        = [0.0 A_temp[2] 0.0 A_temp[4]]
    A[3,:]        = [A_temp[2] A_temp[1] A_temp[4] A_temp[3]]

    B₀            = Bₗ₀ + A*H₀

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

    # ∫ δEᵀ D ΔE + ∇δuᵀ S ∇Δu dΩ 
    kₑ            = ( transpose(B₀)*D*B₀ + transpose(H₀)*R*H₀ )*detJ*t*w[gp]/2 
    return kₑ, fₑ
end

function RobinIntegral(ke,ge,cell,ΓN,fv,uₑ,λ,dₑ,coorde)
    #ge2 = ge.*0;
    #ke2 = ke.*0;
    for face in 1:nfaces(cell)
        if (cellid(cell), face) in ΓN
            Ferrite.reinit!(fv, cell, face)
            for q_point in 1:getnquadpoints(fv)
                dΓ = getdetJdV(fv, q_point)
                u_n = function_value(fv,q_point,uₑ)
                d_n = function_value(fv,q_point,dₑ)
                if (cellid(cell), face) in Γ1
                    for i in 1:2:11
                        Ni  = shape_value(fv, q_point, i)
                        #ge[i]   += Ni ⋅ (u_n - λ * d_n) * dΓ 
                        ge[i]   += Ni ⋅ [1;1] 
                        for j in 1:2:11
                            Nj = shape_value(fv, q_point, j)
                            ke[i,j] += Ni ⋅ Nj * dΓ 
                        end
                    end
                elseif (cellid(cell), face) in Γ2
                    for i in 2:2:12
                        Ni  = shape_value(fv, q_point, i)
                        #ge[i]   += Ni ⋅ (u_n - λ * d_n) * dΓ 
                        ge[i]   += Ni ⋅ [1;1]  
                        for j in 2:2:12
                            Nj = shape_value(fv, q_point, j)
                            ke[i,j] += Ni ⋅ Nj * dΓ 
                        end
                    end
                end 
            end
            #if (cellid(cell), face) in Γ1 
            #    for gp in 1:2
            #        @inbounds Jᵀ[:,2:3]     = transpose(dNf[:,index[gp,:]]) * coorde
            #        J⁻                      = inv(Jᵀ)
            #        detJ                    = det(Jᵀ)
            #        dΓ                      = w[gp] * detJ /2
            #        #@inbounds ge[1:2:11]             += Nf[:,1:2:11,gp]' * (Nf[:,1:2:11,gp]*uₑ[1:2:11] - λ * Nf[:,1:2:11,gp] * dₑ[1:2:11]) * dΓ
            #        @inbounds ke[1:2:11,1:2:11]      += Nf[:,1:2:11,gp]' * Nf[:,1:2:11,gp] * dΓ
            #    end
            #elseif  (cellid(cell), face) in Γ2
            #    for gp in 1:2
            #        @inbounds Jᵀ[:,2:3]     = transpose(dNf[:,index[gp,:]]) * coorde
            #        J⁻                      = inv(Jᵀ)
            #        detJ                    = det(Jᵀ)
            #        dΓ                      = w[gp] * detJ /2
            #        #@inbounds ge[2:2:end]            += Nf[:,2:2:end,gp]' * (Nf[:,2:2:end,gp]*uₑ[2:2:12] - λ * Nf[:,2:2:end,gp] * dₑ[2:2:12]) * dΓ
            #        @inbounds ge            += Nf[:,:,gp]' * [1;1] 
            #        println(ge)
            #        @inbounds ke[2:2:end,2:2:end]    += Nf[:,2:2:end,gp]' * Nf[:,2:2:end,gp] * dΓ
            #    end
            #end
            #println("ke = ", ke)
            #println("ge = ", ge)
        end
    end
    #println("ge = ", ge)
    return ke,ge
end

function assemElem(coord,ed,mp,t)
    kₑ = zeros(12,12)
    fₑ = zeros(12)
    for gp ∈ 1 : ngp
        ke, fe= assemGP(coord,ed,gp,mp,t)
        kₑ += ke
        fₑ += fe
    end
    return kₑ, fₑ
end

function assem_dr(coord,ed,mp,t)
    drₑ = zeros(12,12)
    for gp ∈ 1 : ngp
        dre  = dr_GP(coord,ed,gp,mp,t)
        drₑ += dre
    end
    return drₑ
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
    dX = init_∂X();
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


        dre[:,dof]              = ( transpose(dB₀)*S + transpose(B₀)*∂S_∂x ) *detJ*t*w[gp]/2 + transpose(B₀)*S*∂J_∂x*t*w[gp]/2
    end
    return dre
end

function Robin(coorde,ue,de,λ)
    L1 = norm(coorde[1,:] - coorde[3,:])
    L2 = norm(coorde[1,:] - coorde[2,:])
    Kc   = zeros(6,6)

    #N1N1 =L2/3 - ((L1*L2^2)/6 - L2^3/30)/L1^2
    #N1N2 = - (L2^3*(5*L1 - 2*L2))/(60*L1^2*(L1 - L2))
    #N1N3 = L2/6 + L2^3/(20*L1*(L1 - L2))
    #N2N2 = L2^5/(30*L1^2*(L1 - L2)^2)
    #N2N3 =  - (L2^3*(5*L1 - 3*L2))/(60*L1*(L1 - L2)^2)
    #N3N3 = L2/3 + ((L1*L2^2)/6 - (2*L2^3)/15)/(L1 - L2)^2

    #Kc[1,1:2:5] = [N1N1 N1N2 N1N3]
    #Kc[2,2:2:6] = [N1N1 N1N2 N1N3]
    #Kc[3,1:2:5] = [N1N3 N2N3 N3N3]
    #Kc[4,2:2:6] = [N1N3 N2N3 N3N3]#
    #Kc[5,1:2:5] = [N1N2 N2N2 N2N3]
    #Kc[6,2:2:6] = [N1N2 N2N2 N2N3]

    N1N1 = L2/3 - ((L1*L2^2)/6 - L2^3/30)/L1^2
    N1N3 = - (L2^3*(5*L1 - 2*L2))/(60*L1^2*(L1 - L2))
    N1N2 = L2/6 + L2^3/(20*L1*(L1 - L2))
    N3N3 = L2^5/(30*L1^2*(L1 - L2)^2)
    N2N3 =  - (L2^3*(5*L1 - 3*L2))/(60*L1*(L1 - L2)^2)
    N2N2 = L2/3 + ((L1*L2^2)/6 - (2*L2^3)/15)/(L1 - L2)^2
    Kc[1,1:2:5] = [N1N1 N1N2 N1N3]
    Kc[2,2:2:6] = [N1N1 N1N2 N1N3]
    Kc[3,1:2:5] = [N1N2 N2N2 N2N3]
    Kc[4,2:2:6] = [N1N2 N2N2 N2N3]
    Kc[5,1:2:5] = [N1N3 N2N3 N3N3]
    Kc[6,2:2:6] = [N1N3 N2N3 N3N3]


    println("Stämmer detta om randen är krokig????")
    #∫N1         = L2/2 - L2^2/(6*L1)
    #∫N2         = -L2^3/(6*L1*(L1 - L2))
    #∫N3         = L2/2 + L2^2/(6*(L1 - L2)) 
    #println(∫N1,∫N2,∫N3)

return  Kc, Kc*(ue-λ*de) # [0;∫N1;0;∫N2;0;∫N3]*0.5
end

function dΩ(X)
    dVol = 0.0
    for gp ∈ 1 : 3
        Jᵀ[:,2:3]   = transpose(dNᵣ[:,index[gp,:]]) * X
        detJ        = det(Jᵀ)
        dVol       += detJ * w[gp]*t/2
    end
    return dVol
end

function ∂dΩ∂x(X)
    dVole = zeros(12)
    for gp ∈ 1 : 3
        @inbounds Jᵀ[:,2:3]     = transpose(dNᵣ[:,index[gp,:]]) * X
        detJ                    = det(Jᵀ)
        J⁻                      = inv(Jᵀ)
        dNₓ                     = P₀ * J⁻ * transpose(dNᵣ[:,index[gp,:]])
        dX = init_∂X();
        for dof in 1:12
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
        #println(∂dΩ∂x(coord[enod[ie][2:7],:]))
        dVol[cell_dofs]   += ∂dΩ∂x(coord[enod[ie][2:7],:])
    end
    return dVol
end