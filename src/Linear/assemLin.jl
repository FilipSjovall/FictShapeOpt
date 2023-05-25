# ------------------------ # 
# Assemble global matrices #
# ------------------------ #

function assemGlobal!(K,Fᵢₙₜ,dh,mp,t,a,coord,enod)
    assembler = start_assemble(K,Fᵢₙₜ)
    ie = 0
    kₑ = zeros(6,6)
    fₑ = zeros(6)
    for cell in CellIterator(dh)
        fill!(kₑ,0.0)
        fill!(fₑ,0.0)
        ie += 1
        cell_dofs= celldofs(cell)
        kₑ, fₑ = assemElem(coord[enod[ie][2:end],:],a[cell_dofs],mp,t)
        assemble!(assembler, cell_dofs, kₑ, fₑ)
    end            
end

function assemGlobal!(K,Fᵢₙₜ,dh,mp,t,a,coord,enod,Γt,τ)
    assembler = start_assemble(K,Fᵢₙₜ)
    ie        = 0
    kₑ        = zeros(6,6)
    fₑ        = zeros(6)
    for cell in CellIterator(dh)
        fill!(kₑ,0.0)
        fill!(fₑ,0.0)
        ie       += 1
        cell_dofs = celldofs(cell)
        kₑ, fₑ    = assemElem(coord[enod[ie][2:end],:],a[cell_dofs],mp,t)
        for face in 1:nfaces(cell)
            if (cellid(cell), face) in Γt
                face_nods      = [Ferrite.facedof_indices(ip)[face][1]; Ferrite.facedof_indices(ip)[face][2]]
                face_dofs      = [face_nods[1]*2-1; face_nods[1]*2; face_nods[2]*2-1; face_nods[2]*2]
                X              = coord[enod[ie][face_nods.+1] ,:]
                fₑ[face_dofs] += tractionLoad(X,τ)
            end    
        end
        assemble!(assembler, cell_dofs, kₑ, fₑ)
    end            
end

function assemGlobal!(K,Fᵢₙₜ,dh,mp,t,Ψ,coord,enod,λ,d,Γ_robin)
    assembler = start_assemble(K,Fᵢₙₜ)
    ie = 0
    for cell in CellIterator(dh)
        ie += 1
        cell_dofs = celldofs(cell)
        kₑ        = zeros(6,6)
        fₑ        = zeros(6)
        kₑ, fₑ    = assemElem(coord[enod[ie][2:end],:],Ψ[cell_dofs],mp,t)
        ke        = zeros(6,6)
        fe        = zeros(6)
        for face in 1:nfaces(cell)
            if (cellid(cell), face) in Γ_robin 
                face_nods = [ Ferrite.facedof_indices(ip)[face][1]; Ferrite.facedof_indices(ip)[face][2] ]
                face_dofs = [ face_nods[1]*2-1; face_nods[1]*2; face_nods[2]*2-1; face_nods[2]*2 ]
                X         = coord[ enod[ie][face_nods.+1] ,: ]
                ke[face_dofs,face_dofs],fe[face_dofs]   = Robin(X,Ψ[cell_dofs[face_dofs]],d[cell_dofs[face_dofs]],λ)
            end
        end
        assemble!(assembler, cell_dofs, kₑ+ke, fₑ+fe)
    end            
end

function assemGlobal!(Fₑₓₜ,dh,t,a,coord,enod,Γt,τ)
    #assembler = start_assemble(Fₑₓₜ)
    #Fₑₓₜ      = zeros(size(coord,1)*2)
    ie        = 0
    kₑ        = zeros(6,6)
    fₑ        = zeros(6)
    for cell in CellIterator(dh)
        fill!(fₑ,0.0)
        ie       += 1
        cell_dofs = celldofs(cell)
        for face in 1:nfaces(cell)
            if (cellid(cell), face) in Γt
                face_nods      = [Ferrite.facedof_indices(ip)[face][1]; Ferrite.facedof_indices(ip)[face][2]]
                face_dofs      = [face_nods[1]*2-1; face_nods[1]*2; face_nods[2]*2-1; face_nods[2]*2]
                X              = coord[enod[ie][face_nods.+1] ,:]
                fₑ[face_dofs]  = tractionLoad(X,τ)
            end    
        end
        #assemble!(assembler, cell_dofs, fₑ)
        Fₑₓₜ[cell_dofs] -= fₑ
    end            
end

function assemGlobal!(K,Fᵢₙₜ,dh,mp,t,a,coord,enod,Γₘ,Γₛ,ϵ)
    assembler = start_assemble(K,Fᵢₙₜ)
    ie        = 0
    kₑ        = zeros(6,6)
    fₑ        = zeros(6)
    for cell in CellIterator(dh)
        fill!(kₑ,0.0)
        fill!(fₑ,0.0)
        ie       += 1
        cell_dofs = celldofs(cell)
        kₑ, fₑ    = assemElem(coord[enod[ie][2:end],:],a[cell_dofs],mp,t)
        # assemble into global
        assemble!(assembler, cell_dofs, kₑ, fₑ)
    end
    
    # Contact
    X_ordered                      = getXfromCoord(coord)
    rc                             = contact_residual(X_ordered,a,ϵ)
    Kc                             = ForwardDiff.jacobian( u -> contact_residual(X_ordered,u,ϵ), a)
    K[contact_dofs, contact_dofs] -= Kc[contact_dofs, contact_dofs]

    Fᵢₙₜ[contact_dofs]            -= rc[contact_dofs]
    #assemble!(assembler, contact_dofs, Kc[contact_dofs, contact_dofs], rc[contact_dofs])
  
end

function volume(dh,coord,enod)
    Ω   = 0.0
    ie  = 0
    for cell in CellIterator(dh)
        ie  += 1
        #println(cellid(cell))
        Ω   += dΩ(coord[enod[ie][2:end],:])
    end
    return Ω
end

function StressExtract(dh,a,mp)
    σ = zeros(dh.ndofs.x)
    for cell in CellIterator(dh)
        # Compute stresses in gauss points - convert to Cauchy
        # Extract GP-stresses to nodes
        # 
        eff    = 
        stress = mises(eff,mp)
    end
end