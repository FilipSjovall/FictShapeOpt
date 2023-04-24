f1dofs = [1,2,3,4] # ????
f1  = [1,2] # ???

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
    ie = 0
    kₑ = zeros(6,6)
    fₑ = zeros(6)
    for cell in CellIterator(dh)
        fill!(kₑ,0.0)
        fill!(fₑ,0.0)
        ie += 1
        cell_dofs= celldofs(cell)
        kₑ, fₑ = assemElem(coord[enod[ie][2:end],:],a[cell_dofs],mp,t)
        for face in 1:nfaces(cell)
            if (cellid(cell), face) in Γt
                face_nods = [ Ferrite.facedof_indices(ip)[face][1]; Ferrite.facedof_indices(ip)[face][2] ]
                face_dofs = [ face_nods[1]*2-1; face_nods[1]*2; face_nods[2]*2-1; face_nods[2]*2 ]
                X = coord[ enod[ie][face_nods.+1] ,: ]
                fₑ[face_dofs] += tractionLoad( X, τ )
            end    
        end
        assemble!(assembler, cell_dofs, kₑ, fₑ)
    end            
end

function assemGlobal!(K,Fᵢₙₜ,dh,mp,t,a,coord,enod,fv,λ,d,ΓN)
    assembler = start_assemble(K,Fᵢₙₜ)
    ie = 0
    for cell in CellIterator(dh)
        ie += 1
        cell_dofs= celldofs(cell)
        kₑ = zeros(6,6)
        fₑ = zeros(6)
        kₑ, fₑ    = assemElem(coord[enod[ie][2:end],:],a[cell_dofs],mp,t)
        ke = zeros(6,6)
        fe = zeros(6)
        for face in 1:nfaces(cell)
            if (cellid(cell), face) in Γt  #|| (cellid(cell), face) in Γ1
                ke[f1dofs,f1dofs],fe[f1dofs]   = Robin(coord[enod[ie][f1.+1],:],a[cell_dofs[f1dofs]],d[cell_dofs[f1dofs]],λ)
            end
        end
        assemble!(assembler, cell_dofs, kₑ+ke, fₑ+fe)
    end            
end

function volume(dh,coord)
    Ω   = 0.0
    ie  = 0
    for cell in CellIterator(dh)
        ie  += 1
        #println(cellid(cell))
        Ω   += dΩ(coord[enod[ie][2:end],:])
    end
    return Ω
end

