function drψ(dr_dd,dh,a,λ,d,Γ_robin)
    assembler = start_assemble(dr_dd)
    ie = 0
    for cell in CellIterator(dh)
        ie += 1
        cell_dofs= celldofs(cell)
        dfe = zeros(6,6) ## 12x1 eller 12x12 ?? 
        for face in 1:nfaces(cell)
            if (cellid(cell), face) in Γ_robin
                # Fixa
                face_nods = [ Ferrite.facedof_indices(ip)[face][1]; Ferrite.facedof_indices(ip)[face][2] ]
                face_dofs = [ face_nods[1]*2-1; face_nods[1]*2; face_nods[2]*2-1; face_nods[2]*2 ]
                X         = coord[ enod[ie][face_nods.+1] ,: ]
                dfe[face_dofs,face_dofs],_ = Robin(X,a[cell_dofs[face_dofs]],d[cell_dofs[face_dofs]],λ)
            end
        end
        assemble!(assembler, cell_dofs, -dfe)
    end
    return dr_dd
end

function drᵤ_dx(dr,dh,mp,t,a,coord,enod)
    assembler = start_assemble(dr)
    ie = 0
    drₑ = zeros(6,6)
    for cell in CellIterator(dh)
        ie += 1
        cell_dofs= celldofs(cell)
        drₑ = assem_dr(coord[enod[ie][2:end],:],a[cell_dofs],mp,t)
        assemble!(assembler, cell_dofs, drₑ)
    end 
    return dr
end



