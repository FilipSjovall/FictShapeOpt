# Funktion "assemGlobal"
f1dofs = [1,2,3,4,7,8]
f1  = [1,2,4]

#f1 = [1,6,4]
#f1dofs = [1,2,11,12,7,8]

#f1dofs = [1,2,7,8,3,4]
#f1  = [1,4,2]

function assemGlobal!(K,Fᵢₙₜ,dh,mp,t,a,coord,enod)
    assembler = start_assemble(K,Fᵢₙₜ)
    ie = 0
    kₑ = zeros(12,12)
    fₑ = zeros(12)
    for cell in CellIterator(dh)
        fill!(kₑ,0.0)
        fill!(fₑ,0.0)
        ie += 1
        cell_dofs= celldofs(cell)
        kₑ, fₑ = assemElem(coord[enod[ie][2:7],:],a[cell_dofs],mp,t)
        assemble!(assembler, cell_dofs, kₑ, fₑ)
    end            
end

function assemGlobal!(K,Fᵢₙₜ,dh,mp,t,a,coord,enod,fv,λ,d,ΓN)
    assembler = start_assemble(K,Fᵢₙₜ)
    ie = 0
    for cell in CellIterator(dh)
        ie += 1
        cell_dofs= celldofs(cell)
        kₑ = zeros(12,12)
        fₑ = zeros(12)
        kₑ, fₑ    = assemElem(coord[enod[ie][2:7],:],a[cell_dofs],mp,t)
        ke = zeros(12,12)
        fe = zeros(12)
        for face in 1:nfaces(cell)
            if (cellid(cell), face) in Γt  #|| (cellid(cell), face) in Γ1
                ke[f1dofs,f1dofs],fe[f1dofs]   = Robin(coord[enod[ie][f1.+1],:],a[cell_dofs[f1dofs]],d[cell_dofs[f1dofs]],λ)
            end
        end
        assemble!(assembler, cell_dofs, kₑ+ke, fₑ+fe)
    end            
end

