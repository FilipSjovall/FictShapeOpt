# Funktion "assemGlobal"
function assemGlobal!(K,Fᵢₙₜ,dh,mp,t,a,coord,enod)
    assembler = start_assemble(K,Fᵢₙₜ)
    ie = 0
    kₑ = zeros(12,12)
    fₑ = zeros(12)
    for cell in CellIterator(dh)
        #fill!(kₑ,0.0)
        #fill!(fₑ,0.0)
        ie += 1
        cell_dofs= celldofs(cell)
        kₑ, fₑ = assemElem(coord[enod[ie][2:7],:],a[cell_dofs],mp,t)
        
        assemble!(assembler, cell_dofs, kₑ, fₑ)
    end            
end




function assemGlobal!(K,Fᵢₙₜ,dh,mp,t,a,coord,enod,fv,λ,d,ΓN)
    assembler = start_assemble(K,Fᵢₙₜ)
    ie = 0
    kₑ = zeros(12,12)
    fₑ = zeros(12)
    for cell in CellIterator(dh)
        ie += 1
        cell_dofs= celldofs(cell)
        kₑ, fₑ = assemElem(coord[enod[ie][2:7],:],a[cell_dofs],mp,t)
        ke = zeros(12,12)
        fe = zeros(12)
        ke, fe = RobinIntegral(ke,fe,cell,ΓN,fv,a[cell_dofs],λ,d[cell_dofs])
        assemble!(assembler, cell_dofs, kₑ+ke, fₑ+fe)
    end            
end