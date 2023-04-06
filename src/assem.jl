# Funktion "assemGlobal"
f1dofs = [1,2,3,4,7,8]
f1  = [1,2,4]

#f1dofs = [1,2,7,8,3,4]
#f1  = [1,4,2]

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
        
        #kₑ,fₑ = RobinIntegral(kₑ,fₑ,cell,ΓN,fv,a[cell_dofs],λ,d[cell_dofs],coord[enod[ie][2:7],:])
        #for face in 1:nfaces(cell)
        #    if (cellid(cell), face) in Γ2  #|| (cellid(cell), face) in Γ1
        #        ke[f1dofs,f1dofs],fe[f1dofs]   = Robin(coord[enod[ie][f1.+1],:],a[cell_dofs[f1dofs]],d[cell_dofs[f1dofs]],λ)
        #        #ke[1:2:end-1,1:2:end-1]       .= 0.0
        #        #fe[1:2:end-1]                 .= 0.0
        #        #println(fe[f1dofs])
        #        println("Element ", ie )
        #        println(ke[2,2], " ", fe[2]," ", fₑ[2])
        #        println(ke[4,4], " ", fe[4]," ", fₑ[4])
        #        println(ke[8,8], " ", fe[8]," ", fₑ[8])
        #    end
        #end
        tn = 1
        ge = zeros(12)
        for face in 1:nfaces(cell)
            if (cellid(cell), face) in Γ2
                Ferrite.reinit!(fv, cell, face)
                for q_point in 1:getnquadpoints(fv)
                    tt = tn * getnormal(fv, q_point)
                    dΓ = getdetJdV(fv, q_point)
                    for i in 1:12
                        δui = shape_value(fv, q_point, i)
                        ge[i] -= (δui ⋅ tt) * dΓ
                    end
                end
            end
        end
        assemble!(assembler, cell_dofs, kₑ+ke, fₑ+fe+ge)
    end            
end

