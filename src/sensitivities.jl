function compliance(F_ext,u)
    C = transpose(F_ext)*u
    return C
end

function drψ(dr_dd,dh,a,fv,λ,d,ΓN)
    assembler = start_assemble(dr_dd)
    ie = 0
    kₑ = zeros(12,12)
    fₑ = zeros(12)
    for cell in CellIterator(dh)
        ie += 1
        cell_dofs= celldofs(cell)
        dfe = zeros(12,12) ## 12x1 eller 12x12 ?? 
        #dfe = RobinIntegral(dfe,cell,ΓN,fv,a[cell_dofs],λ,d[cell_dofs])
        dfe = d_RobinIntegral(dfe,cell,ΓN,fv,a[cell_dofs],λ,d[cell_dofs])
        assemble!(assembler, cell_dofs, dfe)
    end
    return dr_dd
end

function d_RobinIntegral(dfe,cell,ΓN,fv,uₑ,λ,dₑ)
    for face in 1:nfaces(cell)
        if (cellid(cell), face) in ΓN
            Ferrite.reinit!(fv, cell, face)
            for q_point in 1:getnquadpoints(fv)
                #t = 1 * getnormal(fv, q_point)
                dΓ = getdetJdV(fv, q_point)
                for i in 1:12
                    Ni = shape_value(fv, q_point, i)
                    for j in 1:12
                        Nj = shape_value(fv, q_point, j)
                        dfe[i,j]   += Ni ⋅ Nj * ( - λ ) * dΓ
                    end
                end
            end
        end
    end
    return dfe
end

function dcompliance_dx(F_ext)
    return F_ext # Detta är en implicit derivata i kedjeregeln F_ext * du/dx
end

function drᵤ_dx(dr,dh,mp,t,a,coord,enod)
    assembler = start_assemble(dr)
    ie = 0
    drₑ = zeros(12,12)
    for cell in CellIterator(dh)
        ie += 1
        cell_dofs= celldofs(cell)
        drₑ = assem_dr(coord[enod[ie][2:7],:],a[cell_dofs],mp,t)
        assemble!(assembler, cell_dofs, drₑ)
    end 
    return dr
end

function assem_dr(coord,ed,mp,t)
    drₑ = zeros(12,12)
    dre = zeros(12,12)
    for gp ∈ 1 : 3
        dre= dr_GP(coord,ed,gp,mp,t)
        drₑ += dre
    end
    return dre
end

