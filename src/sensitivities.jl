# Funktion "assemGlobal"
f1dofs = [1,2,3,4]
f1  = [1,2]

function compliance(F_ext,u)
    C = transpose(F_ext)*u
    return C
end

function drψ(dr_dd,dh,a,fv,λ,d,ΓN)
    assembler = start_assemble(dr_dd)
    ie = 0
    for cell in CellIterator(dh)
        @show ie += 1
        cell_dofs= celldofs(cell)
        dfe = zeros(12,12) ## 12x1 eller 12x12 ??
        for face in 1:nfaces(cell)
            if (cellid(cell), face) in ΓN
                #@show Ferrite.facedof_indices(ip)[face]
                #face_nods = [ Ferrite.facedof_indices(ip)[face][1]; Ferrite.facedof_indices(ip)[face][2] ]
                #face_dofs = [ face_nods[1]*2-1; face_nods[1]*2; face_nods[2]*2-1; face_nods[2]*2 ]
                #X         = coord[ enod[ie][face_nods.+1] ,: ]
                #dfe[face_dofs,face_dofs],_ = Robin(X,a[cell_dofs[face_dofs]],d[cell_dofs[face_dofs]],λ)
                dfe[f1dofs, f1dofs], _ = Robin(coord[enod[ie][f1.+1], :], a[cell_dofs[f1dofs]], d[cell_dofs[f1dofs]], λ)
            end
        end
        assemble!(assembler, cell_dofs, -dfe)
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
                if (cellid(cell), face) in Γ1
                    for i in 1:2:11
                        Ni  = shape_value(fv, q_point, i)
                        for j in 1:2:11
                            Nj = shape_value(fv, q_point, j)
                            dfe[i,j]   += Ni ⋅ Nj * ( - λ ) * dΓ
                        end
                    end
                elseif (cellid(cell), face) in Γ2
                    for i in 2:2:12
                        Ni  = shape_value(fv, q_point, i)
                        for j in 2:2:12
                            Nj = shape_value(fv, q_point, j)
                            dfe[i,j]   += Ni ⋅ Nj * ( - λ ) * dΓ
                        end
                    end
                end
                #for i in 1:12
                #    Ni = shape_value(fv, q_point, i)
                #    for j in 1:12
                #        Nj = shape_value(fv, q_point, j)
                #        dfe[i,j]   += Ni ⋅ Nj * ( - λ ) * dΓ
                #    end
                #end
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
        drₑ = assem_dr(coord[enod[ie][2:end],:],a[cell_dofs],mp,t)
        assemble!(assembler, cell_dofs, drₑ)
    end
    return dr
end
