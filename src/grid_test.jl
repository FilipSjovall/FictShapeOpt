
using FEMSparse
using SparseArrays

include("mesh_reader.jl")
filename = "mesh_fine_fine.txt"
coord, enod, edof = readAscii(filename);

ndof = size(coord,1)*2

sparse_pattern = zeros(ndof,ndof)

sparse_pattern[edof,edof] .= 1

K = sparse(sparse_pattern)

assembler = FEMSparse.start_assemble(K)

ke = zeros(ndof, ndof)
ge = zeros(ndof)

using FerriteMeshParser
using Ferrite
# För rätt format kan z-koordinat tas bort i notepad++ med specialsök: ", 0\n" - replace
grid = get_ferrite_grid("data/mesh_fine_fine.inp")
dh = DofHandler(grid)
add!(dh, :u, 2)
close!(dh)
K = create_sparsity_pattern(dh)


ie = 0
for cell in CellIterator(dh)
    ie +=1
    
    println(coord)
end

coord  = getcoordinates(grid,1)


function setBC(bc_load)
    ## Find bc cell_dofs
    bc_dof = Vector{Int64}()
    bc_val = Vector{Float64}()
    for cell in CellIterator(dh)
        c = getcoordinates(cell)
        for i in 1:6
            x = c[i][1]
            y = c[i][2]
            if y == 0.0
                idx = celldofs(cell)[2*i] 
                push!(bc_dof,idx)
                push!(bc_val,0.0)
            elseif x == 0.0
                idx = celldofs(cell)[2*i-1] 
                push!(bc_dof,idx)
                push!(bc_val,0.0)
            elseif x == 0.5
                idx = celldofs(cell)[2*i-1] 
                push!(bc_dof,idx)
                push!(bc_val,bc_load)
            end
        end
    end
    return bc_dof, bc_val
end