
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

