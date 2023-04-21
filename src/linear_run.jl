using Gmsh
using FerriteGmsh
using FerriteMeshParser
include("..//src//mesh_reader.jl")
# Box 1
grid1 = createBoxMesh("box_1",0.0,0.0,1.0,1.0,0.1)

grid2 = createBoxMesh("box_2",0.0,1.0,1.0,1.0,0.05)

using Tensors, Ferrite#JuAFEM



#new_grid = merge_grids(grid1::G, grid2::G, grids::G...) where G <: Grid =
#    merge_grids(merge_grids(grid1, grid2), grids...)

new_grid = merge_grids(grid1, grid2; tol=0.01)
dh = DofHandler(new_grid)
add!(dh, :u, 2)
close!(dh)