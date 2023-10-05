# Start gmsh
Gmsh.initialize()

# Open initial file
gmsh.open(filename1*".msh")

# Get nodes
nodeTags, nodeCoords, parametric = gmsh.model.mesh.getNodes()
# mesh_nodes = Dict(zip(nodeTags, nodeCoords))

# # Format displacement field
# Ψ_sorted = getDisplacementsOrdered(dh, Ψ) # reshape_to_nodes(dh, Ψ, :u)

# # Apply displacement field
# parametricCoord = [0.0]
for (nodeTag, displacement) in zip(nodeTags, Ψ_sorted )
    #if nodeTag % 3 != 0
            #=
            # Get the original coordinates of the node
            original_coords = mesh_nodes[nodeTag]

            # Apply the displacement to the original coordinates
            displaced_coords = original_coords .+ displacement

            # Update the node's coordinates in Gmsh
            gmsh.model.mesh.setNode(nodeTag, displaced_coords)
            =#
        gmsh.model.mesh.setNode(nodeTag, [dh.grid.nodes[nodeTag].x; 0.0],[])
    #end
end

# Optimize mesh
#gmsh.model.mesh.optimize()

# Write to file
grid = mktempdir() do dir
    path = joinpath(dir, filename * ".msh")
    gmsh.write(path)
    togrid(path)
end
Gmsh.finalize()
# return grid
