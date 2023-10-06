function remesh(filename)
    # Start gmsh
    Gmsh.initialize()

    # Open initial file
    gmsh.open(filename1 * ".msh")

    # Get nodes
    nodeTags, nodeCoords, parametric = gmsh.model.mesh.getNodes();

    # # Format displacement field
    Ψ_sorted = getDisplacementsOrdered(dh, Ψ) # reshape_to_nodes(dh, Ψ, :u)

    # # Apply displacement field
    for (i,nodeTag) in enumerate(nodeTags)
        gmsh.model.mesh.setNode(nodeTag, [dh.grid.nodes[i].x; 0.0], [])
    end

    # Specify optimization parameters (adjust as needed)
    gmsh.model.mesh.reclassifyNodes()
    gmsh.model.mesh.createGeometry()
    gmsh.model.mesh.generate(2)
    algorithm = "Netgen" # "Relocate2D" # "HighOrderFastCurving" # "Netgen" #

    # Run the optimization
    gmsh.model.mesh.optimize(algorithm)
    # Write to file
    grid = mktempdir() do dir
        path = joinpath( filename1 * "optimized" * ".msh")
        gmsh.write(path)
        togrid(path)
    end
    Gmsh.finalize()
    return grid
end

testgrid = remesh(filename1)
# return grid
