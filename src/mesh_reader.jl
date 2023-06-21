#
# Read mesh data from file
#

using  Ferrite
using  FerriteGmsh
using  FerriteMeshParser
using  Gmsh
#file = open(mesh_path,"r")
function getMesh_ASCII(filename)

    mesh_path = "data//Quadratic//"*filename
    n    = 0
    file    = open(mesh_path,"r")
    n       = 0
    nstart  = 0
    nend    = 0
    elstart = 0
    eled    = 0

    coord   = Array{Float64}[]
    enod    = Array{Int64}[]

    indices = [1 6 7 8 9 10 11]

    # Mark where to reset to
    mark(file)
    for line in readlines(file)
        n += 1
        if chomp(line)=="\$NOD"
            nstart  = n
        elseif chomp(line) == "\$ENDNOD"
            nend    = n
        elseif chomp(line) == "\$ELM"
            elstart = n
        elseif chomp(line) == "\$ENDELM"
            eled    = n
        end
    end
    reset(file)


    n = 0
    while ! eof(file)
        n += 1
        line = readline(file)
        if n > nstart + 1 && n <= nend - 1
            numb = parse.(Float64,split(line," ")[2:3])
            push!(coord, numb)
        elseif n > elstart + 1 && n < eled
            println(line)
            numb = parse.(Int64, split(line, " ")[indices])
            push!(enod,numb)
        else

        end
    end

    close(file)

    return mapreduce(permutedims, vcat, coord), enod
end

function getEdof(enod)
    edof = Array{Int64,2}(undef,length(enod),2*(size(enod[1],2)-1))
    for el = 1 : enod[end][1]
        edof[el,1:2:11] = enod[el][2:7]*2 - ones(6)
        edof[el,2:2:12] = enod[el][2:7]*2
    end

    return edof
end

function readAscii(filename)
    coord,enod = getMesh_ASCII(filename)
    edof       = getEdof(enod)
    return coord, enod, edof
end

function modify_msh(filename)
    mesh_path = "data//"*filename
    n    = 0
    file    = open(mesh_path,"w")

    while ! eof(file) && (flag_1 == false )
        line = readline(file)
        if chomp(line)=="\$NOD"
            flag_2 = true
        elseif chomp(line) == "\$ENDNOD"
            flag_1 == true
        end
        if flag_2 == true && flag_1 == false
            print(line[1:3])
        end
    end
end

function setBC(bc_load,dh)
    ## Find bc cell_dofs
    bc_dof = Vector{Int64}()
    bc_val = Vector{Float64}()
    for cell in CellIterator(dh)
        c = getcoordinates(cell)
        for i in 1:3
            x = c[i][1]
            y = c[i][2]
            if y == 0.0
                idx = celldofs(cell)[2*i]
                if idx ∉ bc_dof
                    push!(bc_dof,idx)
                    push!(bc_val,0.0)
                end
            end
            if y == 1.0
                idx = celldofs(cell)[2*i]
                if idx ∉ bc_dof
                push!(bc_dof,idx)
                push!(bc_val,0.0)
                end
            end
            if x == 0.0
                idx = celldofs(cell)[2*i-1]
                if idx ∉ bc_dof
                    push!(bc_dof,idx)
                    push!(bc_val,0.0)
                end
            end
            #if x == 0.0
            #    idx = celldofs(cell)[2*i]
            #    if idx ∉ bc_dof
            #        push!(bc_dof,idx)
            #        push!(bc_val,bc_load)
            #    end
            #end
            if x == 1.0
                idx = celldofs(cell)[2*i-1]
                if idx ∉ bc_dof
                    push!(bc_dof,idx)
                    push!(bc_val,bc_load)
                end
            end
        end
    end
    return bc_dof, bc_val
end

function setup_grid(h)
    # Initialize gmsh
    Gmsh.initialize()
    gmsh.option.set_number("General.Verbosity", 2)
    # Add the points
    o  = gmsh.model.geo.add_point(0.0, 0.0, 0.0, h)
    p1 = gmsh.model.geo.add_point(0.5, 0.0, 0.0, h)
    p2 = gmsh.model.geo.add_point(1.0, 0.0, 0.0, h)
    p3 = gmsh.model.geo.add_point(1.0, 1.0, 0.0, h)
    p4 = gmsh.model.geo.add_point(0.0, 1.0, 0.0, h)
    p5 = gmsh.model.geo.add_point(0.0, 0.5, 0.0, h)
    p6 = gmsh.model.geo.add_point(0.5, 0.5, 0.0, h)

    # Add the lines
    l1 = gmsh.model.geo.add_line(p1, p2)
    l2 = gmsh.model.geo.add_line(p2, p3)
    l3 = gmsh.model.geo.add_line(p3, p4)
    l4 = gmsh.model.geo.add_line(p4, p5)
    l5 = gmsh.model.geo.add_line(p5, p6)
    l6 = gmsh.model.geo.add_line(p6, p1)

    # Create the closed curve loop and the surface
    loop = gmsh.model.geo.add_curve_loop([l1, l2, l3, l4, l5, l6])
    surf = gmsh.model.geo.add_plane_surface([loop])

    # Synchronize the model
    gmsh.model.geo.synchronize()

    # Create the physical domains
    gmsh.model.add_physical_group(1, [l1], -1, "Γ1")
    gmsh.model.add_physical_group(1, [l2], -1, "Γ2")
    gmsh.model.add_physical_group(1, [l3], -1, "Γ3")
    gmsh.model.add_physical_group(1, [l4], -1, "Γ4")
    gmsh.model.add_physical_group(1, [l5], -1, "Γ5")
    gmsh.model.add_physical_group(1, [l6], -1, "Γ6")
    gmsh.model.add_physical_group(2, [surf])

    # Generate a 2D mesh
    gmsh.model.mesh.generate(2)

    # Save the mesh, and read back in as a Ferrite Grid
    grid = mktempdir() do dir
        path = joinpath(dir, "mesh.msh")
        gmsh.write(path)
        togrid(path)
    end

    # Finalize the Gmsh library
    Gmsh.finalize()

    return grid
end

function updateCoords!(dh,Ψ)
    c        = similar(dh.grid.nodes)
    Ψ_sorted = sortNodalDisplacements(dh,Ψ)
    x_new    = Ferrite.Vec{2,Float64}
    for i in 1:length(c)
        x_new = Ferrite.Vec{2}(dh.grid.nodes[i].x + Ψ_sorted[:,i])
        c[i] = Node(x_new)
    end
    copyto!(dh.grid.nodes, c)
end

function sortNodalDisplacements(dh,a)
    return reshape_to_nodes(dh, a, :u)[1:2,:]
end

function getX(dh)
    # Size dh.ndofs.x
    x = Float64[]
    for node in dh.grid.nodes
        append!(x,node.x)
    end
    return x
end

function getXordered(dh)
    # Size dh.ndofs.x
    x = Float64[]
    y = Float64[]
    for node in dh.grid.nodes
        append!(x,node.x[1])
        append!(y,node.x[2])
    end
    return [x;y]
end

function getDisplacementsOrdered(dh,a::AbstractVector{T}) where T
    # This orders a as
    ax = Real[]
    ay = Real[]
    anew = Real[]
    for (idx,node) in enumerate(dh.grid.nodes)
        xdof, ydof = register[idx, :]#[2idx-1 2idx]#
        axn       = a[xdof]
        ayn       = a[ydof]
        append!(ax,axn)
        append!(ay,ayn)
        append!(anew,axn)
        append!(anew,ayn)
    end
    return anew
    #return [ax;ay]
end

function getXorderedDict(xDictionary)
    x = Float64[]
    y = Float64[]
    for (key,value) in xDictionary
        append!(x,value[1])
        append!(y,value[2])
    end
    return [x;y]
end

function getCoord(X::AbstractVector{T},dh) where T
    coord = zeros(length(dh.grid.nodes),2)
    coord[:,1] = X[1:2:end-1]
    coord[:,2] = X[2:2:end]
    return coord
end

function getCoordVector(X::AbstractVector{T}) where T<:Number
    n = length(X)
    coord0 = reshape(X, (Int(n/2), 2))
    return coord0
end

function createBoxMesh(filename,x₀,y₀,Δx,Δy,h)

    # Initialize gmsh
    Gmsh.initialize()
    gmsh.option.set_number("General.Verbosity", 2)


    # Add the points
    p1 = gmsh.model.geo.add_point(x₀,    y₀,    0.0, h)
    p2 = gmsh.model.geo.add_point(x₀+Δx, y₀,    0.0, h)
    p3 = gmsh.model.geo.add_point(x₀+Δx, y₀+Δy, 0.0, h)
    p4 = gmsh.model.geo.add_point(x₀,    y₀+Δy, 0.0, h)

    # Add the lines
    l1 = gmsh.model.geo.add_line(p1, p2)
    l2 = gmsh.model.geo.add_line(p2, p3)
    l3 = gmsh.model.geo.add_line(p3, p4)
    l4 = gmsh.model.geo.add_line(p4, p1)

    # Create the closed curve loop and the surface
    loop = gmsh.model.geo.add_curve_loop([l1, l2, l3, l4])
    surf = gmsh.model.geo.add_plane_surface([loop])

    # Synchronize the model
    gmsh.model.geo.synchronize()

    # Create the physical domains
    gmsh.model.add_physical_group(1, [l1], -1, "Γ")
    gmsh.model.add_physical_group(1, [l2], -1, "Γ2")
    gmsh.model.add_physical_group(1, [l3], -1, "Γ3")
    gmsh.model.add_physical_group(1, [l4], -1, "Γ4")
    gmsh.model.add_physical_group(2, [surf])

    gmsh.model.mesh.generate(2)

    # Save the mesh, and read back in as a Ferrite Grid
    grid = mktempdir() do dir
        path = joinpath(dir, filename*".msh")
        gmsh.write(path)
        togrid(path)
    end

    # Finalize the Gmsh library
    Gmsh.finalize()

    return grid
end


function createBoxMeshRev(filename, x₀, y₀, Δx, Δy, h)

    # Initialize gmsh
    Gmsh.initialize()
    gmsh.option.set_number("General.Verbosity", 2)


    # Add the points
    p1 = gmsh.model.geo.add_point(x₀, y₀, 0.0, h)
    p2 = gmsh.model.geo.add_point(x₀ + Δx, y₀, 0.0, h)
    p3 = gmsh.model.geo.add_point(x₀ + Δx, y₀ + Δy, 0.0, h)
    p4 = gmsh.model.geo.add_point(x₀, y₀ + Δy, 0.0, h)

    # Add the lines
    l1 = gmsh.model.geo.add_line(p1, p4)
    l2 = gmsh.model.geo.add_line(p4, p3)
    l3 = gmsh.model.geo.add_line(p3, p2)
    l4 = gmsh.model.geo.add_line(p2, p1)

    # Create the closed curve loop and the surface
    loop = gmsh.model.geo.add_curve_loop([l1, l2, l3, l4])
    surf = gmsh.model.geo.add_plane_surface([loop])

    # Synchronize the model
    gmsh.model.geo.synchronize()

    # Create the physical domains
    gmsh.model.add_physical_group(1, [l1], -1, "Γ")
    gmsh.model.add_physical_group(1, [l2], -1, "Γ2")
    gmsh.model.add_physical_group(1, [l3], -1, "Γ3")
    gmsh.model.add_physical_group(1, [l4], -1, "Γ4")
    gmsh.model.add_physical_group(2, [surf])

    gmsh.model.mesh.generate(2)

    # Save the mesh, and read back in as a Ferrite Grid
    grid = mktempdir() do dir
        path = joinpath(dir, filename * ".msh")
        gmsh.write(path)
        togrid(path)
    end

    # Finalize the Gmsh library
    Gmsh.finalize()

    return grid
end

function createCircleMesh(filename, x₀, y₀, r, h)

    # Initialize gmsh
    Gmsh.initialize()
    gmsh.option.set_number("General.Verbosity", 2)


    p1 = gmsh.model.geo.add_point(x₀, y₀, 0.0, h)
    p2 = gmsh.model.geo.add_point(x₀ + r, y₀, 0.0, h/2)
    p3 = gmsh.model.geo.add_point(x₀ - r, y₀, 0.0, h/2)
    p4 = gmsh.model.geo.add_point(x₀, y₀ - r, 0.0, h/4 )

    # Add lines
    #    # Start - Center - End
    l1 = gmsh.model.geo.add_circle_arc(p2, p1, p4)
    l2 = gmsh.model.geo.add_circle_arc(p4, p1, p3)
    l3 = gmsh.model.geo.add_line(p3, p1)
    l4 = gmsh.model.geo.add_line(p1, p2)

    # Create the closed curve loop and the surface
    loop = gmsh.model.geo.add_curve_loop([l4, l3, l2, l1])
    surf = gmsh.model.geo.add_plane_surface([loop])

    # Synchronize the model
    gmsh.model.geo.synchronize()

    # Create the physical domains
    gmsh.model.add_physical_group(1, [l1], -1, "Γ")
    gmsh.model.add_physical_group(1, [l2], -1, "Γ2")
    gmsh.model.add_physical_group(1, [l3], -1, "Γ3")
    gmsh.model.add_physical_group(1, [l4], -1, "Γ3")
    gmsh.model.add_physical_group(2, [surf])

    gmsh.model.mesh.generate(2)

    # Save the mesh, and read back in as a Ferrite Grid
    grid = mktempdir() do dir
        path = joinpath(dir, filename * ".msh")
        gmsh.write(path)
        togrid(path)
    end

    # Finalize the Gmsh library
    Gmsh.finalize()

    return grid
end

function merge_grids(grid1::Grid{dim,CellType}, grid2::Grid{dim,CellType}; tol=1e-6) where {N, dim, CellType <: Cell{<:Any, N}}
    cells′ = copy(grid1.cells)
    nodes′ = copy(grid1.nodes)
    nodemap = Dict{Int,Int}()
    next = getnnodes(grid1) + 1
    for (i2, n2) in enumerate(grid2.nodes)
        found = false
        for (i1, n1) in enumerate(grid1.nodes)
            if norm(n1.x - n2.x) < tol
                nodemap[i2] = i1
                found = true
                break
            end
        end
        if !found
            push!(nodes′, n2)
            nodemap[i2] = next
            next += 1
        end
    end
    for c in grid2.cells
        t = ntuple(N) do i
            return nodemap[c.nodes[i]]
        end
        cell′ = CellType(t)
        push!(cells′, cell′)
    end
    return Grid(cells′, nodes′)
end

function getTopology(dh)
    coord   = zeros( length(dh.grid.nodes), 2)
    enod    = Array{Int64}[]
    for cell in CellIterator(dh)
        push!(enod,[ [cellid(cell)]; cell.nodes])
    end
    nod_nr = 0
    for node in dh.grid.nodes
        nod_nr += 1
        coord[nod_nr,:] = node.x
    end
    return coord, enod
end

function setBCLin(bc_load,dh)
    ## Find bc cell_dofs
    bc_dof = Vector{Int64}()
    bc_val = Vector{Float64}()
    for cell in CellIterator(dh)
        c = getcoordinates(cell)
        for i in 1:3
            x = c[i][1]
            y = c[i][2]
            if x == 0.0
                idx = celldofs(cell)[2*i-1]
                if idx ∉ bc_dof
                    push!(bc_dof,idx)
                    push!(bc_val,bc_load)
                end
                if y == 0.0
                    #idx = celldofs(cell)[2*i-1]
                    #if idx ∉ bc_dof
                    #    push!(bc_dof,idx)
                    #    push!(bc_val,bc_load)
                    #end
                    idx = celldofs(cell)[2*i]
                    if idx ∉ bc_dof
                        push!(bc_dof,idx)
                        push!(bc_val,bc_load)
                    end
                end
            end
        end
    end
    return bc_dof, bc_val
end

function setBC(bc_load,dh,nodes)
    ## Find bc cell_dofs
    bc_dof   = Vector{Int64}()
    bc_val   = Vector{Float64}()
    nodeDofs = getNodeDofs(dh)
    for node in nodes
        idx1 = nodeDofs[node,1]
        idx2 = nodeDofs[node,2]
        #append!(bc_dof,idx1)
        append!(bc_dof,idx2)
        append!(bc_val,bc_load)
    end
    return bc_dof, bc_val
end

function setBCXY(bc_load, dh, nodes)
    ## Find bc cell_dofs
    bc_dof = Vector{Int64}()
    bc_val = Vector{Float64}()
    nodeDofs = getNodeDofs(dh)
    for node in nodes
        xdof = nodeDofs[node, 1]
        ydof = nodeDofs[node, 2]
        append!(bc_dof,xdof)
        append!(bc_dof, ydof)
        append!(bc_val, 0.0)
        append!(bc_val, bc_load)
    end
    return bc_dof, bc_val
end

function getNodeDofs(dh)
    node_dofs = Matrix{Int64}(undef,length(dh.grid.nodes),2)
    for cell in CellIterator(dh)
        element_dofs = celldofs(cell)
        for (i,node) in enumerate(cell.nodes)
            node_dofs[node,:] = [element_dofs[2i-1] element_dofs[2i]]
        end
    end
    return node_dofs
end

function getXfromCoord(coord::AbstractVector{T}) where T
    X = Real[]
    for row ∈ 1:size(coord,1) # eachindex(coord)
        append!(X,coord[row,1])
        append!(X,coord[row,2])
    end
    return X
end

function getCoordfromX(X::AbstractVector{T}) where T
    #coord = zeros(length(dh.grid.nodes),2)
    coordz = Array{Real,2}(undef,length(dh.grid.nodes),2)
    for incr in 2:2:size(coordz,1)*2
        coordz[Int(incr/2),:] = [X[incr-1] X[incr]]
    end
    return coordz
end

function getContactDofs(nₛ, nₘ)
    contact_dofs = Int64[]
    register = getNodeDofs(dh)
    for node_nbr in nₛ
        for dof in 1:2
            append!(contact_dofs, register[node_nbr, dof])
        end
    end
    for node_nbr in nₘ
        for dof in 1:2
            append!(contact_dofs, register[node_nbr, dof])
        end
    end
    return contact_dofs
end

function getContactNods(nₛ, nₘ)
    contact_nods = Int64[]
    for node_nbr in nₛ
            append!(contact_nods, node_nbr)
    end
    for node_nbr in nₘ
            append!(contact_nods, node_nbr)
    end
    return contact_nods
end

function getXfromCoord(coord)
    X = Real[]
    for row ∈ 1:size(coord, 1) # eachindex(coord)
        append!(X, coord[row, 1])
        append!(X, coord[row, 2])
    end
    return X
end

function getXinDofOrder(dh,X,coord)
    X_ordered = similar(X)
    for node ∈ eachindex(coord[:, 1])
        dofs = register[node,:]
        X_ordered[dofs] = coord[node,:]
    end
    return X_ordered
end

function getX_from_Dof_To_Node_order(dh,X::AbstractVector{T}) where T
    X_ordered = similar(X) # Real[]
    for (node, dofs) ∈ enumerate( eachrow(register) )
        X_ordered[2node - 1] = X[dofs[1]]
        X_ordered[2node]     = X[dofs[2]]
    end
    return X_ordered
end

function index_nod_to_grid(dh, coord)
    coord = getCoord(getX(dh), dh)
    X = getX(dh)
    X_nods = reshape_to_nodes(dh, X, :u)[1:2, :]
    index_register = zeros(Int, length(dh.grid.nodes), 2)
    for ii in 1:length(coord[:, 1])
        temp2 = X_nods[:, ii]
        for jj in 1:length(coord[:, 1])
            temp1 = coord[jj, :]
            if temp1 == temp2
                index_register[ii, :] = [ii, jj]
            end
        end
    end
    return index_register
end
