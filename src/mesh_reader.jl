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
    c       = similar(dh.grid.nodes)
    Ψ_sorted = sortNodalDisplacements(dh,Ψ)
    x_new = Vec{2,Float64}
    for i in 1:length(c)
        x_new = Vec{2}(dh.grid.nodes[i].x + Ψ_sorted[:,i])
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

function getCoord(X,dh)
    coord = zeros(length(dh.grid.nodes),2)
    coord[:,1] = X[1:2:end-1]
    coord[:,2] = X[2:2:end]
    return coord
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
    #p4 = gmsh.model.geo.add_point(0.0, 0.5, 0.0, h)

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

function merge_grids(grid1::Grid{dim,CellType}, grid2::Grid{dim,CellType}; tol=0.01) where {N, dim, CellType <: Cell{<:Any, N}}
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
        #println(cellid(cell))
        #println(cell.nodes)  
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