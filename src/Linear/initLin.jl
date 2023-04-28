#load_files()

grid1 = createBoxMesh("box_1",0.0,0.0,1.0,1.0,0.1)

dh = DofHandler(grid1)
add!(dh, :u, 2)
close!(dh)

mp₀    = [1.0 1.0]
t      = 1.0

ip      = Lagrange{2, RefTetrahedron, 1}()
qr      = QuadratureRule{2, RefTetrahedron}(1)
qr_face = QuadratureRule{1, RefTetrahedron}(1)
cv      = CellVectorValues(qr, ip)
fv      = FaceVectorValues(qr_face, ip)

coord, enod = getTopology(dh)


addfaceset!(dh.grid,"Γₜ", x -> x[1] ≈ 1.0)
Γt = getfaceset(dh.grid, "Γₜ")

## Robin faces
#addfaceset!(dh.grid,"Γ_1", x -> x[1] ≈ 1.0)
addfaceset!(dh.grid,"Γ_2", x -> x[2] ≈ 1.0)
#addfaceset!(dh.grid,"Γ3", x -> x[2] ≈ 0.0)
Γ_robin =union(
#    getfaceset(dh.grid, "Γ_1"),
    getfaceset(dh.grid, "Γ_2"),
#    getfaceset(dh.grid, "Γ3"),
)
# Robin nodes
#addnodeset!(dh.grid,"n1", x -> x[1] ≈ 1.0)
addnodeset!(dh.grid,"n2", x -> x[2] ≈ 1.0)
#addnodeset!(dh.grid,"n3", x -> x[2] ≈ 0.0)
n_robin =union(
#    getnodeset(dh.grid, "n1"),
    getnodeset(dh.grid, "n2"),
#    getnodeset(dh.grid, "n3"),
)
coord₀   = deepcopy(coord)
dh0      = deepcopy(dh)





function index_nod_to_grid(dh,coord)
    coord = getCoord(getX(dh),dh);
    X = getX(dh)
    X_nods = reshape_to_nodes(dh, X, :u)[1:2,:] 
    index_register = zeros(Int,length(dh.grid.nodes),2)
    for ii in 1:length(coord[:,1])
        temp2 = X_nods[:,ii]
        for jj in 1:length(coord[:,1])
            temp1 = coord[jj,:]
            if temp1 == temp2
                index_register[ii,:] = [ii,jj]
            end
        end
    end
    return index_register
end

register = index_nod_to_grid(dh,coord)

free_d   = [] # borde vara Vector{Int64}[] men då funkar inte append....

#nodx = Ferrite.getnodeset(dh.grid,"n1")
nody = Ferrite.getnodeset(dh.grid,"n2")

#for inod in nodx
#   append!(free_d,register[inod,2]*2-1)
#end

for jnod in nody
   append!(free_d,register[jnod,2]*2)
end

d         = 0*ones(size(coord,1)*2) 
d[free_d].= 0.1

function postprocess_opt(Ψ, dh, str)
    begin
        vtk_grid(str, dh) do vtkfile
            vtk_point_data(vtkfile, dh, Ψ)
        end
    end
end

function postprocess(a,dh)
        begin
        vtk_grid("hyperelasticity_2", dh) do vtkfile
            vtk_point_data(vtkfile, dh, a)
        end
    end
end