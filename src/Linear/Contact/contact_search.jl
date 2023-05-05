# ------------------- #
## Small test        ##
# ------------------- #
using Ferrite, FerriteGmsh, FerriteMeshParser
using Mortar2D
#using FerriteViz, WGLMakie
include("..//..//mesh_reader.jl")
include("..//initLin.jl")
#include("Mortar2D//Mortar2D.jl")


# Create two grids and a common dofhandler for them
grid1    = createBoxMesh("box_1",0.0,0.0,1.0,1.0,0.3)
grid2    = createBoxMesh("box_2",0.0,1.01,1.0,1.0,0.4)

grid_tot = merge_grids(grid1, grid2; tol=0.01)

dh = DofHandler(grid_tot)
add!(dh, :u, 2)
close!(dh)
coord, enod = getTopology(dh)

# -------------------------------- #
# Polynomial and integration order #
# -------------------------------- #
ip      = Lagrange{2, RefTetrahedron, 1}()
qr      = QuadratureRule{2, RefTetrahedron}(1)
qr_face = QuadratureRule{1, RefTetrahedron}(1)
cv      = CellVectorValues(qr, ip)
fv      = FaceVectorValues(qr_face, ip)

# ------------------ #
# Create master sets #
# ------------------ #
addfaceset!(dh.grid,"Γ_master", x -> x[2] ≈ 1.0)
Γm = getfaceset(dh.grid, "Γ_master")

addnodeset!(dh.grid,"nₘ", x -> x[2] ≈ 1.0)
nₘ = getnodeset(dh.grid,"nₘ")

# ----------------- #
# Create slave sets #
# ----------------- #
addfaceset!(dh.grid,"Γ_s", x -> x[2] ≈ 1.01)
Γs = getfaceset(dh.grid, "Γ_s")

addnodeset!(dh.grid,"nₛ", x -> x[2] ≈ 1.01)
nₛ = getnodeset(dh.grid,"nₛ")

# ------------------- #
# Create Mortar Dicts #
# ------------------- #
i                  = 0 
elements           = Dict{Int64, Vector{Int64}}()
coords             = Dict{Int64, Vector{Float64}}()
slave_elements     = Dict{Int64, Vector{Int64}}()
element_types      = Dict{Int64, Symbol}()
slave_element_ids  = Vector{Int64}()
master_element_ids = Vector{Int64}()
for face in  Γs 
   i        += 1
   face_el   = face[1]
   face_nods = Ferrite.faces(dh.grid.cells[face_el])[1]
   push!(elements,face_el=>[face_nods[1],face_nods[2]])
   push!(coords,face_nods[1]=>dh.grid.nodes[face_nods[1]].x)
   push!(coords,face_nods[2]=>dh.grid.nodes[face_nods[2]].x)
   #push!(coords,i=>dh.grid.nodes[face_nods[1]].x)
   #push!(coords,(i+1)=>dh.grid.nodes[face_nods[2]].x)
   push!(slave_element_ids,face_el)
   push!(slave_elements,face_el=>[face_nods[1],face_nods[2]])
   push!(element_types, face_el => :Seg2)
   #println("face element ", face_el, " with nodes ", face_nods, " that have the coordinates ", dh.grid.nodes[face_nods[1]].x, " ", dh.grid.nodes[face_nods[2]].x)
   #println("nod nbrs: ", i , " ", i+1)
   #println("facenods ", face_nods )
end

for face in  Γm
   i        += 1
   face_el   = face[1]
   face_nods = Ferrite.faces(dh.grid.cells[face_el])[1]
   push!(elements,face_el=>[face_nods[1],face_nods[2]])
   push!(coords,face_nods[1]=>dh.grid.nodes[face_nods[1]].x)
   push!(coords,face_nods[2]=>dh.grid.nodes[face_nods[2]].x)
   #push!(coords,i=>dh.grid.nodes[face_nods[1]].x)
   #push!(coords,(i+1)=>dh.grid.nodes[face_nods[2]].x)
   push!(master_element_ids,face_el)
   push!(element_types, face_el => :Seg2)
end


# ------------------- #
# Create normal field #
# ------------------- #
normals = Mortar2D.calculate_normals(slave_elements, element_types, coords)


Mortar2D.calculate_segments(slave_element_ids, master_element_ids, elements,
                            element_types, coords, normals)

s, m, D, M = Mortar2D.calculate_mortar_assembly(elements, element_types, coords,
                                       slave_element_ids, master_element_ids)

# --------------------------- #
# Help functions - conversion #
# --------------------------- #
# Stämmer inte...
Xm            = zeros(length(master_element_ids)*2)
Xm[1:2:end-1] = coord[master_element_ids,1]
Xm[2:2:end]   = coord[master_element_ids,2]

Xs            = zeros(length(slave_element_ids)*2)
Xs[1:2:end-1] = coord[slave_element_ids,1]
Xs[2:2:end]   = coord[slave_element_ids,2]

g = zeros(length(s),2)
for (j,A) in (enumerate(s))
   slave  = [0;0]
   for B in s
      slave  += D[A,B]*coords[B]  
   end
   master = [0;0]
   for C in m
      master += M[A,C]*coords[C]
   end
   #println(A, " " , j)
   #println(g[j][1])
   g[j,:] = slave - master

end

g_vec = zeros(length(s))

# g_vec måste skalas med 1/n_AD
for (j,A) in (enumerate(s))
   g_vec[j] = g[j,:]' * normals[A]
end

function penalty(g,penalty)
   λ = penalty * max(g,0)
end


using ForwardDiff

function contact_residual()
   rc = zeros(length(s)*2)
   for (j,A) in enumerate(s)
      λ = penalty(g_vec[j],1)*normal[A]
      λ*M[A,:]
      λ*D[A,:]
   end
end
# coord[master_element_ids,:]
# reshape(coord[master_element_ids,:],(8,1))