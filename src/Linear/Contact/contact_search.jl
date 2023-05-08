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

gap = zeros(length(s))

# gap måste skalas med 1/n_AD
for (j,A) in (enumerate(s))
   gap[j] = g[j,:]' * normals[A]
end

function penalty(g,penalty)
   #λ = penalty * max(g,0)
   if g < 0
      λ = penalty * g
   else
      λ = 0
   end
   return λ
end

using ForwardDiff
X   = getX(dh)
gap = gap_function(X)

ForwardDiff.jacobian(x -> gap_function(x),real(X))

ForwardDiff.gradient(x -> gap_function(x), real(X) )

ϵ = 10.0
using BenchmarkTools
@benchmark ForwardDiff.derivative(x -> penalty(x,ϵ),-0.1)

function getCoordVector(X::AbstractVector{T},dh) where T<:Number
   n = length(X)
   coord0 = reshape(X, (Int(n/2), 2))
   return coord0
end

J_getCoord = ForwardDiff.jacobian(x -> getCoordVector(x, dh), X)



normals      = Mortar2D.calculate_normals(elements, element_types, coords)
segmentation = Mortar2D.calculate_segments(slave_element_ids, master_element_ids,
                                      elements, element_types, coords, normals)

X = getXordered(dh)
elements,element_types, slave_elements, slave_element_ids, master_element_ids, coords = create_contact_list(dh,Γs,Γm, coord)
slave_dofs, master_dofs, D, M                 = Mortar2D.calculate_mortar_assembly(elements, element_types, coords, slave_element_ids, master_element_ids,segmentation)

using Zygote;

J_gap = ForwardDiff.jacobian( x -> gap_function(x,segmentation), X)

gap_function(X, segmentation)

function gap_function(X::AbstractVector{T},segmentation) where T
   X_float = real.(X) # Convert input to Vector{Float64}
   coord  = getCoordVector(X_float,dh)
   elements,element_types, slave_elements, slave_element_ids, master_element_ids, coords = create_contact_list(dh,Γs,Γm, coord)
   slave_dofs, master_dofs, D, M                 = Mortar2D.calculate_mortar_assembly(elements, element_types, coords, slave_element_ids, master_element_ids,segmentation)

   g0      = zeros(eltype(X_float),length(slave_dofs),2)
   for (j,A) in (enumerate(slave_dofs))
      slave  = [0;0]
      for B in slave_dofs
         slave  += D[A,B]*coords[B]  
      end
      master = [0;0]
      for C in master_dofs
         master += M[A,C]*coords[C]
      end
      g0[j,:] = slave - master
      #setindex!(g0,slave-master,j,[1,2])
   end
   return reshape(g0,8,1)
end

function finite_diff_jacobian(f, x, h)
   n = length(x)
   J = zeros(eltype(x), length(f(x)), n)
   for i = 1:n
       x_plus = copy(x)
       x_plus[i] += h
       x_minus = copy(x)
       x_minus[i] -= h
       J[:,i] = (f(x_plus) - f(x_minus)) / (2h)
   end
   return J
end

# Test the AD-computed Jacobian against finite differences
X = rand(10)
segmentation = Dict(1 => [1, 2, 3], 2 => [4, 5, 6], 3 => [7, 8, 9, 10])
J_ad = ForwardDiff.jacobian(x -> gap_function(x, segmentation), X)
J_fd = finite_diff_jacobian(x -> gap_function(x, segmentation), X, 1e-6)
fel = norm(J_ad - J_fd) ./ norm(J_fd)
println("Error: $error")

for i in eachindex(J_ad)
   if J_ad[i] != 0 && J_fd[i] != 0
       error = abs(J_ad[i] - J_fd[i])
       quotient = J_ad[i]/J_fd[i]
       num1     = J_ad[i] 
       num2     = J_fd[i]
       #println("Element $i error: $error ", "quotient: $quotient ")
       println("Element $i quotient: $quotient with values $num1 $num2 ")
   end
end


function create_contact_list(dh,Γs,Γm, coord_dual)
   i                  = 0 
   elements           = Dict{Int64, Vector{Int64}}()
   coords             = Dict{Int64, Vector{Real}}()
   slave_elements     = Dict{Int64, Vector{Int64}}()
   element_types      = Dict{Int64, Symbol}()
   slave_element_ids  = Vector{Int64}()
   master_element_ids = Vector{Int64}()
   for face in  Γs 
      i        += 1
      face_el   = face[1]
      face_nods = Ferrite.faces(dh.grid.cells[face_el])[1]
      push!(elements,face_el=>[face_nods[1],face_nods[2]])
      push!(coords,face_nods[1]=>coord_dual[face_nods[1],:])
      push!(coords,face_nods[2]=>coord_dual[face_nods[2],:])
      push!(slave_element_ids,face_el)
      push!(slave_elements,face_el=>[face_nods[1],face_nods[2]])
      push!(element_types, face_el => :Seg2)
   end
   
   for face in  Γm
      i        += 1
      face_el   = face[1]
      face_nods = Ferrite.faces(dh.grid.cells[face_el])[1]
      push!(elements,face_el=>[face_nods[1],face_nods[2]])
      push!(coords,face_nods[1]=>coord_dual[face_nods[1],:])
      push!(coords,face_nods[2]=>coord_dual[face_nods[2],:])
      push!(master_element_ids,face_el)
      push!(element_types, face_el => :Seg2)
   end
   return elements,element_types, slave_elements, slave_element_ids, master_element_ids, coords
end



