
using Mortar2D, ForwardDiff
using Ferrite, FerriteGmsh, FerriteMeshParser
using LinearSolve, SparseArrays # LinearSolvePardiso
using IterativeSolvers, IncompleteLU    # AlgebraicMultigrid

include("..//mesh_reader.jl")
include("initLin.jl") # initieras massa skit
include("Contact//contact_help.jl")
include("assemLin.jl")
include("assemElemLin.jl")
include("..//material.jl")
include("..//fem.jl")


# Create two grids 
grid1 = createBoxMesh("box_1", 0.0, 0.0, 1.0, 1.0, 0.1)
grid2 = createBoxMesh("box_2", 0.33, 0.99, 0.33, 0.5, 0.05)

# Merge into one grid
grid_tot = merge_grids(grid1, grid2; tol=0.01)

# Create dofhandler with displacement field u
dh = DofHandler(grid_tot)
add!(dh, :u, 2)
close!(dh)

# Extract CALFEM-style matrices
coord, enod = getTopology(dh)

register = index_nod_to_grid(dh, coord)

# --------------------------------------------------------------------------- #
# These are useful if Shape functions/gradients defined by Ferrite are needed #
# --------------------------------------------------------------------------- #
ip = Lagrange{2,RefTetrahedron,1}();
qr = QuadratureRule{2,RefTetrahedron}(1);
qr_face = QuadratureRule{1,RefTetrahedron}(1);
cv = CellVectorValues(qr, ip);
fv = FaceVectorValues(qr_face, ip);

# ------------------ #
# Create master sets #
# ------------------ #
addfaceset!(dh.grid, "Γ_master", x -> x[2] ≈ 1.0)
Γm = getfaceset(dh.grid, "Γ_master")

addnodeset!(dh.grid, "nₘ", x -> x[2] ≈ 1.0)
nₘ = getnodeset(dh.grid, "nₘ")

# ----------------- #
# Create slave sets #
# ----------------- #
addfaceset!(dh.grid, "Γ_s", x -> x[2] ≈ 0.99)
Γs = getfaceset(dh.grid, "Γ_s")

addnodeset!(dh.grid, "nₛ", x -> x[2] ≈ 0.99)
nₛ = getnodeset(dh.grid, "nₛ")

contact_dofs = getContactDofs(nₛ, nₘ)

# Define top nodeset for displacement controlled loading
addnodeset!(dh.grid, "Γ_top", x -> x[2] ≈ 1.49)
Γ_top = getnodeset(dh.grid, "Γ_top")

# Define bottom nodeset subject to  u(X) = 0 ∀ X ∈ Γ_bot
addnodeset!(dh.grid, "Γ_bot", x -> x[2] ≈ 0.0)
Γ_bot = getnodeset(dh.grid, "Γ_bot")




register = getNodeDofs(dh)
X = getX(dh)

coord = getCoordfromX(X)

# Solve the nonlinear equillibrium problem .. 
function solver(dh, coord)

    # ---------- #
    # Set params # // Kanske som input till solver???
    # ---------- # // definiera mp här? och kanske ε ? iofs snyggare utanför!  
    t = 1.0

    # Penalty parameter
    ε = 100.0

    # Define material parameters
    mp = [210 0.3] # [E ν]

    # ------------- #
    # Init-stuff    #
    # ------------- #
    imax = 25
    TOL = 1e-8
    residual = 0.0
    iter = 1
    # ------------- #.0
    # ------------- #
    K = create_sparsity_pattern(dh)

    # ------ #
    #  Init  #
    # ------ #
    global Fᵢₙₜ = zeros(dh.ndofs.x)
    global Fₑₓₜ = zeros(dh.ndofs.x)
    global a    = zeros(dh.ndofs.x)
    global Δa   = zeros(dh.ndofs.x)
    global res  = zeros(dh.ndofs.x)

    # ---------- #
    # Set BCS    # 
    # ---------- # 
    # Set bcs - should be moved outside this function
    bcdof_top, bcval_top = setBCXY(-0.01, dh, Γ_top)
    bcdof_bot, bcval_bot = setBCXY(0.0, dh, Γ_bot)
    bcdof = [bcdof_top; bcdof_bot]
    bcval = [bcval_top; bcval_bot]

    ϵᵢⱼₖ = sortperm(bcdof)
    bcdof = bcdof[ϵᵢⱼₖ]
    bcval = bcval[ϵᵢⱼₖ]

    # - For Linear solver..
    pdofs = bcdof
    fdofs = setdiff(1:dh.ndofs.x, pdofs)

    bcval₀ = bcval

    for loadstep ∈ 1:35
        res = res .* 0
        bcval = bcval₀
        residual = 0 * residual
        iter = 0
        fill!(Δa, 0.0)
        println("Starting equillibrium iteration at loadstep: ", loadstep)

        # # # # # # # # # #
        # Newton solve.   #
        # # # # # # # # # #
        while (iter < imax && residual > TOL) || iter < 2
            iter += 1
            a += Δa
            assemGlobal!(K, Fᵢₙₜ, dh, mp, t, a, coord, enod, ε)

            solveq!(Δa, K, -Fᵢₙₜ, bcdof, bcval)
            bcval = 0 * bcval
            res = Fᵢₙₜ - Fₑₓₜ
            res[bcdof] = 0 * res[bcdof]
            residual = norm(res, 2)
            println("Iteration: ", iter, " Residual: ", residual)

            postprocess_opt(a, dh, "contact_mesh" * string(loadstep))
            #postprocess_opt(Fᵢₙₜ, dh, "contact_mesh" * string(loadstep))
            σx,σy = StressExtract(dh,a,mp)
            vtk_grid("contact"*string(loadstep), dh) do vtkfile
            vtk_point_data(vtkfile, dh, a) # displacement field
            vtk_point_data(vtkfile, σx, "σx")
            vtk_point_data(vtkfile, σy, "σy")
            end
        end
    end
    fill!(Fₑₓₜ, 0.0)
    Fₑₓₜ[bcdof] = -Fᵢₙₜ[bcdof]
    return a, dh, Fₑₓₜ, Fᵢₙₜ, K
end

a = zeros(dh.ndofs.x)
postprocess_opt(a, dh, "contact_mesh" * string(1))
# Create dictionaries that are needed for the Mortar2D package
elements, element_types, slave_elements, slave_element_ids, master_element_ids, coords = create_contact_list(dh, Γs, Γm, coord)
normals  = Mortar2D.calculate_normals(elements, element_types, coords)
segments = Mortar2D.calculate_segments(slave_element_ids, master_element_ids, elements, element_types, coords, normals)
# Assemble D and M matrices and the slave and master dofs corresponding to the mortar segmentation
slave_dofs, master_dofs, D, M = Mortar2D.calculate_mortar_assembly(elements, element_types, coords, slave_element_ids, master_element_ids)


g = gap_function(X)

# Solve
a, dh, Fₑₓₜ, Fᵢₙₜ, K = solver(dh, coord);

# Visualize using e.g. paraview (import image) or FerriteViz
# postprocess_opt(Fᵢₙₜ, dh, "contact_mesh");

# function to extract stressfield..