#load_files()

grid1 = createBoxMesh("box_1",0.0,0.0,1.0,1.0,0.1)

dh = DofHandler(grid1)
add!(dh, :u, 2)
close!(dh)

ip      = Lagrange{2, RefTetrahedron, 1}()
qr      = QuadratureRule{2, RefTetrahedron}(1)
qr_face = QuadratureRule{1, RefTetrahedron}(1)
cv      = CellVectorValues(qr, ip)
fv      = FaceVectorValues(qr_face, ip)

coord, enod = getTopology(dh)


addfaceset!(dh.grid,"Γₜ", x -> x[1] ≈ 1.0)
Γt = getfaceset(dh.grid, "Γₜ")


