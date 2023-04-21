//+
SetFactory("OpenCASCADE");
// Set the meshing algorithm
Mesh.Algorithm = 8;

// Define the geometries to be meshed
Rectangle(1) = {0, 0, 0, 1, 1, 0};
Circle(2) = {0.5, 1.3, 0, 0.25, 0, 2*Pi};

// Define the element type to be used in the mesh
Mesh.ElementOrder = 1;
Mesh.ElementType = 2;

// Define physical groups for the geometries
Physical Point("Rectangle") = {1};
Physical Point("Circle") = {2};

// Generate the mesh for the rectangle
Mesh.SubdivisionAlgorithm = 1;
Mesh 1;
Save "rectangle.inp";

// Generate the mesh for the circle
Mesh.SubdivisionAlgorithm = 1;
Mesh 2;
Save "circle.inp";

