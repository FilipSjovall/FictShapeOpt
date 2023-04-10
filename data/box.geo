// Gmsh project created on Sat Apr 08 01:13:21 2023
SetFactory("OpenCASCADE");
//+
Point(1) = {1, 1, 0, 1.0};
//+
Point(2) = {0, 1, 0, 1.0};
//+
Point(3) = {0, 0, 0, 1.0};
//+
Point(4) = {1, 0, 0, 1.0};
//+
Line(1) = {2, 1};
//+
Line(2) = {1, 4};
//+
Line(3) = {4, 3};
//+
Line(4) = {2, 3};
//+
Curve Loop(1) = {4, -3, -2, -1};
//+
Plane Surface(1) = {1};
//+
Physical Surface(5) = {1};
