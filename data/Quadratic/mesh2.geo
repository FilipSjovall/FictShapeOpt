// Gmsh project created on Thu Mar 23 09:46:50 2023
SetFactory("Built-in");
//+
Point(1) = {0.5, 0, 0, 0.1};
//+
Point(2) = {1, 0, 0, 0.1};
//+
Point(3) = {1, 1, 0, 0.1};
//+ 
Point(4) = {0, 1, 0, 0.1};
//+ 
Point(5) = {0, 0.5, 0, 0.1};
//+ 
Point(6) = {0.5, 0.5, 0, 0.1};
//+
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};
//+
Line loop(1)     = {1,2,3,4,5,6};

//+
Plane Surface(1) = {1};

Recursive Delete {
  Curve{1}; Curve{2}; Curve{3}; Curve{4}; Curve{5}; Curve{6}; 
}

//+
Physical Surface("domain") = {1};
/*
Physical Line("a") = {1};
Physical Line("b") = {2};
Physical Line("c") = {3};
Physical Line("d") = {4};
Physical Line("e") = {5};
Physical Line("f") = {6};
*/

