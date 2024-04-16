lc = 0.002;

Point(1) = {0,  0,   0,  lc};
Point(2) = {0.083,  0.048,  0,  lc};
Point(3) = {-0.083,   0.048,  0,  lc};

Line(1) = {1,  2};
Line(2) = {2,  3};
Line(3) = {3,  1};

Curve Loop(1) = {1,  2,  3};

Plane Surface(1) = {1};

Mesh.SaveAll = 1; 
Mesh 2;
Save "wedge.vtk";