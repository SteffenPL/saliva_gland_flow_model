// Gmsh project created on Sat Mar  3 13:21:35 2018
//+
lc = DefineNumber[ 0.1, Name "Parameters/lc" ];
//+
Point(1) = {1, 0.2, 0, lc};
//+
Point(2) = {1, 0, 0, lc};
//+
Point(3) = {0.7, 0.3, 0, lc};
//+
Point(4) = {0.7, 0.7, 0, lc};
//+
Point(5) = {0.5, 0.8, 0, lc};
//+
Point(6) = {-0.3, 0.8, 0, lc};
//+
Point(7) = {-0.4, 0.7, 0, lc};
//+
Point(8) = {-0.2, -0.7, 0, lc};
//+
Point(9) = {0.4, -0.9, 0, lc};
//+
Point(10) = {0.8, -0.9, 0, lc};
//+
Point(11) = {0.9, -0.7, 0, lc};
//+
Point(12) = {0.7, -0.4, 0, lc};
//+
Point(13) = {0.8, -0.1, 0, lc};
//+
BSpline(1) = {1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 2};
//+
Line(2) = {2, 1};
//+
Line Loop(1) = {1, 2};
//+
Plane Surface(1) = {1};
