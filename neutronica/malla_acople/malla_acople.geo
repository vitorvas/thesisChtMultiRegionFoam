// Gmsh project created on Fri Sep 12 10:30:40 2014
Point(1) = {-0.05, -0.05, 0, 1.0};
Point(2) = {0, -0.05, 0, 1.0};
Point(3) = {-0.05, 0, 0, 1.0};
Point(4) = {0, 0, 0, 1.0};
Line(1) = {1, 2};
Line(2) = {2, 4};
Line(3) = {4, 3};
Line(4) = {3, 1};
Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = {5};
Translate {0, 0, 0.4} {
  Duplicata { Surface{6}; }
}
Line(12) = {1, 5};
Line(13) = {2, 6};
Line(14) = {4, 10};
Line(15) = {3, 14};
Line Loop(16) = {15, 11, -12, -4};
Plane Surface(17) = {16};
Line Loop(18) = {2, 14, -9, -13};
Plane Surface(19) = {18};
Line Loop(20) = {1, 13, -8, -12};
Plane Surface(21) = {20};
Line Loop(22) = {3, 15, -10, -14};
Plane Surface(23) = {22};
Surface Loop(24) = {7, 21, 6, 19, 23, 17};
Volume(25) = {24};
Translate {0.05, 0, 0} {
  Duplicata { Volume{25}; }
}
Translate {0, 0.05, 0} {
  Duplicata { Volume{25, 26}; }
}
Translate {0, 0, 0.4} {
  Duplicata { Volume{26, 25, 79, 48}; }
}
Translate {0, 0, -0.08} {
  Point{384, 192, 188, 187, 196, 388, 292, 283, 484};
}
Physical Surface("<-x>") = {220, 74, 158, 17};
Physical Surface("<+x>") = {179, 95, 117, 42};
Physical Surface("<-y>") = {138, 21, 107, 32};
Physical Surface("<+y>") = {215, 69, 184, 100};
Physical Surface("<-z>") = {37, 6, 90, 59};
Physical Surface("<+z>") = {102, 133, 164, 195};
Physical Volume("<geom>") = {101, 163, 132, 194, 26, 79, 25, 48};
Characteristic Length {192, 188, 187, 196, 283, 292, 384, 388, 484, 62, 58, 20, 16, 14, 10, 6, 5, 154, 94, 90, 186, 52, 32, 4, 3, 2, 1} = 0.025;
