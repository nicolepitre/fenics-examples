cl__1 = 1;
Point(1) = {0, 0, 0, 1};
Point(2) = {4.211, 0, 0, 1};
Point(3) = {4.211, 0.15, 0, 1};
Point(4) = {0, 0.15, 0, 1};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = {5};
Physical Surface(7) = {6};

Transfinite Line{1} = 100 Using Progression 1;
Transfinite Line{2} = 2 Using Progression 1;
Transfinite Line{3} = 100 Using Progression 1;
Transfinite Line{4} = 2 Using Progression 1;
Transfinite Surface{6} = {1,2,3,4};
Recombine Surface {6};

Physical Line(111) = {1};
Physical Line(112) = {2};
Physical Line(113) = {3};
Physical Line(114) = {4};
