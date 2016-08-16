// Cube Geometry File - periodic in x, y, and z directions

gridsize = 0.2;

Point(1) = {0, 0, 0, gridsize};
Point(2) = {1, 0, 0, gridsize};
Point(3) = {0, 1, 0, gridsize};
Point(4) = {1, 1, 0, gridsize};
Point(5) = {1, 1, 1, gridsize};
Point(6) = {1, 0, 1, gridsize};
Point(7) = {0, 1, 1, gridsize};
Point(8) = {0, 0, 1, gridsize};

Line(1) = {3, 7};
Line(2) = {7, 5};
Line(3) = {5, 4};
Line(4) = {4, 3};
Line(5) = {3, 1};
Line(6) = {2, 4};
Line(7) = {2, 6};
Line(8) = {6, 8};
Line(9) = {8, 1};
Line(10) = {1, 2};
Line(11) = {8, 7};
Line(12) = {6, 5};

Line Loop(13) = {7, 8, 9, 10};
bottom = news ;
Plane Surface(bottom) = {13};

Line Loop(15) = {-10, -5, -4, -6};
back = news ;
Plane Surface(back) = {15};

Line Loop(17) = {1, 2, 3, 4};
top = news ;
Plane Surface(top) = {17};

Line Loop(19) = {12, -2, -11, -8};
front = news ;
Plane Surface(front) = {19};

Line Loop(21) = {6, -3, -12, -7};
right = news ;
Plane Surface(right) = {21};

Line Loop(23) = {11, -1, 5, -9};
left = news;
Plane Surface(left) = {23};

Surface Loop(25) = {bottom, right, front, top, back, left};
vol = newv;
Volume(vol) = {25};

Physical Line(0) = {1};
Physical Line(1) = {2};
Physical Line(2) = {3};
Physical Line(3) = {4};
Physical Line(4) = {5};
Physical Line(5) = {6};
Physical Line(6) = {7};
Physical Line(7) = {8};
Physical Line(8) = {9};
Physical Line(9) = {10};
Physical Line(10) = {11};
Physical Line(11) = {12};

Physical Surface("Surface 0", 111) = {bottom};
Physical Surface("Surface 1", 222) = {back};
Physical Surface("Surface 2", 333) = {top};
Physical Surface("Surface 3", 444) = {front};
Physical Surface("Surface 4", 555) = {right};
Physical Surface("Surface 5", 666) = {left};
  
Physical Volume ("Volume", 01) = {vol};

Periodic Surface bottom {-9, -8, -7, -10} = top   {1, 2, 3, 4} ;
Periodic Surface back   {6, 4, 5, 10}     = front {12, -2, -11, -8} ;
Periodic Surface left   {-5, 1, -11, 9}   = right {6, -3, -12, -7} ;
