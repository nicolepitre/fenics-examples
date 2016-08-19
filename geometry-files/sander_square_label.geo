lc = 1e-1;

Point(1) = {0, 0, 0, lc};
Point(2) = {1, 0,  0, lc} ;
Point(3) = {1, 1, 0, lc} ;
Point(4) = {0, 1, 0, lc} ;

Line(1) = {1,2} ;
Line(2) = {3,2} ;
Line(3) = {3,4} ;
Line(4) = {4,1} ;

Line Loop(5) = {4,1,-2,3} ;
Plane Surface(6) = {5} ;

Physical Point(1) = {1,2} ;

Physical Line(200) = {1,2} ;
Physical Line(300) = {3,4} ;
Physical Surface(1000) = {6} ;
