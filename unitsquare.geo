   /* Unit Square Geometry File */

    gridsize = 0.05; // prescribed mesh element size
    
    Lx = 1.0; // length measured along bottom edge
    Ly = 1.0;   // height
    
    Point(1) = {0, 0, 0, gridsize};
    Point(2) = {0, 1, 0, gridsize};
    Point(3) = {1, 1, 0, gridsize};
    Point(4) = {1, 0, 0, gridsize};
 
    Line(1) = {1, 2};
    Line(2) = {2, 3};
    Line(3) = {3, 4};
    Line(4) = {4, 1};

    Line Loop(5) = {1, 2, 3, 4};
    Plane Surface(6) = {5};