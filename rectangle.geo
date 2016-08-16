   /* Rectangle Geometry File */

    gridsize = 2.0; // prescribed mesh element size
    
    Lx = 10.0; // width
    Ly = 10.0;   // height
    
    Point(1) = {0, 0, 0, gridsize};
    Point(2) = {0, Ly, 0, gridsize};
    Point(3) = {Lx, Ly, 0, gridsize};
    Point(4) = {Lx, 0, 0, gridsize};
 
    Line(1) = {1, 2}; 
    Line(2) = {2, 3};
    Line(3) = {3, 4};
    Line(4) = {4, 1};

    Line Loop(5) = {1, 2, 3, 4};
    rectangle = news;
    Plane Surface(rectangle) = {5};
  
    Physical Line("Left", 1) = {1};
    Physical Line("Top", 2) = {2};
    Physical Line("Right", 3) = {3};
    Physical Line("Bottom", 4) = {4};

    Physical Surface("Rectangle", 5) = {rectangle};
    