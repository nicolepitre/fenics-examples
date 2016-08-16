   /* Cape Geometry File */

    gridsize = 0.02; // prescribed mesh element size
    
    Lx = 1.0; // length measured along bottom edge
    Ly = 1.0;   // height
    
    Point(1) = {0, 0, 0, gridsize};
    Point(2) = {0, 1, 0, gridsize};
    Point(3) = {1, 1, 0, gridsize};
    Point(4) = {1, 0, 0, gridsize};
 
    Line(1) = {3, 4};
    Line(2) = {4, 1};
    Line(3) = {1, 2};

    x0 = 0.5;   // centre of the Gaussian
    W  = 0.2; // Width of the cape
    H  = 0.1; // Height of the cape
    Np = 51;  // Number of points across Northern Wall
    
    pList[0] = 2; // First point label (top-left point of the inlet region)
    nPoints = Np; // # of discretization points
    For i In {1 : nPoints}
      x = Lx*i/(nPoints + 1);
      pList[i] = newp;
      Point(pList[i]) = {x,
                    ( Ly - H*Exp((-((x-x0)/W)^2))), 
                    0, gridsize};
    EndFor
    pList[nPoints+1] = 3; // Last point label (top-right point of the outlet region)
 
    Spline(newl) = pList[]; // Interpolate the points with a spline
 
    Line Loop(5) = {1, 2, 3, 4};
    Plane Surface(6) = {5};

    Physical Line("Right", 1) = {1};
    Physical Line("Bottom", 2) = {2};
    Physical Line("Left", 3) = {3};
    Physical Line("Top", 4) = {4};
    Physical Surface("Cape", 5) = {6};
