   /* 3D Cape Geometry File */

    gridsize = 0.05; // prescribed mesh element size
    meshThickness = 0.1;
    
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
    
    surfaceVector[] = Extrude {0, 0, meshThickness} {
     Surface{6};
     Layers{1};
     Recombine;
    };
    /* surfaceVector contains in the following order:
    [0]	- front surface (opposed to source surface)
    [1] - extruded volume
    [2] - bottom surface (belonging to 1st line in "Line Loop (6)")
    [3] - right surface (belonging to 2nd line in "Line Loop (6)")
    [4] - top surface (belonging to 3rd line in "Line Loop (6)")
    [5] - left surface (belonging to 4th line in "Line Loop (6)") */
    Physical Surface("front", 1) = surfaceVector[0];
    Physical Volume("internal", 2) = surfaceVector[1];
    Physical Surface("bottom", 3) = surfaceVector[2];
    Physical Surface("right", 4) = surfaceVector[3];
    Physical Surface("top", 5) = surfaceVector[4];
    Physical Surface("left", 6) = surfaceVector[5];
    Physical Surface("back", 7) = {6}; // from Plane Surface (6) ...
    
