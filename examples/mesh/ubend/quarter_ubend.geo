// Quarter U-bend duct: straight inlet -> 90 degree bend -> straight outlet.
// Structured (transfinite) all-hex mesh throughout.

SetFactory("Built-in");

// parameters
R      = 3.0;  // bend radius to duct centerline
w      = 1.0;  // duct cross-section side length
L1     = 2.0;  // inlet straight length
L2     = 2.0;  // outlet straight length

nCross  = 5;   // divisions across each cross-section edge
nInlet  = 6;   // divisions along inlet leg
nBend   = 18;  // divisions along the bend
nOutlet = 6;   // divisions along outlet leg

// --- inlet cross-section (at y = -L1) ---
Point(1) = {R - w/2, -L1, 0, 1.0};
Point(2) = {R + w/2, -L1, 0, 1.0};
Point(3) = {R + w/2, -L1, w, 1.0};
Point(4) = {R - w/2, -L1, w, 1.0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Transfinite Line {1, 3} = nCross + 1;
Transfinite Line {2, 4} = nCross + 1;
Transfinite Surface {1};
Recombine Surface {1};

// --- inlet straight leg ---
inlet[] = Extrude {0, L1, 0} { Surface{1}; Layers{nInlet}; Recombine; };

// --- 90 degree bend ---
bend[] = Extrude { {0,0,1}, {0,0,0}, Pi/2 } { Surface{inlet[0]}; Layers{nBend}; Recombine; };

// --- outlet straight leg ---
outlet[] = Extrude {-L2, 0, 0} { Surface{bend[0]}; Layers{nOutlet}; Recombine; };

Physical Surface("inlet") = {1};
Physical Surface("outlet") = {outlet[0]};
Physical Surface("wall") = {inlet[2], inlet[3], inlet[4], inlet[5],
                             bend[2], bend[3], bend[4], bend[5],
                             outlet[2], outlet[3], outlet[4], outlet[5]};
Physical Volume("duct") = {inlet[1], bend[1], outlet[1]};

Mesh.ElementOrder = 1;
