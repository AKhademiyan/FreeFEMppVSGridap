// Wavy channel geometry with sinusoidal walls
// Parameters
L = 1.0;      // Length of the channel
H = 0.5;      // Height of the channel
A = 0.1;      // Amplitude of the sine wave
lambda = 0.5; // Wavelength
n_waves = 2;  // Number of complete waves
n_points = 20; // Number of points per wave (reduced for simplicity)

// Mesh parameters
lc = 0.0055;    // Characteristic length for mesh
lc_wall = lc/2; // Finer mesh near walls

// Create points for the bottom wall (sine wave)
For i In {0:n_points}
    x = i * L/n_points;
    y = -H/2 + A * Sin(2*Pi*x/lambda);
    Point(i+1) = {x, y, 0, lc_wall};
EndFor

// Create points for the top wall (negative sine wave)
For i In {0:n_points}
    x = i * L/n_points;
    y = H/2 - A * Sin(2*Pi*x/lambda);
    Point(i+n_points+2) = {x, y, 0, lc_wall};
EndFor

// Create lines for bottom wall
For i In {1:n_points}
    Line(i) = {i, i+1};
EndFor

// Create lines for top wall
For i In {1:n_points}
    Line(i+n_points) = {i+n_points+1, i+n_points+2};
EndFor

// Create inlet and outlet lines (straight)
Line(2*n_points+1) = {1, n_points+2};    // Inlet
Line(2*n_points+2) = {n_points+1, 2*n_points+2};  // Outlet

// Create surface
Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 42, -21, -22, -23, -24, -25, -26, -27, -28, -29, -30, -31, -32, -33, -34, -35, -36, -37, -38, -39, -40, -41};
Plane Surface(1) = {1};

// Define physical groups
Physical Curve("wall_down") = {1:n_points};
Physical Curve("wall_up") = {n_points+1:2*n_points};
Physical Curve("inlet") = {2*n_points+1};
Physical Curve("outlet") = {2*n_points+2};
Physical Surface("fluid") = {1};

// Generate mesh
Mesh 2;
SetOrder 2; 