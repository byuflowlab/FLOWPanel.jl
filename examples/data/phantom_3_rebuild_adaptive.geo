SetFactory("OpenCASCADE");
Merge "phantom_3_rebuild_r2.step";   // <-- your STEP filename

// =====================================================================
// PARAMETERS
// =====================================================================
S0_LE  = 0.8;   // element size right at the leading edge
S0_TE  = 0.4;   // element size right at the trailing edge
S1     = 3.0;    // background element size away from refined edges
D0_LE  = 4.0;    // transition distance over which LE size grows S0->S1
D0_TE  = 4.0;    // transition distance over which TE size grows S0->S1

// =====================================================================
// PHYSICAL SURFACE  (replace with the surface IDs of your two blades)
// =====================================================================
Physical Surface("rotor") = {5, 7, 9, 11, 6, 10};
//Physical Surface("rotor") = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 /* surface IDs from STEP, e.g. 1, 2, ... */ };

// =====================================================================
// DISTANCE FIELDS
// =====================================================================
// Distance to LE curves (BSplines 11 and 17, one per blade)
Field[1] = Distance;
Field[1].CurvesList = {11, 17};
Field[1].NumPointsPerCurve = 1000;

// Distance to TE curves (BSplines 13 and 18, one per blade)
Field[2] = Distance;
Field[2].CurvesList = {13, 18};
Field[2].NumPointsPerCurve = 1000;

// =====================================================================
// SMOOTH-STEP SIZE TRANSITION
//
//   size(d) = S0 + (S1 - S0) * H(min(1, d/D0))
//   H(t)    = 3 t^2 - 2 t^3        (smoothstep, C^1 at both ends)
//
//   d=0      -> size = S0  (refined)
//   d>=D0    -> size = S1  (background)
// =====================================================================

// LE size field (uses Field[1])
Field[3] = MathEval;
Field[3].F = Sprintf(
  "%g + (%g) * (3*(min(1, F1/%g))^2 - 2*(min(1, F1/%g))^3)",
  S0_LE, S1 - S0_LE, D0_LE, D0_LE);

// TE size field (uses Field[2])
Field[4] = MathEval;
Field[4].F = Sprintf(
  "%g + (%g) * (3*(min(1, F2/%g))^2 - 2*(min(1, F2/%g))^3)",
  S0_TE, S1 - S0_TE, D0_TE, D0_TE);

// =====================================================================
// COMBINE: take the minimum size from LE and TE refinements
// (so a panel near both LE and TE picks up the tighter constraint)
// =====================================================================
Field[5] = Min;
Field[5].FieldsList = {3, 4};
Background Field = 5;

// =====================================================================
// MESHER SETTINGS
// =====================================================================
Mesh.Algorithm = 6;          // Frontal-Delaunay (good for surface meshes)
Mesh.MeshSizeFromCurvature = 0;
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;

// Optional: write ASCII v4 .msh on save
// Mesh.MshFileVersion = 4.1;
// Mesh.Binary = 0;

