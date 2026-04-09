SetFactory("OpenCASCADE");
Merge "phantom_3_mod2.step";
//+
Physical Surface(134) = {61, 34, 27, 29, 35, 39, 40, 26, 46, 38, 41, 33, 32, 31, 15, 36, 5, 30, 37, 7, 6, 14, 47, 23, 19, 59, 22, 13, 45, 44, 10, 11, 12, 16, 21, 18, 28, 56, 4, 57, 58, 48, 24, 17, 49, 3, 20, 25, 8, 43, 42, 60, 1};
//+
Field[1] = Distance;
//+
Field[1].CurvesList = {88, 83, 82, 95, 78, 79, 87, 86, 73, 74};
//+
Field[1].CurvesList = {88, 83, 82, 95, 78, 79, 87, 86, 73, 74, 18, 15, 51, 56};
//+
Field[1].CurvesList = {88, 83, 82, 95, 78, 79, 87, 86, 73, 74, 18, 15, 51, 56, 55, 43, 14, 8, 7, 9};
//+
Field[2] = Distance;
//+
Field[2].CurvesList = {97, 42};
//+
Field[3] = MathEval;
//+
Field[3].F = "F1^1.2*0.75 + 30";
//+
Background Field = 3;
