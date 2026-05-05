SetFactory("OpenCASCADE");
Merge "phantom_3_rebuild_r2.STEP";
//+
Field[1] = Distance;
//+
Field[1].CurvesList = {11, 13, 17, 18};
//+
Field[1].Sampling = 100;
//+
Field[2] = Param;
//+
Delete Field [2];
//+
Field[2] = MathEval;
