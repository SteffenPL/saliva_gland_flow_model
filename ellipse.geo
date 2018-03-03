//+
r_1 = DefineNumber[ 0.1, Name "Parameters/r_1" ];
//+
r_2 = DefineNumber[ 0.01, Name "Parameters/r_2" ];
//+
SetFactory("OpenCASCADE");
Ellipse(1) = {0, 0, 0, r_1, r_2, 0, 2*Pi};
