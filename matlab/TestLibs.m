clear;
M1 = rand(3);
M2 = rand(3);

A = func_m_GenMatrixStruct(M1);
B = func_m_GenMatrixStruct(M2);

C = func_m_MatrixProd(A,B);

V1.e1 = [1;2];
V1.e2 = [3;4];
V1.e3 = [5;6];

V2.e1 = [1;2];
V2.e2 = [3;4];
V2.e3 = [5;6];

