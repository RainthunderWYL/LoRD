function C = func_m_VectorCross(A,B)
% for 3 array
% C.ei = eijk .* A.ej .* B.ek;

C.e1 = A.e2 .* B.e3 - A.e3 .* B.e2;
C.e2 = A.e3 .* B.e1 - A.e1 .* B.e3;
C.e3 = A.e1 .* B.e2 - A.e2 .* B.e1;

end