function C = func_m_VectorSum(A,B)
% C.ei = A.ei + B.ei;

nfields = length(fieldnames(A));
for i = 1:nfields
        eval(sprintf('C.e%d = A.e%d + B.e%d;',i,i,i));
end
end