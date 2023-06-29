function Ac = func_m_VectorScaleProd(A,c)
% C.ei = A.ei .* B.ei;

nfields = length(fieldnames(A));
for i = 1:nfields
        eval(sprintf('Ac.e%d = A.e%d .* c;',i,i));
end
end