function C = func_m_VectorDot(A,B)
% C.ei = A.ei .* B.ei;

nfields = length(fieldnames(A));
C = zeros(size(A.e1));
for i = 1:nfields
        eval(sprintf('C = C + A.e%d .* B.e%d;',i,i));
end

end