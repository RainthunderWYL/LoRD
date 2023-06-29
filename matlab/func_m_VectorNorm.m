function Anorm = func_m_VectorNorm(A)
% C.ei = A.ei .* B.ei;

nfields = length(fieldnames(A));
Anorm = cast(zeros(size(A.e1)),'like',A.e1);
for i = 1:nfields
        eval(sprintf('Anorm = Anorm + A.e%d .* A.e%d;',i,i));
end
Anorm = sqrt(Anorm);

end