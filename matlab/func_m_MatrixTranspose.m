function B = func_m_MatrixTranspose(A)
% B.eij = A.eji;

nfields = length(fieldnames(A));
N = fix(sqrt(nfields));

for i = 1:N
    for j = 1:N
        eval(sprintf('B.e%d%d = A.e%d%d;',i,j,j,i));
    end
end

end