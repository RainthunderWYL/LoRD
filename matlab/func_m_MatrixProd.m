function C = func_m_MatrixProd(A,B)
% Calculate Cij = Aim*Bmj
% A and B are nxn matrices.
% Aij is a structure with fields: eij

nfields = length(fieldnames(A));
N = fix(sqrt(nfields));
datasize = size(A.e11);

for i = 1:N
    for j = 1:N
        eval(sprintf('C.e%d%d = zeros(datasize);',i,j));
        for m = 1:N
            eval(sprintf('C.e%d%d = C.e%d%d + A.e%d%d .* B.e%d%d;',i,j,i,j,i,m,m,j));
        end
    end
end

end