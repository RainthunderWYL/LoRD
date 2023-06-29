function A = func_m_GenMatrixStruct(M)
% Aij is a structure with fields: eij

[m,n,l] = size(M);

for i = 1:m
    for j = 1:n
        eval(sprintf('A.e%d%d = reshape(M(i,j,:),[l,1]);',i,j));
    end
end

end