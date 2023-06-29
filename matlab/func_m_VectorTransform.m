function Y = func_m_VectorTransform(M,V)
% Y.ei = M.eij .* V.ej

nfields = length(fieldnames(V));
e_size = size(V.e1);

for i = 1:nfields
    eval(sprintf('Y.e%d = zeros(e_size);',i));
    for j = 1:nfields
        eval(sprintf('Y.e%d = Y.e%d + M.e%d%d .* V.e%d;',i,i,i,j,j));
    end
end
end