function B = func_m_TensorSubset(A,b)

fname = fieldnames(A);
nfields = length(fname);

for i = 1:nfields
    eval(sprintf('B.%s = A.%s(b);',fname{i},fname{i}));
end

end