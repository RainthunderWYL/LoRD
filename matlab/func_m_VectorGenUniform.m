function V = func_m_VectorGenUniform(v,SIZE)
len = length(v);
for i = 1:len
    eval(sprintf('V.e%d = v(%d)*ones(SIZE);',i,i));
end

end