function b = func_iszero(X)
b = abs(X)<=eps('single');
% b = abs(X) < 1e-3;
end