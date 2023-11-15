function J = func_m_VectorJacobian(V,x,y,z,FixTrace)

dx = mean(diff(x));
dy = mean(diff(y));
dz = mean(diff(z));

% J
[J.e11,J.e12,J.e13] = gradient(V.e1,dx,dy,dz);
[J.e21,J.e22,J.e23] = gradient(V.e2,dx,dy,dz);
[J.e31,J.e32,J.e33] = gradient(V.e3,dx,dy,dz);

% Fix trace
if FixTrace ~= 0
    trJ = 1/3*(J.e11 + J.e22 + J.e33);
    J.e11 = J.e11 - trJ;
    J.e22 = J.e22 - trJ;
    J.e33 = J.e33 - trJ;
end
end
