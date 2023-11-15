function dV = func_m_VectorCurl(V,x,y,z)

[X,Y,Z] = meshgrid(x,y,z);
[dV.e1,dV.e2,dV.e3,~] = curl(X,Y,Z,V.e1,V.e2,V.e3);

end
