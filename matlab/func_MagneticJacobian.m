function DB = func_MagneticJacobian(Bx,By,Bz,x,y,z,FixTrace)

dx = mean(diff(x));
dy = mean(diff(y));
dz = mean(diff(z));

% DB
[DB.e11,DB.e12,DB.e13] = gradient(Bx,dx,dy,dz);
[DB.e21,DB.e22,DB.e23] = gradient(By,dx,dy,dz);
[DB.e31,DB.e32,DB.e33] = gradient(Bz,dx,dy,dz);

% Fix trace
if FixTrace ~= 0
    trDB = 1/3*(DB.e11 + DB.e22 + DB.e33);
    DB.e11 = DB.e11 - trDB;
    DB.e22 = DB.e22 - trDB;
    DB.e33 = DB.e33 - trDB;
end
end
