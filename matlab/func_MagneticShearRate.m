function S = func_MagneticShearRate(Bx,By,Bz,x,y,z,FixTrace)

DB = func_MagneticJacobian(Bx,By,Bz,x,y,z,FixTrace);

% Symmetric part of DB
S_a_12 = 0.5*(DB.e12+DB.e21);
S_b_13 = 0.5*(DB.e13+DB.e31);
S_c_23 = 0.5*(DB.e23+DB.e32);

S = sqrt(2*(S_a_12.^2+S_b_13.^2+S_c_23.^2)+DB.e11.^2+DB.e22.^2+DB.e33.^2);


end