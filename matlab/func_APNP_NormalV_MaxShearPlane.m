function [ex,ey,ez] = func_APNP_NormalV_MaxShearPlane(Bx,By,Bz,x,y,z,Parameters)
% DB_{ij} = \partial_jB_i
% DBS_{ij} = 0.5(DB_{ij} + DB_{ji})
% DBA_{ij} = 0.5(DB_{ij} - DB_{ji})

% UseCurrent = Parameters.UseCurrent;

dx = mean(diff(x));
dy = mean(diff(y));
dz = mean(diff(z));

% DB
[DB_11,DB_12,DB_13] = gradient(Bx,dx,dy,dz);
[DB_21,DB_22,DB_23] = gradient(By,dx,dy,dz);
[DB_31,DB_32,DB_33] = gradient(Bz,dx,dy,dz);

% Symmetric part of DB
S_a_12 = 0.5*(DB_12+DB_21);
S_b_13 = 0.5*(DB_13+DB_31);
S_c_23 = 0.5*(DB_23+DB_32);

% Snorm = sqrt(2*(S_a_12.^2+S_b_13.^2+S_c_23.^2)+DB_11.^2+DB_22.^2+DB_33.^2);

clear DB_12 DB_13 DB_21 DB_23 DB_31 DB_32;

% "Guide field direction: median eigen-value of S"
[~,l2,~,infuncV1,infuncV2,infuncV3] = func_Eig_3DSymm(S_a_12,S_b_13,S_c_23,DB_11,DB_22,DB_33);

eG_x = infuncV1(l2,S_a_12,S_b_13,S_c_23,DB_11,DB_22,DB_33);
eG_y = infuncV2(l2,S_a_12,S_b_13,S_c_23,DB_11,DB_22,DB_33);
eG_z = infuncV3(l2,S_a_12,S_b_13,S_c_23,DB_11,DB_22,DB_33);
eG_norm = sqrt(eG_x.^2+eG_y.^2+eG_z.^2);

ex = eG_x./eG_norm;
ey = eG_y./eG_norm;
ez = eG_z./eG_norm;
end