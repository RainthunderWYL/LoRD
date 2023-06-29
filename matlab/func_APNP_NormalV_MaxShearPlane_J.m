function [ex,ey,ez] = func_APNP_NormalV_MaxShearPlane_J(Bx,By,Bz,x,y,z,Parameters)
% For Real Roots, e3 is the direction of slowest changing direction: Median eigenvalue
% For Complex roots, e3 is the direction of real eigenvector

FixTrace = Parameters.APNP_FixTrace;

dx = mean(diff(x));
dy = mean(diff(y));
dz = mean(diff(z));

NX = size(Bx);
ndata = numel(Bx);

% DB
% DB = func_MagneticGradient(Bx,By,Bz,x,y,z,1);
[DB_11,DB_12,DB_13] = gradient(Bx,dx,dy,dz);
[DB_21,DB_22,DB_23] = gradient(By,dx,dy,dz);
[DB_31,DB_32,DB_33] = gradient(Bz,dx,dy,dz);

% Reshape to 1D array
DB_11 = reshape(DB_11,[1,ndata]); % a
DB_22 = reshape(DB_22,[1,ndata]); % b
DB_33 = reshape(DB_33,[1,ndata]); % -a-b
DB_12 = reshape(DB_12,[1,ndata]); % c
DB_21 = reshape(DB_21,[1,ndata]); % d
DB_13 = reshape(DB_13,[1,ndata]); % e
DB_31 = reshape(DB_31,[1,ndata]); % f
DB_23 = reshape(DB_23,[1,ndata]); % g
DB_32 = reshape(DB_32,[1,ndata]); % h

% Fix trace
if FixTrace ~= 0
    trDB = 1/3*(DB_11 + DB_22 + DB_33);
    DB_11 = DB_11 - trDB;
    DB_22 = DB_22 - trDB;
    DB_33 = DB_33 - trDB;
    clear trDB;
end

[l1,~,~,infuncV1,infuncV2,infuncV3] = ...
    func_Eig_3DTraceLess([DB_11;DB_21;DB_31;DB_12;DB_22;DB_32;DB_13;DB_23;DB_33]);

ex = infuncV1(l1,DB_11,DB_22,DB_12,DB_21,DB_13,DB_31,DB_23,DB_32);
ey = infuncV2(l1,DB_11,DB_22,DB_12,DB_21,DB_13,DB_31,DB_23,DB_32);
ez = infuncV3(l1,DB_11,DB_22,DB_12,DB_21,DB_13,DB_31,DB_23,DB_32);
enorm = sqrt(ex.^2+ey.^2+ez.^2);

ex = reshape(ex./enorm,NX);
ey = reshape(ey./enorm,NX);
ez = reshape(ez./enorm,NX);

end
