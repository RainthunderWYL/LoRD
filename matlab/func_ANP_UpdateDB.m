function NPInfo = func_ANP_UpdateDB(NPInfo_Raw,B1,B2,B3,x1,x2,x3)
% Update DB by interpol DB

nNPs = size(NPInfo_Raw,1);
NPInfo = NPInfo_Raw;

[X1,X2,X3] = meshgrid(x1,x2,x3);
dx = mean(diff(x1));
dy = mean(diff(x2));
dz = mean(diff(x3));
% DB
[DB_11,DB_12,DB_13] = gradient(B1,dx,dy,dz);
[DB_21,DB_22,DB_23] = gradient(B2,dx,dy,dz);
[DB_31,DB_32,DB_33] = gradient(B3,dx,dy,dz);

DB_ALL = zeros(nNPs,9);
x_all = NPInfo_Raw(:,1);
y_all = NPInfo_Raw(:,2);
z_all = NPInfo_Raw(:,3);

DB_All(:,1) = interp3(X1,X2,X3,DB_11,x_all,y_all,z_all,'linear');
DB_All(:,2) = interp3(X1,X2,X3,DB_21,x_all,y_all,z_all,'linear');
DB_All(:,3) = interp3(X1,X2,X3,DB_31,x_all,y_all,z_all,'linear');
DB_All(:,4) = interp3(X1,X2,X3,DB_12,x_all,y_all,z_all,'linear');
DB_All(:,5) = interp3(X1,X2,X3,DB_22,x_all,y_all,z_all,'linear');
DB_All(:,6) = interp3(X1,X2,X3,DB_32,x_all,y_all,z_all,'linear');
DB_All(:,7) = interp3(X1,X2,X3,DB_13,x_all,y_all,z_all,'linear');
DB_All(:,8) = interp3(X1,X2,X3,DB_23,x_all,y_all,z_all,'linear');
DB_All(:,9) = interp3(X1,X2,X3,DB_33,x_all,y_all,z_all,'linear');

NPInfo(:,5:13) = DB_All;
end