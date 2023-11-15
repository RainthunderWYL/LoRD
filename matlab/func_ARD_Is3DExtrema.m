function Is3DExtrema = func_ARD_Is3DExtrema(...
                           R0,Npara, x1,x2,x3,NparaAll)

dL = mean([mean(diff(x1)),mean(diff(x3)),mean(diff(x3))]);
[X1,X2,X3] = meshgrid(x1,x2,x3);

% Relative coord of 6 points surounding R0
e_size = size(R0.e1);
dr_1m2030 = func_m_VectorGenUniform([-dL,0,0],e_size);
dr_1p2030 = func_m_VectorGenUniform([ dL,0,0],e_size);

dr_102m30 = func_m_VectorGenUniform([0,-dL,0],e_size);
dr_102p30 = func_m_VectorGenUniform([0, dL,0],e_size);

dr_10203m = func_m_VectorGenUniform([0,0,-dL],e_size);
dr_10203p = func_m_VectorGenUniform([0,0, dL],e_size);

% Coord of 6 points surounding R0
R_1m2030 = func_m_VectorSum(R0,dr_1m2030);
R_1p2030 = func_m_VectorSum(R0,dr_1p2030);
R_102m30 = func_m_VectorSum(R0,dr_102m30);
R_102p30 = func_m_VectorSum(R0,dr_102p30);
R_10203m = func_m_VectorSum(R0,dr_10203m);
R_10203p = func_m_VectorSum(R0,dr_10203p);

% Interpolate value of data_thres at 6 points
DT_1m2030 = interp3(X1,X2,X3,NparaAll,R_1m2030.e1,R_1m2030.e2,R_1m2030.e3,'linear');
DT_1p2030 = interp3(X1,X2,X3,NparaAll,R_1p2030.e1,R_1p2030.e2,R_1p2030.e3,'linear');
DT_102m30 = interp3(X1,X2,X3,NparaAll,R_102m30.e1,R_102m30.e2,R_102m30.e3,'linear');
DT_102p30 = interp3(X1,X2,X3,NparaAll,R_102p30.e1,R_102p30.e2,R_102p30.e3,'linear');
DT_10203m = interp3(X1,X2,X3,NparaAll,R_10203m.e1,R_10203m.e2,R_10203m.e3,'linear');
DT_10203p = interp3(X1,X2,X3,NparaAll,R_10203p.e1,R_10203p.e2,R_10203p.e3,'linear');

IsPositiveSign = DT_1m2030>=0 & DT_1p2030>=0 & DT_102m30>=0 & DT_102p30>=0 & DT_10203m>=0 & DT_10203p>=0 & Npara>=0;
IsNegativeSign = DT_1m2030<=0 & DT_1p2030<=0 & DT_102m30<=0 & DT_102p30<=0 & DT_10203m<=0 & DT_10203p<=0 & Npara<=0;
IsSameSign = IsPositiveSign | IsNegativeSign;

Is3DExtrema = Npara>DT_1m2030 & Npara>DT_1p2030 & ...
              Npara>DT_102m30 & Npara>DT_102p30 & ...
              Npara>DT_10203m & Npara>DT_10203p & ...
              IsSameSign;
end