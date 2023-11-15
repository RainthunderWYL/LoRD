function Is2DExtrema = func_ARD_Is2DExtrema(...
                                R0,Tinv,Npara, x1,x2,x3,NparaAll)

dL = mean([mean(diff(x1)),mean(diff(x3)),mean(diff(x3))]);
[X1,X2,X3] = meshgrid(x1,x2,x3);

% Local coord of 8 points surounding R0
e_size = size(R0.e1);
r_1m20 = func_m_VectorGenUniform([-dL,0,0],e_size);
r_1p20 = func_m_VectorGenUniform([dL,0,0],e_size);

r_102m = func_m_VectorGenUniform([0,-dL,0],e_size);
r_102p = func_m_VectorGenUniform([0,dL,0],e_size);

% r_1m2m = func_m_VectorGenUniform([-dL,-dL,0],e_size);
% r_1p2p = func_m_VectorGenUniform([dL,dL,0],e_size);
% 
% r_1m2p = func_m_VectorGenUniform([-dL,dL,0],e_size);
% r_1p2m = func_m_VectorGenUniform([dL,-dL,0],e_size);

% Tranform back to x-y-z frame
R_1m20 = func_m_VectorSum(R0,func_m_VectorTransform(Tinv,r_1m20));
R_1p20 = func_m_VectorSum(R0,func_m_VectorTransform(Tinv,r_1p20));
R_102m = func_m_VectorSum(R0,func_m_VectorTransform(Tinv,r_102m));
R_102p = func_m_VectorSum(R0,func_m_VectorTransform(Tinv,r_102p));
% R_1m2m = func_m_VectorSum(R0,func_m_VectorTransform(Tinv,r_1m2m));
% R_1p2p = func_m_VectorSum(R0,func_m_VectorTransform(Tinv,r_1p2p));
% R_1m2p = func_m_VectorSum(R0,func_m_VectorTransform(Tinv,r_1m2p));
% R_1p2m = func_m_VectorSum(R0,func_m_VectorTransform(Tinv,r_1p2m));

% Interpolate value of data_thres at 8 points
DT_1m20 = interp3(X1,X2,X3,NparaAll,R_1m20.e1,R_1m20.e2,R_1m20.e3,'linear');
DT_1p20 = interp3(X1,X2,X3,NparaAll,R_1p20.e1,R_1p20.e2,R_1p20.e3,'linear');
DT_102m = interp3(X1,X2,X3,NparaAll,R_102m.e1,R_102m.e2,R_102m.e3,'linear');
DT_102p = interp3(X1,X2,X3,NparaAll,R_102p.e1,R_102p.e2,R_102p.e3,'linear');
% DT_1m2m = interp3(X1,X2,X3,NparaAll,R_1m2m.e1,R_1m2m.e2,R_1m2m.e3,'linear');
% DT_1p2p = interp3(X1,X2,X3,NparaAll,R_1p2p.e1,R_1p2p.e2,R_1p2p.e3,'linear');
% DT_1m2p = interp3(X1,X2,X3,NparaAll,R_1m2p.e1,R_1m2p.e2,R_1m2p.e3,'linear');
% DT_1p2m = interp3(X1,X2,X3,NparaAll,R_1p2m.e1,R_1p2m.e2,R_1p2m.e3,'linear');

% IsProjectLocalMaxima = DataThres>DT_1m20 & DataThres>DT_1p20 & ...
%     DataThres>DT_102m & DataThres>DT_102p & ...
%     DataThres>DT_1m2m & DataThres>DT_1p2p & ...
%     DataThres>DT_1m2p & DataThres>DT_1p2m;
IsPositiveSign = DT_1m20>=0 & DT_1p20>=0 & DT_102m>=0 & DT_102p>=0 & Npara>=0;
IsNegativeSign = DT_1m20<=0 & DT_1p20<=0 & DT_102m<=0 & DT_102p<=0 & Npara<=0;
IsSameSign = IsPositiveSign | IsNegativeSign;

Is2DExtrema = abs(Npara) > abs(DT_1m20) & ...
              abs(Npara) > abs(DT_1p20) & ...
              abs(Npara) > abs(DT_102m) & ...
              abs(Npara) > abs(DT_102p) & ...
              IsSameSign;
end