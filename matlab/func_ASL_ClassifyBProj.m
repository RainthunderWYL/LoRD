function out = func_ASL_ClassifyBProj(A,B,C,D)
% A = DB11; B = DB22; C = DB12; D = DB21;

SLType = zeros(size(A));

Tr = A + B;
[lambda1,lambda2,funcV1,funcV2] = func_Eig_2D(A,B,C,D);
Jthres = sqrt((A-B).^2+(C+D).^2);
J3 = D - C;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SLType
b_is3D = abs(Tr) > 0;
b_is2D = ~b_is3D;
b_isRealEig = abs(J3) < Jthres;
b_isImagEig = abs(J3) > Jthres;
b_isParaLine = abs(J3) == Jthres;

SLType(b_is2D & b_isParaLine) = 9; % 2D Anti-parallel
SLType(b_is2D & b_isRealEig) = 7; % 2D X
SLType(b_is2D & b_isImagEig) = 8; % 2D O

SLType(b_is3D & b_isParaLine) = 6; % 3D Anti-parallel
SLType(b_is3D & b_isRealEig & lambda1.*lambda2 < 0) = 1; % 3D X
SLType(b_is3D & b_isRealEig & lambda1>0 & lambda2>0) = 4; % 3D Repelling
SLType(b_is3D & b_isRealEig & lambda1<0 & lambda2<0) = 5; % 3D Attracting


SLType(b_is3D & b_isImagEig & Tr>0) = 2; % 3D O (Repelling)
SLType(b_is3D & b_isImagEig & Tr<0) = 3; % 3D O (Attracting)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EigAngle: 0~90
% Calculate EigAngle
EigAngle = zeros(size(SLType));
EigAngle(b_isImagEig) = atan(sqrt((abs(J3(b_isImagEig))./Jthres(b_isImagEig)).^2-1))/pi*180;
EigAngle(~b_isImagEig) = atan(sqrt((Jthres(~b_isImagEig)./abs(J3(~b_isImagEig))).^2-1))/pi*180;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Importance of trace part
Mprime_norm = sqrt((0.5*(A-B)).^2*2+C.^2+D.^2);
T_norm = sqrt((0.5*(A+B)).^2*2);


%%%% Directly calculate angle spanned by eigen-vectors: Not Used 
%%%% For real eigenvectors, their angle is the same as EigAngle;
% infunc_Ve1 = @(l) funcV1(l,A,B,C,D);
% infunc_Ve2 = @(l) funcV2(l,A,B,C,D);;
% V1_tmp.e1 = infunc_Ve1(lambda1);
% V1_tmp.e2 = infunc_Ve2(lambda1);
% V2_tmp.e1 = infunc_Ve1(lambda2);
% V2_tmp.e2 = infunc_Ve2(lambda2);
% V1.e1 = V1_tmp.e1;
% V1.e2 = V1_tmp.e2;
% V2.e1 = V2_tmp.e1;
% V2.e2 = V2_tmp.e2;
% V1.e1(b_isImagEig) = real(V1_tmp.e1(b_isImagEig));
% V1.e2(b_isImagEig) = real(V1_tmp.e2(b_isImagEig));
% V2.e1(b_isImagEig) = imag(V1_tmp.e1(b_isImagEig));
% V2.e2(b_isImagEig) = imag(V1_tmp.e2(b_isImagEig));
% % 
% EigAngleM = acos(func_m_VectorDot(V1,V2))/pi*180;% degree
% b_lt90 = EigAngleM>90;
% EigAngleM(b_lt90) = 180-EigAngleM(b_lt90);

% Output
out.SLType = SLType;
out.EigAngle = EigAngle;
out.RatioMTrace = T_norm./Mprime_norm;% Larger RatioTM means trace effect is more dominate
% out.EigAngleM = EigAngleM;
end