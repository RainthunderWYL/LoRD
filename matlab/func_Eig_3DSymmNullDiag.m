function [L1,L2,L3,funcV1,funcV2,funcV3] = func_Eig_3DSymmNullDiag(a,b,c)
% M = |0,a,b|
%     |a,0,c|
%     |b,c,0|

% Resize inputs into 1-D array
OriginalSize = size(a);
len = numel(a);
a = reshape(a,[1,len]);
b = reshape(b,[1,len]);
c = reshape(c,[1,len]);

Ltmp = zeros(3,len);

a2pb2pc2 = a.^2+b.^2+c.^2;
CSqrt = (9*a.*b.*c+1/3*sqrt(729.*(a.*b.*c).^2-27.*a2pb2pc2.^3)).^(1/3);

Ltmp(1,:) = real((3^(1/3)*a2pb2pc2+CSqrt.^2)./(3^(2/3)*CSqrt));
Ltmp(2,:) = real((-6i*2^(1/3)*(-1i+sqrt(3))*a2pb2pc2+2^(2/3)*(-1+sqrt(3)*1i)*6^(2/3)*CSqrt.^2)./...
    (12*6^(1/3)*CSqrt));
Ltmp(3,:) = real((1i*(1i+sqrt(3))*a2pb2pc2)./(2^(2/3)*6^(1/3)*CSqrt)-...
    ((1+1i*sqrt(3))*6^(1/3)*CSqrt)/(6*2^(1/3)));

Ltmp = sort(Ltmp);

L1 = reshape(Ltmp(1,:),OriginalSize);
L2 = reshape(Ltmp(2,:),OriginalSize);
L3 = reshape(Ltmp(3,:),OriginalSize);

funcV1 = @(L,A,B,C) (L.*B+A.*C)./(L.^2-A.*A);
funcV2 = @(L,A,B,C) (L.*C+A.*B)./(L.^2-A.*A);
funcV3 = @(L,A,B,C) ones(size(A));

% % V1
% V1{1} = real((l1.*b+a.*c)./(l1.^2-a.*a));
% V1{2} = real((l1.*c+a.*b)./(l1.^2-a.*a));
% V1{3} = ones(size(a));
% 
% %V2
% V2{1} = real((l2.*b + a.*c)./(l2.^2-a.*a));
% V2{2} = real((l2.*c + a.*b)./(l2.^2-a.*a));
% V2{3} = ones(size(a));
% 
% %V3
% V3{1} = real((l3.*b + a.*c)./(l3.^2-a.*a));
% V3{2} = real((l3.*c + a.*b)./(l3.^2-a.*a));
% V3{3} = ones(size(a));
% 
% VNorm = sqrt(V1{1}.^2 + V1{2}.^2 +V1{3}.^2);
% V1{1} = V1{1}./VNorm;
% V1{2} = V1{2}./VNorm;
% V1{3} = V1{3}./VNorm;
% 
% VNorm = sqrt(V2{1}.^2 + V2{2}.^2 +V2{3}.^2);
% V2{1} = V2{1}./VNorm;
% V2{2} = V2{2}./VNorm;
% V2{3} = V2{3}./VNorm;
% 
% VNorm = sqrt(V3{1}.^2 + V3{2}.^2 +V3{3}.^2);
% V3{1} = V3{1}./VNorm;
% V3{2} = V3{2}./VNorm;
% V3{3} = V3{3}./VNorm;
end