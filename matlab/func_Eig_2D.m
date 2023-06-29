function [L1,L2,funcV1,funcV2] = func_Eig_2D(A,B,C,D)
% M = |A,C|
%     |D,B|

Tr = A + B;
Det = A .* B - C .* D;
Disc = Tr.^2 - 4.*Det;

L1 = 0.5*(Tr + sqrt(Disc));
L2 = 0.5*(Tr - sqrt(Disc));

funcV1 = @(l,a,b,c,d) (l - b)./sqrt(d.^2+abs(l-b).^2);
funcV2 = @(l,a,b,c,d) d./sqrt(d.^2+abs(l-b).^2);
end