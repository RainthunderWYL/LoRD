function [L1,L2,L3,funcV1,funcV2,funcV3] = func_Eig_3DTraceLess(M)
% M = |a,c,e|
%     |d,b,g|
%     |f,h,-a-b|
% M = [a;d;f;c;b;h;e;g;-a-b]

a = M(1,:);
d = M(2,:);
f = M(3,:);
c = M(4,:);
b = M(5,:);
h = M(6,:);
e = M(7,:);
g = M(8,:);
% x = M(9,:);

Q = -(a.*a+b.*b+a.*b+c.*d+e.*f+g.*h);
R = (a+b).*(a.*b-c.*d)+e.*(b.*f-d.*h)+g.*(a.*h-c.*f);
Delta = -(4*Q.^3+27*R.^2);
Cons = (1/3*sqrt(-27*Delta)-9*R).^(1/3);

Ltmp = zeros(3,length(a));
Ltmp(1,:) = (2.^(1/3)*Cons.^2-2*3^(1/3)*Q)./(6^(2/3)*Cons);
Ltmp(2,:) = 1i*(-2*3^(1/3)*(-1i+sqrt(3))*(-Q)+2^(1/3)*(1i+sqrt(3))*Cons.^2)./(2*6^(2/3)*Cons);
Ltmp(3,:) = 1i*(2*3^(1/3)*(1i+sqrt(3))*(-Q)-2^(1/3)*(-1i+sqrt(3))*Cons.^2)./(2*6^(2/3)*Cons);

bR = func_iszero(imag(Ltmp(1,:))) & func_iszero(imag(Ltmp(2,:))) & func_iszero(imag(Ltmp(3,:)));
L_R = real(Ltmp(:,bR));
L_C = Ltmp(:,~bR);

if ~isempty(L_R)
    [~,IR] = sort(abs(L_R));
    Ltmp(:,bR) = func_MapSortedIndex(L_R,IR);
end

if ~isempty(L_C)
    [~,IC] = sort(abs(imag(L_C)));
    Ltmp(:,~bR) = func_MapSortedIndex(L_C,IC);
end

L1 = real(Ltmp(1,:));
L2 = Ltmp(2,:);
L3 = Ltmp(3,:);

funcV1 = @(L,A,B,C,D,E,F,G,H) (C.*G-E.*(B-L))./((A-L).*(B-L)-C.*D);
funcV2 = @(L,A,B,C,D,E,F,G,H) (D.*E-G.*(A-L))./((A-L).*(B-L)-C.*D);
funcV3 = @(L,A,B,C,D,E,F,G,H) ones(size(A));
end
