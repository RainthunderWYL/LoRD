function [L1,L2,L3,funcV1,funcV2,funcV3] = func_Eig_3D(M)
% Not finish!!!!!
% M = |a,c,e|
%     |d,b,g|
%     |f,h,x|
% M = [a;d;f;c;b;h;e;g;x]

a = M(1,:);
d = M(2,:);
f = M(3,:);
c = M(4,:);
b = M(5,:);
h = M(6,:);
e = M(7,:);
g = M(8,:);
x = M(9,:);



funcV1 = @(L,A,B,C,X,Y,Z) (A.*C-B.*(Y-L))./((X-L).*(Y-L)-A.*A);
funcV2 = @(L,A,B,C,X,Y,Z) (A.*B-C.*(X-L))./((X-L).*(Y-L)-A.*A);
funcV3 = @(L,A,B,C,X,Y,Z) ones(size(A));
end