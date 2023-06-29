function [e1,e2,e3] = func_VectorLocalFrame(V,varargin)
%V.e1 = e31, V.e2 = e32, V.e3 = e33
%varargin = 1: skip normalize V

% e3 is b
if nargin>1
    e3 = V;
else
    e3 = func_m_VectorDirection(V);
end
% Generate e2: e2 = [e3.e2,-e3.e1,0] which is perpendicular to e3
e2.e1 = e3.e2;
e2.e2 = -e3.e1;
e2.e3 = zeros(size(e2.e1));
% Address special case of e3 = ez
b_e2 = e3.e1 == 0 & e3.e2 == 0;
e2.e1(b_e2) = 0;
e2.e2(b_e2) = 1;
e2.e3(b_e2) = 0;

e2norm = sqrt(e2.e1.^2+e2.e2.^2+e2.e3.^2);
e2.e1 = e2.e1./e2norm;
e2.e2 = e2.e2./e2norm;
e2.e3 = e2.e3./e2norm;

% e2 = func_m_VectorDirection(e2);
e1 = func_m_VectorCross(e2,e3);
end