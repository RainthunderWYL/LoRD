function Delta = func_Discriminant4TracelessMatrix(varargin)
% Delta = func_Discriminant4TracelessMatrix(M): M is a 3x3 matrix;
% Delta = func_Discriminant4TracelessMatrix(a,b,c,d,e,f,g,h):
% M = [a c e
%      d b g
%      f h -a=b]

if nargin == 1
    M = varargin{1};
    a = M(1,1);
    b = M(2,2);
    c = M(1,2);
    d = M(2,1);
    e = M(1,3);
    f = M(3,1);
    g = M(2,3);
    h = M(3,2);
elseif nargin >1
    a = varargin{1};
    b = varargin{2};
    c = varargin{3};
    d = varargin{4};
    e = varargin{5};
    f = varargin{6};
    g = varargin{7};
    h = varargin{8};
else
    error('func_Discriminant4TracelessMatrix: No input!');
end

Q = -(a.^2+b.^2+a.*b+c.*d+e.*f+g.*h);
R = (a+b).*(a.*b-c.*d)+e.*(b.*f-d.*h)+g.*(a.*h-c.*f);
Delta = -(4*Q.^3+27*R.^2);
end