function [x,y,z] = func_GenRandCorInSphere(X0,R,N)
% https://ww2.mathworks.cn/help/matlab/math/numbers-placed-randomly-within-volume-of-sphere.html
rng(0);
rvals = 2*rand(N,1)-1;
elevation = asin(rvals);
azimuth = 2*pi*rand(N,1);
radii = R*(rand(N,1).^(1/3));
[x,y,z] = sph2cart(azimuth,elevation,radii);

x = x+X0(1);
y = y+X0(2);
z = z+X0(3);
end