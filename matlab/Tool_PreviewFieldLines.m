function out = Tool_PreviewFieldLines(B1,B2,B3,x1,x2,x3,Type,R0,R1,Dpi,varargin)
% Type: 'l': R0 and R1 are start and end points of a line; Dpi is grid
% number
% Type: 's': R0 is center, R1 is radius, Dpi is sample number
% LineConf.Len
% LineConf.Color
% LineConf.Style
% LineConf.Width

if isempty(varargin)
    LineConf.Len = 10000;
    LineConf.Color = 'k';
    LineConf.Style = '-';
    LineConf.Width = 0.5;
else
    LineConf = varargin{1};
end

% Generate start points
if strcmp(Type,'l')
    R1 = reshape(R1,[1,length(R1)]);
    R0 = reshape(R0,[1,length(R0)]);
    r = R1-R0;
    er = r/norm(r);
    len = linspace(0,norm(r),Dpi);
    SPs = repmat(R0,[Dpi,1]);
    ER = repmat(er,[Dpi,1]);
    LEN = repmat(len',[1,3]);

    SPs = SPs + LEN.*ER;
elseif strcmp(Type,'s')
    xdomain = [R0(1)-R1,R0(1)+R1];
    ydomain = [R0(2)-R1,R0(2)+R1];
    zdomain = [R0(3)-R1,R0(3)+R1];
    SPs = zeros(Dpi,3);

	[SPs(:,1),SPs(:,2),SPs(:,3)] = infunc_GenRandCorInSphere(R0,R1,Dpi);
else
    error('In function ASP_PreviewFieldLines: Wrong Type: l or s!');
end
[X,Y,Z] = meshgrid(x1,x2,x3);
out=infunc_DrawSymmStreamLine3D( X,Y,Z,B1,B2,B3,...
    SPs(:,1),SPs(:,2),SPs(:,3),[.1,LineConf.Len],...
    LineConf.Color,LineConf.Style,LineConf.Width);
end

function output_args = infunc_DrawSymmStreamLine3D( X,Y,Z,B1,B2,B3,...
    startx,starty,startz,option,color,style,width,varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if isempty(varargin)
    hfl1 = streamline(X,Y,Z,B1,B2,B3,startx,starty,startz,option);
    set(hfl1,'Color',color,'LineStyle',style,'LineWidth',width);

    hfl2 = streamline(X,Y,Z,-B1,-B2,-B3,startx,starty,startz,option);
    set(hfl2,'Color',color,'LineStyle',style,'LineWidth',width);
else
    axes1 = varargin{1};
    hfl1 = streamline(axes1,X,Y,Z,B1,B2,B3,startx,starty,startz,option);
    set(hfl1,'Color',color,'LineStyle',style,'LineWidth',width);

    hfl2 = streamline(axes1,X,Y,Z,-B1,-B2,-B3,startx,starty,startz,option);
    set(hfl2,'Color',color,'LineStyle',style,'LineWidth',width);
end
output_args{1} = hfl1;
output_args{2} = hfl2;
end

function [x,y,z] = infunc_GenRandCorInSphere(X0,R,N)
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
