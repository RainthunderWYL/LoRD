function out = Tool_TraceFieldLineByClick(B1,B2,B3,x1,x2,x3,Radius,Dpi,varargin)
% Click a point and seed initial positions in a sphere.

% LineConf.Len
% LineConf.Color
% LineConf.Style
% LineConf.Width

R0 = func_ginput3D(1);

if isempty(varargin)
    hfl = Tool_PreviewFieldLines(B1,B2,B3,x1,x2,x3,'s',R0,Radius,Dpi);
else
    hfl = Tool_PreviewFieldLines(B1,B2,B3,x1,x2,x3,'s',R0,Radius,Dpi,varargin);
end
out.hfl = hfl;
out.R0 = R0;
end
