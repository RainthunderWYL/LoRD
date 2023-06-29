function out = Tool_ASL_PeekLocalProf(haxis,SLInfo,B1,B2,B3,x1,x2,x3,...
    MSP_ScaleRatio,NumRefine,MSPDispScale,varargin)
% Tool_ASL_PeekLocalProf(gca,SPInfo,B1,B2,B3), ginput3D
% Tool_ASL_PeekLocalProf(gca,SPInfo,B1,B2,B3,[x1,x2,x3])
% Tool_ASL_PeekLocalProf(gca,SPInfo,B1,B2,B3,idx)

if nargin == 12
    if length(varargin{1}) > 1
        Rin = varargin{1};
    else
        Rin = SLInfo.Data(varargin{1},1:3);
    end
else
    Rin = func_ginput3D(1);
end
IdxAll = 1:size(SLInfo.Data,1);
b = SLInfo.Data(:,1) == Rin(1) & SLInfo.Data(:,2) == Rin(2) & SLInfo.Data(:,3) == Rin(3);

if sum(b)==0
    error('ERROR: ASP_PeekLocalProf.m: Input position is NOT in SPInfo!');
else
    idx = IdxAll(b);
end

r0i = double(SLInfo.Data(idx,1:3));
e1i = double(SLInfo.ExtraData(idx,5:7));
e2i = double(SLInfo.ExtraData(idx,8:10));
e3i = double(SLInfo.ExtraData(idx,11:13));

out.h = func_ProjectFieldProf(haxis,double(B1),double(B2),double(B3),double(x1),double(x2),double(x3),r0i,e1i,e2i,e3i,...
    MSP_ScaleRatio,NumRefine,MSPDispScale);
out.R0 = Rin;

end