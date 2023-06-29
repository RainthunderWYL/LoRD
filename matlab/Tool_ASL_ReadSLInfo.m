function out = Tool_ASL_ReadSLInfo(SLInfo,Name,varargin)

nvarargin = length(varargin);
if nvarargin == 1
    ExtraData = varargin{1};
    if ExtraData == 0
        out = SLInfo.Data(:,func_ASL_SLInfoStruc(Name));
    else
        out = SLInfo.ExtraData(:,func_ASL_SLInfoStruc(Name,1));
    end
elseif nvarargin == 0
    out = SLInfo.Data(:,func_ASL_SLInfoStruc(Name));
else
    error('Tool_ASL_ReadSLInfo: Wrong number of inputs');
end
end