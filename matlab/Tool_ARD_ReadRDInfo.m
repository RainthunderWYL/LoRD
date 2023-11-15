function out = Tool_ARD_ReadRDInfo(RDInfo,Name,varargin)

nvarargin = length(varargin);
if nvarargin == 1
    ExtraData = varargin{1};
    if ExtraData == 0
        out = RDInfo.Data(:,func_ARD_RDInfoStruc(Name));
    else
        out = RDInfo.ExtraData(:,func_ARD_RDInfoStruc(Name,1));
    end
elseif nvarargin == 0
    out = RDInfo.Data(:,func_ARD_RDInfoStruc(Name));
else
    error('Tool_ARD_ReadRDInfo: Wrong number of inputs');
end
end