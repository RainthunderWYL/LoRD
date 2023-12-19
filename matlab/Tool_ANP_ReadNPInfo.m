function out = Tool_ANP_ReadNPInfo(NPInfo,Name,varargin)

nvarargin = length(varargin);
if nvarargin == 1
    ExtraData = varargin{1};
    if ExtraData == 0
        out = NPInfo.Data(:,func_ANP_NPInfoStruc(Name));
    else
        out = NPInfo.ExtraData(:,func_ANP_NPInfoStruc(Name,1));
    end
elseif nvarargin == 0
    out = NPInfo.Data(:,func_ANP_NPInfoStruc(Name));
else
    error('Tool_ANP_ReadNPInfo: Wrong number of inputs');
end
end