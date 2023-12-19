function out = func_ANP_NPInfoStruc(varargin)
% out = func_ANP_NPInfoStruc(): NPInfoStruc.Data, NPInfoStruc.ExtraData
% idx = func_ANP_NPInfoStruc('DataName')
% idx = func_ANP_NPInfoStruc('DataName',OutExtraData)

%%%%%%%%%%
NPInfoStruc.Data.fields{1} = 'x1';
NPInfoStruc.Data.fields{2} = 'x2';
NPInfoStruc.Data.fields{3} = 'x3';
NPInfoStruc.Data.fields{4} = 'NPType';

%%%%%%%%%%
NPInfoStruc.ExtraData.fields{1} = 'DB11';
NPInfoStruc.ExtraData.fields{2} = 'DB21';
NPInfoStruc.ExtraData.fields{3} = 'DB31';
NPInfoStruc.ExtraData.fields{4} = 'DB12';
NPInfoStruc.ExtraData.fields{5} = 'DB22';
NPInfoStruc.ExtraData.fields{6} = 'DB32';
NPInfoStruc.ExtraData.fields{7} = 'DB13';
NPInfoStruc.ExtraData.fields{8} = 'DB23';
NPInfoStruc.ExtraData.fields{9} = 'DB33';
NPInfoStruc.ExtraData.fields{10} = 'DetM';
NPInfoStruc.ExtraData.fields{11} = 'J';
NPInfoStruc.ExtraData.fields{12} = 'Jperp';
NPInfoStruc.ExtraData.fields{13} = 'DM';
NPInfoStruc.ExtraData.fields{14} = 'TrM';

NPInfoStruc.Data.n = length(NPInfoStruc.Data.fields);
NPInfoStruc.ExtraData.n = length(NPInfoStruc.ExtraData.fields);

StrAllDataField = NPInfoStruc.Data.fields{1};
for i = 2:length(NPInfoStruc.Data.fields)
    StrAllDataField = strcat(StrAllDataField,sprintf(', %s',NPInfoStruc.Data.fields{i}));
end

StrAllExtraDataField = NPInfoStruc.ExtraData.fields{1};
for i = 2:length(NPInfoStruc.ExtraData.fields)
    StrAllExtraDataField = strcat(StrAllExtraDataField,sprintf(', %s',NPInfoStruc.ExtraData.fields{i}));
end

if nargin == 0
    out = NPInfoStruc;
elseif nargin == 1
    Data_Str = varargin{1};
    IdxAll = 1:length(NPInfoStruc.Data.fields);
    b = strcmp(NPInfoStruc.Data.fields,Data_Str);
    out = IdxAll(b);
    if isempty(out)
        error('func_ANP_NPInfoStruc: Available Data Names are:\n%s',StrAllDataField);
    end
elseif nargin == 2
    Data_Str = varargin{1};
    OutExtraData = varargin{2};

    if OutExtraData == 0
        IdxAll = 1:length(NPInfoStruc.Data.fields);
        b = strcmp(NPInfoStruc.Data.fields,Data_Str);
        out = IdxAll(b);
        if isempty(out)
            error('func_ANP_NPInfoStruc: Available Data Names are:\n%s',StrAllDataField);
        end
    else
        IdxAll = 1:length(NPInfoStruc.ExtraData.fields);
        b = strcmp(NPInfoStruc.ExtraData.fields,Data_Str);
        out = IdxAll(b);
        if isempty(out)
            error('func_ANP_NPInfoStruc: Available ExtraData Names are:\n%s',StrAllExtraDataField);
        end
    end
else
    error('func_ANP_NPInfoStruc: input parameters should <= 2');
end

end