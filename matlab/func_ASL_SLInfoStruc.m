function out = func_ASL_SLInfoStruc(varargin)
% out = func_ASL_SLInfoStruc(): SLInfoStruc.Data, SLInfoStruc.ExtraData
% idx = func_ASL_SLInfoStruc('DataName')
% idx = func_ASL_SLInfoStruc('DataName',OutExtraData)

%%%%%%%%%%
SLInfoStruc.Data.fields{1} = 'x1';
SLInfoStruc.Data.fields{2} = 'x2';
SLInfoStruc.Data.fields{3} = 'x3';
SLInfoStruc.Data.fields{4} = 'SLType';
SLInfoStruc.Data.fields{5} = 'Is2DExtrema';
SLInfoStruc.Data.fields{6} = 'Is3DExtrema';
SLInfoStruc.Data.fields{7} = 'EigAngle';
SLInfoStruc.Data.fields{8} = 'RatioMTrace';


%%%%%%%%%%
SLInfoStruc.ExtraData.fields{1} = 'DB11';
SLInfoStruc.ExtraData.fields{2} = 'DB12';
SLInfoStruc.ExtraData.fields{3} = 'DB21';
SLInfoStruc.ExtraData.fields{4} = 'DB22';
SLInfoStruc.ExtraData.fields{5}  = 'e11';
SLInfoStruc.ExtraData.fields{6}  = 'e12';
SLInfoStruc.ExtraData.fields{7}  = 'e13';
SLInfoStruc.ExtraData.fields{8}  = 'e21';
SLInfoStruc.ExtraData.fields{9}  = 'e22';
SLInfoStruc.ExtraData.fields{10} = 'e23';
SLInfoStruc.ExtraData.fields{11} = 'e31';
SLInfoStruc.ExtraData.fields{12} = 'e32';
SLInfoStruc.ExtraData.fields{13} = 'e33';

SLInfoStruc.Data.n = length(SLInfoStruc.Data.fields);
SLInfoStruc.ExtraData.n = length(SLInfoStruc.ExtraData.fields);

StrAllDataField = SLInfoStruc.Data.fields{1};
for i = 2:length(SLInfoStruc.Data.fields)
    StrAllDataField = strcat(StrAllDataField,sprintf(', %s',SLInfoStruc.Data.fields{i}));
end

StrAllExtraDataField = SLInfoStruc.ExtraData.fields{1};
for i = 2:length(SLInfoStruc.ExtraData.fields)
    StrAllExtraDataField = strcat(StrAllExtraDataField,sprintf(', %s',SLInfoStruc.ExtraData.fields{i}));
end

if nargin == 0
    out = SLInfoStruc;
elseif nargin == 1
    Data_Str = varargin{1};
    IdxAll = 1:length(SLInfoStruc.Data.fields);
    b = strcmp(SLInfoStruc.Data.fields,Data_Str);
    out = IdxAll(b);
    if isempty(out)
        error('func_ASL_SLInfoStruc: Available Data Names are:\n%s',StrAllDataField);
    end
elseif nargin == 2
    Data_Str = varargin{1};
    OutExtraData = varargin{2};

    if OutExtraData == 0
        IdxAll = 1:length(SLInfoStruc.Data.fields);
        b = strcmp(SLInfoStruc.Data.fields,Data_Str);
        out = IdxAll(b);
        if isempty(out)
            error('func_ASL_SLInfoStruc: Available Data Names are:\n%s',StrAllDataField);
        end
    else
        IdxAll = 1:length(SLInfoStruc.ExtraData.fields);
        b = strcmp(SLInfoStruc.ExtraData.fields,Data_Str);
        out = IdxAll(b);
        if isempty(out)
            error('func_ASL_SLInfoStruc: Available ExtraData Names are:\n%s',StrAllExtraDataField);
        end
    end
else
    error('func_ASL_SLInfoStruc: input parameters should <= 2');
end

end