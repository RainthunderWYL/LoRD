function out = func_ARD_RDInfoStruc(varargin)
% out = func_ARD_RDInfoStruc(): RDInfoStruc.Data, RDInfoStruc.ExtraData
% idx = func_ARD_RDInfoStruc('DataName')
% idx = func_ARD_RDInfoStruc('DataName',OutExtraData)

%%%%%%%%%%
RDInfoStruc.Data.fields{1} = 'x1';
RDInfoStruc.Data.fields{2} = 'x2';
RDInfoStruc.Data.fields{3} = 'x3';
RDInfoStruc.Data.fields{4} = 'RDType';
RDInfoStruc.Data.fields{5} = 'Is2DExtrema';
RDInfoStruc.Data.fields{6} = 'EigAngle';
RDInfoStruc.Data.fields{7} = 'RatioMTrace';
% RDInfoStruc.Data.fields{8} = 'Is3DExtrema';

%%%%%%%%%%
RDInfoStruc.ExtraData.fields{1} = 'B0';
RDInfoStruc.ExtraData.fields{2} = 'DB11';
RDInfoStruc.ExtraData.fields{3} = 'DB12';
RDInfoStruc.ExtraData.fields{4} = 'DB21';
RDInfoStruc.ExtraData.fields{5} = 'DB22';
RDInfoStruc.ExtraData.fields{6} = 'DB31';
RDInfoStruc.ExtraData.fields{7} = 'DB32';
RDInfoStruc.ExtraData.fields{8} = 'DB13';
RDInfoStruc.ExtraData.fields{9} = 'DB23';
RDInfoStruc.ExtraData.fields{10} = 'DB33';
RDInfoStruc.ExtraData.fields{11} = 'e11';
RDInfoStruc.ExtraData.fields{12} = 'e12';
RDInfoStruc.ExtraData.fields{13} = 'e13';
RDInfoStruc.ExtraData.fields{14} = 'e21';
RDInfoStruc.ExtraData.fields{15} = 'e22';
RDInfoStruc.ExtraData.fields{16} = 'e23';
RDInfoStruc.ExtraData.fields{17} = 'e31';
RDInfoStruc.ExtraData.fields{18} = 'e32';
RDInfoStruc.ExtraData.fields{19} = 'e33';
RDInfoStruc.ExtraData.fields{20} = 'DNpara1';
RDInfoStruc.ExtraData.fields{21} = 'DNpara2';
RDInfoStruc.ExtraData.fields{22} = 'Gamma1';
RDInfoStruc.ExtraData.fields{23} = 'Gamma2';
RDInfoStruc.ExtraData.fields{24} = 'CurlN';
RDInfoStruc.ExtraData.fields{25} = 'BxCurlN';

RDInfoStruc.Data.n = length(RDInfoStruc.Data.fields);
RDInfoStruc.ExtraData.n = length(RDInfoStruc.ExtraData.fields);

StrAllDataField = RDInfoStruc.Data.fields{1};
for i = 2:length(RDInfoStruc.Data.fields)
    StrAllDataField = strcat(StrAllDataField,sprintf(', %s',RDInfoStruc.Data.fields{i}));
end

StrAllExtraDataField = RDInfoStruc.ExtraData.fields{1};
for i = 2:length(RDInfoStruc.ExtraData.fields)
    StrAllExtraDataField = strcat(StrAllExtraDataField,sprintf(', %s',RDInfoStruc.ExtraData.fields{i}));
end

if nargin == 0
    out = RDInfoStruc;
elseif nargin == 1
    Data_Str = varargin{1};
    IdxAll = 1:length(RDInfoStruc.Data.fields);
    b = strcmp(RDInfoStruc.Data.fields,Data_Str);
    out = IdxAll(b);
    if isempty(out)
        error('func_ARD_RDInfoStruc: Available Data Names are:\n%s',StrAllDataField);
    end
elseif nargin == 2
    Data_Str = varargin{1};
    OutExtraData = varargin{2};

    if OutExtraData == 0
        IdxAll = 1:length(RDInfoStruc.Data.fields);
        b = strcmp(RDInfoStruc.Data.fields,Data_Str);
        out = IdxAll(b);
        if isempty(out)
            error('func_ARD_RDInfoStruc: Available Data Names are:\n%s',StrAllDataField);
        end
    else
        IdxAll = 1:length(RDInfoStruc.ExtraData.fields);
        b = strcmp(RDInfoStruc.ExtraData.fields,Data_Str);
        out = IdxAll(b);
        if isempty(out)
            error('func_ARD_RDInfoStruc: Available ExtraData Names are:\n%s',StrAllExtraDataField);
        end
    end
else
    error('func_ARD_RDInfoStruc: input parameters should <= 2');
end

end