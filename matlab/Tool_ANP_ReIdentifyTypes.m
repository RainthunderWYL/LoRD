function out = Tool_ANP_ReIdentifyTypes(FileFullPathTemplate,IdxMinMax,Stride)
% Need Test!!!
num_extradata = 6;
for idx = IdxMinMax(1):Stride:IdxMinMax(2)
    filename = sprintf(FileFullPathTemplate,idx);

    data = load(filename);
    
    Parameters = data.Parameters;
    NPInfo = data.NPInfo;
    NPInfo_Raw = NPInfo(:,1:end-num_extradata);

    NPInfo = func_ANP_IdentifyNullPointType(NPInfo_Raw,Parameters.ANP_FixTrace);

    data.NPInfo = NPInfo;
    save(filename,'-struct','data');
end
end