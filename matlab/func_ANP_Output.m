function status = func_ANP_Output(NPInfo,Parameters)
% Parameters.OutputType = -1: don't output; 0: .mat; 1: csv;
OutputType = Parameters.OutputType;
OutputDir = Parameters.OutputDir;
OutputLabel = Parameters.OutputLabel;

if OutputType<0 % Do not output
    status = 0;
else
    OutputName_Header = sprintf('%s/NPInfo',OutputDir);
end

if OutputType == 0
    OutputName = sprintf('%s_%s.mat',OutputName_Header,OutputLabel);
    save(OutputName,'NPInfo','Parameters');
end

if OutputType == 1
    AllFields = func_ANP_NPInfoStruc();
    %Output NPInfo.Data
    OutputName = sprintf('%s_Data%s.csv',OutputName_Header,OutputLabel);
    tmp_table = array2table(NPInfo.Data);
    tmp_table.Properties.VariableNames = AllFields.Data.fields;
    writetable(tmp_table,OutputName);

    %Output NPInfo.ExtraData
    if Parameters.OutputExtraData == 1
        OutputName = sprintf('%s_ExtraData%s.csv',OutputName_Header,OutputLabel);
        tmp_table = array2table(NPInfo.ExtraData);
        tmp_table.Properties.VariableNames = AllFields.ExtraData.fields;
        writetable(tmp_table,OutputName);
    end
end

end
