function status = func_IO_Output(RDInfo,Parameters)
% Parameters.OutputType = -1: don't output; 0: .mat; 1: csv;
OutputType = Parameters.OutputType;
OutputDir = Parameters.OutputDir;
OutputLabel = Parameters.OutputLabel;

if OutputType<0 % Do not output
    status = 0;
else
    OutputName_Header = sprintf('%s/RDInfo',OutputDir);
end

if OutputType == 0
    OutputName = sprintf('%s_%s.mat',OutputName_Header,OutputLabel);
    save(OutputName,'RDInfo','Parameters');
end

if OutputType == 1
    AllFields = func_ARD_RDInfoStruc();
    %Output RDInfo.Data
    OutputName = sprintf('%s_Data%s.csv',OutputName_Header,OutputLabel);
    tmp_table = array2table(RDInfo.Data);
    tmp_table.Properties.VariableNames = AllFields.Data.fields;
    writetable(tmp_table,OutputName);

    %Output RDInfo.ExtraData
    if Parameters.OutputExtraData == 1
        OutputName = sprintf('%s_ExtraData%s.csv',OutputName_Header,OutputLabel);
        tmp_table = array2table(RDInfo.ExtraData);
        tmp_table.Properties.VariableNames = AllFields.ExtraData.fields;
        writetable(tmp_table,OutputName);
    end
end

end
