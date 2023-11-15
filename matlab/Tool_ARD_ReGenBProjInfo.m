function out = Tool_ARD_ReGenBProjInfo(FileFullPathTemplate,IdxMinMax,Stride)
% Need Test!!!
for idx = IdxMinMax(1):Stride:IdxMinMax(2)
    filename = sprintf(FileFullPathTemplate,idx);

    data = load(filename);
    
    SLInfo = data.SLInfo;

    DB11 = Tool_ARD_ReadSLInfo(SLInfo,'DB11',1);
    DB22 = Tool_ARD_ReadSLInfo(SLInfo,'DB22',1);
    DB12 = Tool_ARD_ReadSLInfo(SLInfo,'DB12',1);
    DB21 = Tool_ARD_ReadSLInfo(SLInfo,'DB21',1);

    BProjInfo = func_ARD_ClassifyBProj(DB11,DB22,DB12,DB21);

    SLInfo.Data(:,func_ARD_SLInfoStruc('SLType')) = BProjInfo.SLType;
    SLInfo.Data(:,func_ARD_SLInfoStruc('EigAngle')) = BProjInfo.EigAngle;
    SLInfo.Data(:,func_ARD_SLInfoStruc('RatioMTrace')) = BProjInfo.RatioMTrace;

    data.SLInfo = SLInfo;
    save(filename,'-struct','data');
end
end