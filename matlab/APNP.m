function PNPInfo = APNP(B1,B2,B3,x1,x2,x3,Parameters)
fprintf('######## APNP: Analyzing Projected Null Points ...\n');

fprintf('    >>>> Generating Locale Frames ...\n');
e3.e1 = single(zeros(size(B1)));
e3.e2 = single(zeros(size(B1)));
e3.e3 = single(zeros(size(B1)));
% Determine the direction of "guide field"
NumSplitBlocks = Parameters.NumRAMBlock;
[~,~,nz] = size(B1);
if nz>1 % Split into NumSplitBlocks Blocks on z direction to avoid RAM damage
    if NumSplitBlocks<=0
        NumSplitBlocks = 1;
    end
    nz_block = fix(nz/NumSplitBlocks);
    for ib = 1:NumSplitBlocks
        if ib == NumSplitBlocks
            idx_z_b = ((ib-1)*nz_block + 1):nz;
        else
            idx_z_b = ((ib-1)*nz_block + 1):(ib*nz_block);
        end

        [e3.e1(:,:,idx_z_b),e3.e2(:,:,idx_z_b),e3.e3(:,:,idx_z_b)] = ...
            func_APNP_ProjectPlaneDirection(...
            B1(:,:,idx_z_b),B2(:,:,idx_z_b),B3(:,:,idx_z_b),...
            x1,x2,x3(idx_z_b),Parameters);
    end
else
    [e3.e1,e3.e2,e3.e3] = func_APNP_ProjectPlaneDirection(...
        B1,B2,B3,x1,x2,x3,Parameters);
end

[e1,e2,e3] = func_VectorLocalFrame(e3,1);

%%% Locate Singular Points
SPInfo_All = func_APNP_LocateProjectSingluarPoints(B1,B2,B3,x1,x2,x3,e1,e2,e3,Parameters);

fprintf('\n    >>>> Cleaning Adjacent Points...\n');
SPInfo = func_ANP_CleanAdjacentPoints(SPInfo_All,x1,x2,x3);

% Output
PNPInfo.Data = [];
PNPInfo.ExtraData = [];
if ~isempty(SPInfo)
    PNPInfo.Data = SPInfo(:,1:4);
    if Parameters.OutputExtraData ~= 0
        PNPInfo.ExtraData = SPInfo(:,5:end);
    end
end
fprintf('    >>>> APNP: Done!!!\n');
end