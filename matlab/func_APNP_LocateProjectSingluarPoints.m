function SPInfo = func_APNP_LocateProjectSingluarPoints(...
    B1,B2,B3,x,y,z,e1,e2,e3,Parameters)
% SPInfo: [x,y,z,SPType_SCond,M11,M12,M21,M22,e11,e12,e13,e21,e22,e23,e31,e32,e33]

% Load Parameters
NumRefine = Parameters.APNP_MSP_RefineLevel; 
MSP_ScaleRatio = Parameters.APNP_MSP_ScaleRatio;
NumSplitBlocks = Parameters.NumRAMBlock;
NullBThres = Parameters.NullBThres;
NumMaxIter = Parameters.NumMaxIter;
SMALL_NUMBER = Parameters.SMALL_NUMBER;


[X,Y,Z] = meshgrid(x,y,z);
[~,~,nz] = size(B1);
dL = [mean(diff(x)),mean(diff(y)),mean(diff(z))];
MaxCellProjectScale = norm(dL)*MSP_ScaleRatio; %!!!!!!!!!!!!!!!!!!!!!
% MaxCellProjectScale = mean([dx,dy,dz]);

% Generate local grid
if NumRefine<0
    NumRefine = 0;
end
NumLocalGrids = 2*(NumRefine+1)+1;
x1_s = linspace(-0.5*MaxCellProjectScale,0.5*MaxCellProjectScale,NumLocalGrids);
x2_s = x1_s;
[X1_s,X2_s] = meshgrid(x1_s,x2_s);

% Determine NP-Cell Candidates by Tensor Representation
fprintf('    >>>> Selecting PNP Candidates ...\n');
b_CellWithNP = false(size(B1));
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
        b_CellWithNP(:,:,idx_z_b) = ...
            infunc_ExtractCandidates4SigularPoints( ...
            B1(:,:,idx_z_b),B2(:,:,idx_z_b),B3(:,:,idx_z_b), ...
            X(:,:,idx_z_b),Y(:,:,idx_z_b),Z(:,:,idx_z_b), ...
            e1.e1(:,:,idx_z_b),e1.e2(:,:,idx_z_b),e1.e3(:,:,idx_z_b), ...
            e2.e1(:,:,idx_z_b),e2.e2(:,:,idx_z_b),e2.e3(:,:,idx_z_b),X1_s,X2_s);
    end
else
    b_CellWithNP = infunc_ExtractCandidates4SigularPoints(B1,B2,B3,X,Y,Z,...
        e1.e1,e1.e2,e1.e3,e2.e1,e2.e2,e2.e3,X1_s,X2_s);
end

% Extract Candidates
X_WithNP = X(b_CellWithNP);
Y_WithNP = Y(b_CellWithNP);
Z_WithNP = Z(b_CellWithNP);
E11_WithNP = e1.e1(b_CellWithNP);
E12_WithNP = e1.e2(b_CellWithNP);
E13_WithNP = e1.e3(b_CellWithNP);
E21_WithNP = e2.e1(b_CellWithNP);
E22_WithNP = e2.e2(b_CellWithNP);
E23_WithNP = e2.e3(b_CellWithNP);
E31_WithNP = e3.e1(b_CellWithNP);
E32_WithNP = e3.e2(b_CellWithNP);
E33_WithNP = e3.e3(b_CellWithNP);

num_candidates = length(X_WithNP);

% Sequently Processing NP-Cell Candidates
i_disp = floor(linspace(1,num_candidates,11));
idx_i_disp = 1;
fprintf('    >>>> Identifying PNP Candidates ... \n');
SPInfo = [];
% parpool(4);
for i = 1:num_candidates
    if i_disp(idx_i_disp) == i
        fprintf('%d%%...',int32(round(i/num_candidates*100)));
        idx_i_disp = idx_i_disp + 1;
    end

    e1i = [E11_WithNP(i),E12_WithNP(i),E13_WithNP(i)]';
    e2i = [E21_WithNP(i),E22_WithNP(i),E23_WithNP(i)]';
    e3i = [E31_WithNP(i),E32_WithNP(i),E33_WithNP(i)]';
    r0i = [X_WithNP(i),Y_WithNP(i),Z_WithNP(i)]';

    %Lab Cor
    X1_i = r0i(1) + X1_s*e1i(1) + X2_s*e2i(1);
    X2_i = r0i(2) + X1_s*e1i(2) + X2_s*e2i(2);
    X3_i = r0i(3) + X1_s*e1i(3) + X2_s*e2i(3);

    %Extract local B data to improve interp efficiency
    ix1 = infunc_GenCorIdx([min(X1_i(:)),max(X1_i(:))],x);
    ix2 = infunc_GenCorIdx([min(X2_i(:)),max(X2_i(:))],y);
    ix3 = infunc_GenCorIdx([min(X3_i(:)),max(X3_i(:))],z);

    B1_i = interp3(X(ix2(1):ix2(2),ix1(1):ix1(2),ix3(1):ix3(2)),...
        Y(ix2(1):ix2(2),ix1(1):ix1(2),ix3(1):ix3(2)),...
        Z(ix2(1):ix2(2),ix1(1):ix1(2),ix3(1):ix3(2)),...
        B1(ix2(1):ix2(2),ix1(1):ix1(2),ix3(1):ix3(2)),...
        X1_i,X2_i,X3_i,'linear');
    B2_i = interp3(X(ix2(1):ix2(2),ix1(1):ix1(2),ix3(1):ix3(2)),...
        Y(ix2(1):ix2(2),ix1(1):ix1(2),ix3(1):ix3(2)),...
        Z(ix2(1):ix2(2),ix1(1):ix1(2),ix3(1):ix3(2)),...
        B2(ix2(1):ix2(2),ix1(1):ix1(2),ix3(1):ix3(2)),...
        X1_i,X2_i,X3_i,'linear');
    B3_i = interp3(X(ix2(1):ix2(2),ix1(1):ix1(2),ix3(1):ix3(2)),...
        Y(ix2(1):ix2(2),ix1(1):ix1(2),ix3(1):ix3(2)),...
        Z(ix2(1):ix2(2),ix1(1):ix1(2),ix3(1):ix3(2)),...
        B3(ix2(1):ix2(2),ix1(1):ix1(2),ix3(1):ix3(2)),...
        X1_i,X2_i,X3_i,'linear');

    %B in Local Cor
    B1_i_s = B1_i*e1i(1) + B2_i*e1i(2) + B3_i*e1i(3);
    B2_i_s = B1_i*e2i(1) + B2_i*e2i(2) + B3_i*e2i(3);

    %Locate 3D singular points
    Is3D = 1;
    if sum(isnan(B1_i_s)) == 0 % Pass boundary cells with nan
        SPInfo_tmp = func_APNP_ClassifyPNP( x1_s,x2_s,B1_i_s,B2_i_s,NullBThres,NumMaxIter,SMALL_NUMBER,Is3D);
        if ~isempty(SPInfo_tmp)
            num_lsp = size(SPInfo_tmp,1);

            r_sp_i = SPInfo_tmp(:,1:2);
            r_lab_i = zeros(num_lsp,3);

            r_lab_i(:,1) = r0i(1) + r_sp_i(:,1)*e1i(1) + r_sp_i(:,2)*e2i(1);
            r_lab_i(:,2) = r0i(2) + r_sp_i(:,1)*e1i(2) + r_sp_i(:,2)*e2i(2);
            r_lab_i(:,3) = r0i(3) + r_sp_i(:,1)*e1i(3) + r_sp_i(:,2)*e2i(3);

            SPInfo = [SPInfo; ... 
                      [ r_lab_i, SPInfo_tmp(:,3:end), ...
                      repmat(e1i',[num_lsp,1]), ...
                      repmat(e2i',[num_lsp,1]), ...
                      repmat(e3i',[num_lsp,1])]...
                     ];
        end
    end
end
end

function IDX = infunc_GenCorIdx( XLIM, x)

IDX = zeros(2,1);

IdxAll = 1:length(x);
x0 = XLIM(1);
x1 = XLIM(2);

if x0<min(x)
    x0 = min(x);
end
b = x<=x0;
IDX(1) = IdxAll( IdxAll == max(IdxAll(b)) );

if x1>max(x)
    x1 = max(x);
end
b = x>=x1;
IDX(2) = IdxAll( IdxAll == min(IdxAll(b)) );

% if IDX(1) == IDX(2)
%     if x0 == min(x) && x1 ~= max(x)
%         IDX(2) = IDX(2)+1;
%     elseif x1 == max(x) && x0 ~= max(x)
%         IDX(1) = IDX(1)-1;
%     elseif x1 ~= max(x) && x0 ~= max(x)
%         IDX(1) = IDX(1)-1;
%         IDX(2) = IDX(2)+1;
%     else
%         error('infunc_GenCorIdx: max_x == min_x!!!');
%     end
% end
if IDX(1) == IDX(2)
    if IDX(1) == 1
        IDX(2) = IDX(2)+1;
    elseif IDX(2) == IdxAll(end)
        IDX(1) = IDX(1)-1;
    else
        IDX(1) = IDX(1)-1;
        IDX(2) = IDX(2)+1;
    end
end
end

function b_CellWithNP = infunc_ExtractCandidates4SigularPoints(B1,B2,B3,X,Y,Z,e11,e12,e13,e21,e22,e23,X1_s,X2_s)
[ny,nx,nz] = size(B1);
NumLocalGrids = size(X1_s,1);

E11_All = permute(repmat(e11,[1,1,1,NumLocalGrids,NumLocalGrids]),[4,5,1,2,3]);
E12_All = permute(repmat(e12,[1,1,1,NumLocalGrids,NumLocalGrids]),[4,5,1,2,3]);
E13_All = permute(repmat(e13,[1,1,1,NumLocalGrids,NumLocalGrids]),[4,5,1,2,3]);

E21_All = permute(repmat(e21,[1,1,1,NumLocalGrids,NumLocalGrids]),[4,5,1,2,3]);
E22_All = permute(repmat(e22,[1,1,1,NumLocalGrids,NumLocalGrids]),[4,5,1,2,3]);
E23_All = permute(repmat(e23,[1,1,1,NumLocalGrids,NumLocalGrids]),[4,5,1,2,3]);

X1_s_All = repmat(X1_s,[1,1,ny,nx,nz]);
X2_s_All = repmat(X2_s,[1,1,ny,nx,nz]);
R0_1_All = permute(repmat(X,[1,1,1,NumLocalGrids,NumLocalGrids]),[4,5,1,2,3]);
R0_2_All = permute(repmat(Y,[1,1,1,NumLocalGrids,NumLocalGrids]),[4,5,1,2,3]);
R0_3_All = permute(repmat(Z,[1,1,1,NumLocalGrids,NumLocalGrids]),[4,5,1,2,3]);

X1_Lab_All = R0_1_All + X1_s_All.*E11_All + X2_s_All.*E21_All;
X2_Lab_All = R0_2_All + X1_s_All.*E12_All + X2_s_All.*E22_All;
X3_Lab_All = R0_3_All + X1_s_All.*E13_All + X2_s_All.*E23_All;

clear X1_s_All X2_s_All R0_1_All R0_2_All R0_3_All;

B1_Lab_All = interp3(X,Y,Z,B1,X1_Lab_All,X2_Lab_All,X3_Lab_All,'linear');
B2_Lab_All = interp3(X,Y,Z,B2,X1_Lab_All,X2_Lab_All,X3_Lab_All,'linear');
B3_Lab_All = interp3(X,Y,Z,B3,X1_Lab_All,X2_Lab_All,X3_Lab_All,'linear');

B1_s_All = B1_Lab_All.*E11_All + B2_Lab_All.*E12_All + B3_Lab_All.*E13_All;
B2_s_All = B1_Lab_All.*E21_All + B2_Lab_All.*E22_All + B3_Lab_All.*E23_All;
clear E11_All E12_All E13_All E21_All E22_All E23_All;
clear B1_Lab_All B2_Lab_All B3_Lab_All;

iix = 1:NumLocalGrids-1;
iiy = 1:NumLocalGrids-1;
B1_s_Cell = zeros(NumLocalGrids-1,NumLocalGrids-1,ny,nx,nz,4);
B1_s_Cell(:,:,:,:,:,1) = B1_s_All(iiy,iix,:,:,:);
B1_s_Cell(:,:,:,:,:,2) = B1_s_All(iiy,iix+1,:,:,:);
B1_s_Cell(:,:,:,:,:,3) = B1_s_All(iiy+1,iix,:,:,:);
B1_s_Cell(:,:,:,:,:,4) = B1_s_All(iiy+1,iix+1,:,:,:);

b_B1_s_Null = (max(B1_s_Cell,[],6)).*(min(B1_s_Cell,[],6))<=0;
clear B1_s_Cell;

B2_s_Cell = zeros(NumLocalGrids-1,NumLocalGrids-1,ny,nx,nz,4);
B2_s_Cell(:,:,:,:,:,1) = B2_s_All(iiy,iix,:,:,:);
B2_s_Cell(:,:,:,:,:,2) = B2_s_All(iiy,iix+1,:,:,:);
B2_s_Cell(:,:,:,:,:,3) = B2_s_All(iiy+1,iix,:,:,:);
B2_s_Cell(:,:,:,:,:,4) = B2_s_All(iiy+1,iix+1,:,:,:);

b_B2_s_Null = (max(B2_s_Cell,[],6)).*(min(B2_s_Cell,[],6))<=0;
clear B2_s_Cell;

b_CellWithNP_All = b_B1_s_Null & b_B2_s_Null;
clear b_B1_s_Null b_B2_s_Null;

b_CellWithNP = reshape(sum(sum(b_CellWithNP_All,1),2),[ny,nx,nz])>0;
end

function b_indomain = infunc_indomain_Rijk(R,R0,dL)
% R(i,:) = [x1,x2,x3];
% R0 = [x1_0,x2_0,x3_0];
% dL = [dx1,dx2,dx3];
dL_half = 0.5*dL;
x1lim = [R0(1)-dL_half(1),R0(1)+dL_half(1)];
x2lim = [R0(2)-dL_half(2),R0(2)+dL_half(2)];
x3lim = [R0(3)-dL_half(3),R0(3)+dL_half(3)];

b_indomain = R(:,1)<x1lim(2) & R(:,1)>=x1lim(1) & ... 
             R(:,2)<x2lim(2) & R(:,2)>=x2lim(1) & ...
             R(:,3)<x3lim(2) & R(:,3)>=x3lim(1); % NP in domain
end
