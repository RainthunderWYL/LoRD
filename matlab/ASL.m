function SLInfo = ASL(B1,B2,B3,x1,x2,x3,Parameters,varargin)
fprintf('######## ASL: Analyzing Singular Lines ...\n');
% Load Parameters
ShowThresScalorProfile = Parameters.ASL_ShowThresScalorProfile;
ScalorThreshold = Parameters.ASL_ScalorThreshold;
NumSplitBlocks = Parameters.NumRAMBlock;
AnalyzeAllGrids = Parameters.ASL_AnalyzeAllGrids;
FixTrace = Parameters.ASL_FixTrace;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(varargin)>1
    error('ERROR in ASL: length of varargin <= 1');
end
[G.e1,G.e2,G.e3] = meshgrid(x1,x2,x3);

% Calculate DB
fprintf('    >>>> Preparing Cache Data ...\n');
B.e1 = B1; clear B1;
B.e2 = B2; clear B2;
B.e3 = B3; clear B3;
DB = func_MagneticJacobian(B.e1,B.e2,B.e3,x1,x2,x3,FixTrace);

% Get data threshold
fprintf('    >>>> Selecting Subset Data ...\n');
if nargin == 7 % Use default Jpara
    J.e1 = DB.e32 - DB.e23;
    J.e2 = DB.e13 - DB.e31;
    J.e3 = DB.e21 - DB.e12;
    data_thres = abs(func_m_VectorDot(J,B)./func_m_VectorNorm(B)); % |Jpara|
else
    data_thres = varargin{1};
    if prod(size(data_thres) == size(B.e1)) ~= 1
        error('ERROR in ASL: threshold data should be the same size as B1');
    end
end

% Preview histogram of data_thres
if ShowThresScalorProfile == 1
    figure1 = figure;
    h_hist = histogram(data_thres(:),500);
    bin_counts = h_hist.BinCounts;
    bin_center = h_hist.BinEdges;
    bin_center = 0.5*(bin_center(1:end-1)+bin_center(2:end));
    close(figure1);

    figure1 = figure;
    xlabel('Threshold Scalor Value','Interpreter','latex');
    ylabel('Counts','Interpreter','latex');
    plot(bin_center,bin_counts);
    set(gca,'yscale','log');
    SLInfo.Data = [];
    SLInfo.ExtraData = [];
    return;
end

% Select subset
if AnalyzeAllGrids == 0
    b = data_thres >= ScalorThreshold;
else
    b = true(size(B.e1));
end
B = func_m_TensorSubset(B,b);
DB = func_m_TensorSubset(DB,b);
G = func_m_TensorSubset(G,b);
DataThres = data_thres(b);
clear b;

% Main part
fprintf('    >>>> Classifying Local Magnetic Structures ...\n');
n_data = length(B.e1);
SLInfoStruc = func_ASL_SLInfoStruc();
SLInfo.Data = cast(zeros(n_data,SLInfoStruc.Data.n),'like',B.e1);
if Parameters.OutputExtraData ~= 0
    SLInfo.ExtraData = cast(zeros(n_data,SLInfoStruc.ExtraData.n),'like',B.e1);
else
    SLInfo.ExtraData = [];
end

if NumSplitBlocks<=0
    NumSplitBlocks = 1;
end
n_block = fix(n_data/NumSplitBlocks);
for ib = 1:NumSplitBlocks
    if ib == NumSplitBlocks
        idx_b = ((ib-1)*n_block + 1):n_data;
    else
        idx_b = ((ib-1)*n_block + 1):(ib*n_block);
    end

    if Parameters.OutputExtraData ~= 0
        [SLInfo.Data(idx_b,:),SLInfo.ExtraData(idx_b,:)] = infunc_ASL(func_m_TensorSubset(B,idx_b),...
            func_m_TensorSubset(DB,idx_b),func_m_TensorSubset(G,idx_b),DataThres(idx_b),...
            x1,x2,x3,data_thres,Parameters);
    else
        [SLInfo.Data(idx_b,:),~] = infunc_ASL(func_m_TensorSubset(B,idx_b),...
            func_m_TensorSubset(DB,idx_b),func_m_TensorSubset(G,idx_b),DataThres(idx_b),...
            x1,x2,x3,data_thres,Parameters);
    end
end

%Output interface
func_IO_Output(SLInfo,Parameters);

fprintf('    >>>> Done !!!\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [SLInfo, ExData]= infunc_ASL(B,DB,R0,DataThres, ... % Subset
    x1,x2,x3,data_thres,... % Origin data for interp
    Parameters)
% Local frame: e1, e2, e3
[e1,e2,e3] = func_VectorLocalFrame(B);

% Generate transform matrix
[T,Tinv] = func_GenTransMatrixFromE1E2E3(e1,e2,e3);

% Transform DB into local frame
DBL = func_m_MatrixProd(func_m_MatrixProd(T,DB),Tinv);

%%%%%%%%%% Classify Bperp topology
BProjInfo = func_ASL_ClassifyBProj(DBL.e11,DBL.e22,DBL.e12,DBL.e21);

%%%%%%%%%% Determine Projected Local Maximum of data_thres
Is2DExtrema = ...
    func_ASL_Is2DExtrema(R0,Tinv,DataThres,x1,x2,x3,data_thres);

%%%%%%%%%% Determine Local 3D Maximum of data_thres
Is3DExtrema = ...
    func_ASL_Is3DExtrema(R0,DataThres,x1,x2,x3,data_thres);

%%%%%%%%%% Output
SLInfoStruc = func_ASL_SLInfoStruc();
SLInfo = cast(zeros(length(B.e1),SLInfoStruc.Data.n),'like',B.e1);
SLInfo(:,func_ASL_SLInfoStruc('x1')) = R0.e1;
SLInfo(:,func_ASL_SLInfoStruc('x2')) = R0.e2;
SLInfo(:,func_ASL_SLInfoStruc('x3')) = R0.e3;
SLInfo(:,func_ASL_SLInfoStruc('SLType')) = BProjInfo.SLType;
SLInfo(:,func_ASL_SLInfoStruc('Is2DExtrema')) = Is2DExtrema;
SLInfo(:,func_ASL_SLInfoStruc('Is3DExtrema')) = Is3DExtrema;
SLInfo(:,func_ASL_SLInfoStruc('EigAngle')) = BProjInfo.EigAngle;
SLInfo(:,func_ASL_SLInfoStruc('RatioMTrace')) = BProjInfo.RatioMTrace;

% Save Extra Data
if Parameters.OutputExtraData ~= 0
    ExData = cast(zeros(length(B.e1),SLInfoStruc.ExtraData.n),'like',B.e1);
    ExData(:,func_ASL_SLInfoStruc('DB11',1)) = DBL.e11;
    ExData(:,func_ASL_SLInfoStruc('DB12',1)) = DBL.e12;
    ExData(:,func_ASL_SLInfoStruc('DB21',1)) = DBL.e21;
    ExData(:,func_ASL_SLInfoStruc('DB22',1)) = DBL.e22;
    ExData(:,func_ASL_SLInfoStruc('e11',1))  = e1.e1;
    ExData(:,func_ASL_SLInfoStruc('e12',1))  = e1.e2;
    ExData(:,func_ASL_SLInfoStruc('e13',1))  = e1.e3;
    ExData(:,func_ASL_SLInfoStruc('e21',1))  = e2.e1;
    ExData(:,func_ASL_SLInfoStruc('e22',1))  = e2.e2;
    ExData(:,func_ASL_SLInfoStruc('e23',1)) = e2.e3;
    ExData(:,func_ASL_SLInfoStruc('e31',1)) = e3.e1;
    ExData(:,func_ASL_SLInfoStruc('e32',1)) = e3.e2;
    ExData(:,func_ASL_SLInfoStruc('e33',1)) = e3.e3;
else
    ExData = [];
end
end
