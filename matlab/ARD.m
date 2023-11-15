function RDInfo = ARD(B1,B2,B3,x1,x2,x3,Parameters,varargin)
fprintf('######## ARD: Analyzing Reconnection Distribution ...\n');
% Load Parameters
ShowThresScalarProfile = Parameters.ARD_ShowThresScalarProfile;
ScalarThreshold = Parameters.ARD_ScalarThreshold;
NumSplitBlocks = Parameters.NumRAMBlock;
AnalyzeAllGrids = Parameters.ARD_AnalyzeAllGrids;
FixTrace = Parameters.ARD_FixTrace;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check input
if ~(isempty(varargin) || length(varargin) == 3)
    error('ERROR in ARD: length of N should be 0 (default) or 3 for 3 components of non-ideal term!');
end

[G.e1,G.e2,G.e3] = meshgrid(x1,x2,x3);
dx1 = mean(diff(x1));
dx2 = mean(diff(x2));
dx3 = mean(diff(x3));

% Calculate DB
fprintf('    >>>> Preparing Cache Data ...\n');
B.e1 = B1; clear B1;
B.e2 = B2; clear B2;
B.e3 = B3; clear B3;
DB = func_m_VectorJacobian(B,x1,x2,x3,FixTrace);

% Calculate Non-ideal term
if isempty(varargin) % N = J
    N.e1 = DB.e32 - DB.e23;
    N.e2 = DB.e13 - DB.e31;
    N.e3 = DB.e21 - DB.e12;
else
    N.e1 = varargin{1};
    N.e2 = varargin{2};
    N.e3 = varargin{3};
    clear varargin;
    if prod(size(N.e1) == size(B.e1)) ~= 1
        error('ERROR in ARD: N1 should has the same dimension as B');
    end
    if prod(size(N.e2) == size(B.e1)) ~= 1
        error('ERROR in ARD: N2 should has the same dimension as B');
    end
    if prod(size(N.e3) == size(B.e1)) ~= 1
        error('ERROR in ARD: N3 should has the same dimension as B');
    end
end

% Get data threshold
fprintf('    >>>> Selecting Subset Data ...\n');
NparaAll = func_m_VectorDot(N,func_m_VectorDirection(B)); % NparaAll=Epara

% Preview histogram of NparaNorm
if ShowThresScalarProfile == 1
    figure1 = figure;
    h_hist = histogram(abs(NparaAll(:)),500);
    bin_counts = h_hist.BinCounts;
    bin_center = h_hist.BinEdges;
    bin_center = 0.5*(bin_center(1:end-1)+bin_center(2:end));
    close(figure1);

    figure1 = figure;
    xlabel('$\left|E_\parallel\right|$','Interpreter','latex');
    ylabel('Counts','Interpreter','latex');
    plot(bin_center,bin_counts);
    set(gca,'yscale','log');
    RDInfo.Data = [];
    RDInfo.ExtraData = [];
    return;
end

% Select subset
if AnalyzeAllGrids == 0
    b = abs(NparaAll) >= ScalarThreshold;
else
    b = true(size(B.e1));
end

DB = func_m_TensorSubset(DB,b);
G = func_m_TensorSubset(G,b);

if Parameters.ARD_AnalyzeLocalEffects ~= 0
    Parameters.OutputExtraData = 1;
    % Calculate Local frame, Gradients of e1,e2,e3,Nperp/B
    % e1: generate one component to save RAM
    [e1,~,~] = func_VectorLocalFrame(B);
    De1 = func_m_VectorJacobian(e1,x1,x2,x3,0);
    e1 = func_m_TensorSubset(e1,b);
    De1 = func_m_TensorSubset(De1,b);
    % e2
    [~,e2,~] = func_VectorLocalFrame(B);
    De2 = func_m_VectorJacobian(e2,x1,x2,x3,0);
    e2 = func_m_TensorSubset(e2,b);
    De2 = func_m_TensorSubset(De2,b);
    % e3
    [~,~,e3] = func_VectorLocalFrame(B);
    De3 = func_m_VectorJacobian(e3,x1,x2,x3,0);
    De3 = func_m_TensorSubset(De3,b);

    % Nperp/B
    NperpOverB = func_m_VectorScaleProd(...
                      func_m_VectorSum(N,...
                                       func_m_VectorScaleProd(e3,-NparaAll)),...
                      1./func_m_VectorNorm(B)); % (N-NparaAll*b)/|B|
    e3 = func_m_TensorSubset(e3,b);

    DNperpOverB = func_m_VectorJacobian(NperpOverB,x1,x2,x3,0);
    clear NperpOverB;
    DNperpOverB = func_m_TensorSubset(DNperpOverB,b);

    % DNpara
    [DNpara.e1,DNpara.e2,DNpara.e3] = gradient(NparaAll,dx1,dx2,dx3);
    DNpara = func_m_TensorSubset(DNpara,b);

    % CurlN and BxCurlN
    CurlN = func_m_VectorCurl(N,x1,x2,x3);
    BxCurlN = func_m_VectorCross(B,CurlN);
    CurlN = func_m_TensorSubset(CurlN,b);
    BxCurlN = func_m_TensorSubset(BxCurlN,b);
else
    [e1,~,~] = func_VectorLocalFrame(B);
    e1 = func_m_TensorSubset(e1,b);
    [~,e2,~] = func_VectorLocalFrame(B);
    e2 = func_m_TensorSubset(e2,b);
    [~,~,e3] = func_VectorLocalFrame(B);
    e3 = func_m_TensorSubset(e3,b);
end
N = func_m_TensorSubset(N,b);
B = func_m_TensorSubset(B,b);
clear b;

% Main part
fprintf('    >>>> Classifying Local Magnetic Structures ...\n');
n_data = length(B.e1);
RDInfoStruc = func_ARD_RDInfoStruc();
RDInfo.Data = cast(zeros(n_data,RDInfoStruc.Data.n),'like',B.e1);
if Parameters.OutputExtraData ~= 0
    RDInfo.ExtraData = cast(zeros(n_data,RDInfoStruc.ExtraData.n),'like',B.e1);
else
    RDInfo.ExtraData = [];
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

    if Parameters.ARD_AnalyzeLocalEffects ~= 0
        cache_DNpara = func_m_TensorSubset(DNpara,idx_b);
        cache_De1 = func_m_TensorSubset(De1,idx_b);
        cache_De2 = func_m_TensorSubset(De2,idx_b);
        cache_De3 = func_m_TensorSubset(De3,idx_b);
        cache_DNperpOverB = func_m_TensorSubset(DNperpOverB,idx_b);
        cache_CurlN = func_m_TensorSubset(CurlN,idx_b);
        cache_BxCurlN = func_m_TensorSubset(BxCurlN,idx_b);
    else
        cache_DNpara = nan;
        cache_De1 = nan;
        cache_De2 = nan;
        cache_De3 = nan;
        cache_DNperpOverB = nan;
        cache_CurlN = nan;
        cache_BxCurlN = nan;
    end

    [RDInfo_Data,RDInfo_ExtraData] = infunc_ARD(... % The main function of ARD
        func_m_TensorSubset(B,idx_b),...% Subsets of main parameters
        func_m_TensorSubset(DB,idx_b),...
        func_m_TensorSubset(N,idx_b),...
        func_m_TensorSubset(G,idx_b),...
        func_m_TensorSubset(e1,idx_b),...
        func_m_TensorSubset(e2,idx_b),...
        func_m_TensorSubset(e3,idx_b),...
        cache_DNpara,... % Subset for calculating local effects
        cache_De1,...
        cache_De2,...
        cache_De3,...
        cache_DNperpOverB,...
        cache_CurlN,...
        cache_BxCurlN,...
        x1,x2,x3,NparaAll,...% Original data for interpolation
        Parameters);

    if Parameters.OutputExtraData ~= 0
        RDInfo.Data(idx_b,:) = RDInfo_Data;
        RDInfo.ExtraData(idx_b,:) = RDInfo_ExtraData;
    else
        RDInfo.Data(idx_b,:) = RDInfo_Data;
    end
end
clear RDInfo_Data RDInfo_ExtraData;
clear cache_DNpara cache_De1 cache_De2 cache_De3 cache_DNperpOverB cache_CurlN cache_BxCurlN; 

%Output interface
func_IO_Output(RDInfo,Parameters);

fprintf('    >>>> Done !!!\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [RDInfo, ExData]= infunc_ARD(B,DB,N,R0,e1,e2,e3,... % Subsets of main parameters
    DNpara,De1,De2,De3,DNperpOverB,... % Subset for calculating local effects
    CurlN,BxCurlN,... % Subsets only for outputing
    x1,x2,x3,NparaAll,... % Original data for interpolation
    Parameters) % Parameters

% Npara
Bnorm = func_m_VectorNorm(B);
Npara = func_m_VectorDot(N,B)./Bnorm;
% Snorm = func_m_MatrixSymmPart3D(DB);

% Generate transform matrix
[T,Tinv] = func_GenTransMatrixFromE1E2E3(e1,e2,e3);

% Transform DB into local frame
DBL = func_m_MatrixProd(func_m_MatrixProd(T,DB),Tinv);

%%%%%%%%%% Classify Bperp topology
BProjInfo = func_ARD_ClassifyBProj(DBL.e11,DBL.e22,DBL.e12,DBL.e21);

%%%%%%%%%% Determine Projected Local Maximum of NparaNorm
Is2DExtrema = ...
    func_ARD_Is2DExtrema(R0,Tinv,Npara,x1,x2,x3,NparaAll);

%%%%%%%%%% Determine Local 3D Maximum of NparaNorm
% Is3DExtrema = ...
%     func_ARD_Is3DExtrema(R0,Npara,x1,x2,x3,NparaAll);

% %%%%%%%%%% Calculate RatioDB3
% RatioDB3 = sqrt(DBL.e31.^2+DBL.e32.^2)./...
%     sqrt(DBL.e11.^2 + DBL.e22.^2 + DBL.e12.^2 + DBL.e21.^2);

%%%%%%%%%% Analyze local effects: in original frame
if Parameters.ARD_AnalyzeLocalEffects ~= 0
    [DNparaL,GammaL] = func_ARD_AnalyzeLocalEffects(B,N,DNpara,De1,De2,De3,DNperpOverB,T,Tinv);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Output
RDInfoStruc = func_ARD_RDInfoStruc();
RDInfo = cast(zeros(length(B.e1),RDInfoStruc.Data.n),'like',B.e1);
RDInfo(:,func_ARD_RDInfoStruc('x1')) = R0.e1;
RDInfo(:,func_ARD_RDInfoStruc('x2')) = R0.e2;
RDInfo(:,func_ARD_RDInfoStruc('x3')) = R0.e3;
RDInfo(:,func_ARD_RDInfoStruc('RDType')) = BProjInfo.RDType;
RDInfo(:,func_ARD_RDInfoStruc('Is2DExtrema')) = Is2DExtrema;
RDInfo(:,func_ARD_RDInfoStruc('EigAngle')) = BProjInfo.EigAngle;
RDInfo(:,func_ARD_RDInfoStruc('RatioMTrace')) = BProjInfo.RatioMTrace;
% RDInfo(:,func_ARD_RDInfoStruc('Is3DExtrema')) = Is3DExtrema;
% RDInfo(:,func_ARD_RDInfoStruc('Snorm')) = Snorm;
% RDInfo(:,func_ARD_RDInfoStruc('RatioDB3')) = RatioDB3;

% Save Extra Data
if Parameters.OutputExtraData ~= 0
    ExData = cast(zeros(length(B.e1),RDInfoStruc.ExtraData.n),'like',B.e1);
    ExData(:,func_ARD_RDInfoStruc('B0',1)) = Bnorm;
    ExData(:,func_ARD_RDInfoStruc('DB11',1)) = DBL.e11;
    ExData(:,func_ARD_RDInfoStruc('DB12',1)) = DBL.e12;
    ExData(:,func_ARD_RDInfoStruc('DB21',1)) = DBL.e21;
    ExData(:,func_ARD_RDInfoStruc('DB22',1)) = DBL.e22;
    ExData(:,func_ARD_RDInfoStruc('DB31',1)) = DBL.e31;
    ExData(:,func_ARD_RDInfoStruc('DB32',1)) = DBL.e32;
    ExData(:,func_ARD_RDInfoStruc('DB13',1)) = DBL.e13;
    ExData(:,func_ARD_RDInfoStruc('DB23',1)) = DBL.e23;
    ExData(:,func_ARD_RDInfoStruc('DB33',1)) = DBL.e33;
    ExData(:,func_ARD_RDInfoStruc('e11',1)) = e1.e1;
    ExData(:,func_ARD_RDInfoStruc('e12',1)) = e1.e2;
    ExData(:,func_ARD_RDInfoStruc('e13',1)) = e1.e3;
    ExData(:,func_ARD_RDInfoStruc('e21',1)) = e2.e1;
    ExData(:,func_ARD_RDInfoStruc('e22',1)) = e2.e2;
    ExData(:,func_ARD_RDInfoStruc('e23',1)) = e2.e3;
    ExData(:,func_ARD_RDInfoStruc('e31',1)) = e3.e1;
    ExData(:,func_ARD_RDInfoStruc('e32',1)) = e3.e2;
    ExData(:,func_ARD_RDInfoStruc('e33',1)) = e3.e3;
    if Parameters.ARD_AnalyzeLocalEffects ~= 0
        ExData(:,func_ARD_RDInfoStruc('DNpara1',1)) = DNparaL.e1;
        ExData(:,func_ARD_RDInfoStruc('DNpara2',2)) = DNparaL.e2;
        ExData(:,func_ARD_RDInfoStruc('Gamma1',1)) = GammaL.e1;
        ExData(:,func_ARD_RDInfoStruc('Gamma2',1)) = GammaL.e2;
        ExData(:,func_ARD_RDInfoStruc('CurlN',1)) = func_m_VectorNorm(CurlN);
        ExData(:,func_ARD_RDInfoStruc('BxCurlN',1)) = func_m_VectorNorm(BxCurlN);
    end
end
end
