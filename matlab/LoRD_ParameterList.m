% Parameter List

% Public Parameters
Parameters.NumRAMBlock = 4;    % Number of data blocks for saving RAM
Parameters.SMALL_NUMBER = eps('single');    % Value is zero if below this value
Parameters.NullBThres = eps('single');    % Lower threshold of Null magnetic strength
Parameters.NumMaxIter = 50;    % Maximum interation number for root-finding

Parameters.OutputType = -1;
Parameters.OutputDir = '.';
Parameters.OutputLabel = '';
Parameters.OutputExtraData = 1;    % Output Extra Data

% ASL Parameters
Parameters.ASL_AnalyzeAllGrids = 0;    % Analyze all grids
Parameters.ASL_ShowThresScalorProfile = 0;    % Depict histogram of DataThres without running ASL
Parameters.ASL_ScalorThreshold = 150;    % Threshold of DataThres
Parameters.ASL_FixTrace = 0;    % Fix trace of Jacobian when needed by M = M-1/3*trM*I

% ANP Parameters
Parameters.ANP_NPJacobianMethod = 1;    % If 1: First calculate Jacobian by central difference then interp 
Parameters.ANP_FixTrace = 0;    % Fix trace of Jacobian when needed by M = M-1/3*trM*I

% APNP Parameters
Parameters.APNP_E3Type = 0;    % 0: MSP+J; 1: MSP; 2: B; -1: User defined constant
Parameters.APNP_MSP_RefineLevel = 0;    % Refine level of local Maximum Shear Plane
Parameters.APNP_MSP_ScaleRatio = 1.2;    % Refine level of local Maximum Shear Plane
Parameters.APNP_ConstProjectPlaneDirection = [0,0,1];  % Work only for APNP_E3Type = -1
Parameters.APNP_FixTrace = 0;    % Fix trace of Jacobian when needed by M = M-1/3*trM*I