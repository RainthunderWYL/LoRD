%% Instructions
% This script is a demo for usage of the ARD function.
% Please download the input magnetic field files via urls below:
%   1. OneDrive: https://1drv.ms/u/s!AsKfeg1K3uHvix273n6-TUCVV3SI?e=LUvb7a
%   2. Google Drive: https://drive.google.com/file/d/1zIH_6Cr8yErBFqvr2rQvQM4eeS2e1_Cm/view?usp=sharing
%   3. OneDrive for Nanjing University: https://njuedu-my.sharepoint.cn/:u:/g/personal/ky2211909_365_nju_edu_cn/EZMSEyZEkoRCsEUEkqjH77kBOhjeqBapAWiRvfkEpYN2Ww?e=6A7td3
%   4. Baidu Disk (Extraction code: otbq): https://pan.baidu.com/s/14NJ24LhBh6jwCGBFRkpEgQ

clear;
%% Set path of ARD
path('/home/wyl/DATA/Codes/Projects/LoRD/matlab',path); % !!! Use your own path of LoRD

%% Load Data
DIM = [330,260,280]; % [n1,n2,n3]
Precision = 'double';
LIM = [-235,235,-185,185,0,400]; % [xmin,xmax,ymin,ymax,zmin,zmax]
[B1,x1,x2,x3] = Tool_LoadData_Bin('B1.bin',DIM,Precision,LIM);
[B2,~,~,~] = Tool_LoadData_Bin('B2.bin',DIM,Precision,LIM);
[B3,~,~,~] = Tool_LoadData_Bin('B3.bin',DIM,Precision,LIM);

%% Set Parameters
Parameters.NumRAMBlock = 4;    % Number of data blocks for saving RAM
Parameters.OutputType = 1; % -1: No output file; 0: mat-file; 1: csv-file
Parameters.OutputDir = '.';
Parameters.OutputLabel = 'Demo';
Parameters.OutputExtraData = 1;    % Output Extra Data

Parameters.ARD_AnalyzeAllGrids = 0;    % Analyze all grids
Parameters.ARD_ShowThresScalarProfile = 0;    % Depict histogram of DataThres without running ARD
Parameters.ARD_ScalarThreshold = 0.5;    % Threshold of DataThres
Parameters.ARD_FixTrace = 0;    % Fix trace of Jacobian when needed B2 M = M-1/3*trM*I
Parameters.ARD_AnalyzeLocalEffects = 0;

%% Run ARD
RDInfo = ARD(B1,B2,B3,x1,x2,x3,Parameters);

%% Post analyze
% Tool_ARD_ReadRDInfo is a function for reading from RDInfo structure
x_RD = Tool_ARD_ReadRDInfo(RDInfo,'x1');
y_RD = Tool_ARD_ReadRDInfo(RDInfo,'x2');
z_RD = Tool_ARD_ReadRDInfo(RDInfo,'x3');
RDType = Tool_ARD_ReadRDInfo(RDInfo,'RDType');
Is2DExtrema = logical(Tool_ARD_ReadRDInfo(RDInfo,'Is2DExtrema'));
B0_RD = Tool_ARD_ReadRDInfo(RDInfo,'B0','Extra'); % Read from RDInfo.ExtraData
[X1,X2,X3] = meshgrid(x1,x2,x3);

%% Draw figure
figure1 = figure('Units','centimeters','Position',[1 1 15 8.5]);
axes1 = axes('Parent',figure1,'Position',[0.1,0.07,0.85,0.95]);
hold(axes1,'on');
box(axes1,'on');
grid(axes1,'on');
hal = xlabel('$x_1\,\left(\mathrm{Mm}\right)$','Interpreter','latex');
set(hal,'Position',[0,-250,0])
hal = ylabel('$x_2\,\left(\mathrm{Mm}\right)$','Interpreter','latex');
set(hal,'Position',[-270,0,0])
hal = zlabel('$x_3\,\left(\mathrm{Mm}\right)$','Interpreter','latex');
set(hal,'Position',[-240,180,100])

b_RD = RDType == 1 & Is2DExtrema; % X-type: 2D extremal Epara
plot3(x_RD(b_RD),y_RD(b_RD),z_RD(b_RD),'.','MarkerSize',2);
b_RD = RDType == 2 | RDType == 3 & Is2DExtrema; % O-type: 2D extremal Epara
plot3(x_RD(b_RD),y_RD(b_RD),z_RD(b_RD),'.','MarkerSize',1);

hs = slice(X1,X2,X3,B3,[],[],[0]); 
set(hs,'EdgeColor','none');
colormap('gray');

set(axes1,'XLim',[-200,200],'YLim',[-180,180],'ZLim',[0,200],...
    'DataAspectRatio',[1,1,1],...
    'TickDir','in','layer','top','TickLabelInterpreter','latex',...
    'FontSize',8,'Projection','perspective', 'BoxStyle','back');
view(axes1,[-35 12]);

%% Trace several field lines at the flux rope top
% Define initial positions (randomly sampled within a sphere):
R0_O = [-41.47 4.29 131.56]; % Position of sphere center
Radius_O = 20; % Radius of sphere
n_sample_O = 20; % Number of initial samplings
LineConf.Len = 10000;
LineConf.Color = 'k'; % Field-line color
LineConf.Style = '-'; % Field-line style
LineConf.Width = 0.5; % Field-line width

% The field-line preview function:
hflO = Tool_PreviewFieldLines(B1,B2,B3,x1,x2,x3,'s',...
                              R0_O,Radius_O,n_sample_O,LineConf);
hold off;
