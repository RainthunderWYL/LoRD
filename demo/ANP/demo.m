% This script is a demo for usage of the ANP function.
clear;
%% Set path of ANP
path('/home/wyl/DATA/Codes/Projects/LoRD/matlab',path); % !!! Use your own path of LoRD

%% Create discrete magnetic field Data: two potential nulls (positive/negative)
% Parameters
B0 = 1;
L0 = 1;
a = 0.5;
b = 0.75;
c = 0.25;
jsep = 0;
L = 1;
% Field functions
infunc_Bx = @(X,Y,Z) B0/L0*(X+c*X.*Z+b*Y.*Z-0.5*jsep*Y);
infunc_By = @(X,Y,Z) B0/L0*((2*a-c)*Y.*Z-(1+L*a)*Y+b*X.*Z+0.5*jsep*X);
infunc_Bz = @(X,Y,Z) B0/L0*(a*(L*Z-Z.^2)+0.5*c*X.^2+(a-0.5*c)*Y.^2+b*X.*Y);

% Grids
dL = 0.011;
x1 = -L0:dL:L0;
x2 = -L0:dL:L0;
x3 = -L0:dL:L0+L;
[GX,GY,GZ]=meshgrid(x1,x2,x3);

% Magnetic field
B1 = infunc_Bx(GX,GY,GZ);
B2 = infunc_By(GX,GY,GZ);
B3 = infunc_Bz(GX,GY,GZ);

%% Set Parameters
Parameters.NumRAMBlock = 1;    % Number of data blocks for saving RAM
Parameters.OutputType = 1; % -1: No output file; 0: mat-file; 1: csv-file
Parameters.OutputDir = '.';
Parameters.OutputLabel = 'Demo';
Parameters.OutputExtraData = 1;    % Output Extra Data

Parameters.ANP_NullBThres = eps('single');    % Lower threshold of Null magnetic strength
Parameters.ANP_NumMaxIter = 50;    % Maximum interation number for root-finding
Parameters.ANP_NPJacobianMethod = 1;    % If 1: First calculate Jacobian by central difference then interp 
Parameters.ANP_FixTrace = 1;    % Fix trace of Jacobian when needed by M = M-1/3*trM*I

%% Run ANP
NPInfo = ANP(B1,B2,B3,x1,x2,x3,Parameters);

%% Post analyze
% Tool_ANP_ReadRDInfo is a function for reading from NPInfo structure
x_NP = Tool_ANP_ReadNPInfo(NPInfo,'x1');
y_NP = Tool_ANP_ReadNPInfo(NPInfo,'x2');
z_NP = Tool_ANP_ReadNPInfo(NPInfo,'x3');
NPType = Tool_ANP_ReadNPInfo(NPInfo,'NPType');

%% Draw figure
figure1 = figure('Units','centimeters','Position',[1 1 12 15]);
axes1 = axes('Parent',figure1);
hold(axes1,'on');
box(axes1,'on');
grid(axes1,'on');
xlabel('$x_1\,\left(\mathrm{Mm}\right)$','Interpreter','latex');
ylabel('$x_2\,\left(\mathrm{Mm}\right)$','Interpreter','latex');
zlabel('$x_3\,\left(\mathrm{Mm}\right)$','Interpreter','latex');

% Plot Nulls 
plot3(x_NP,y_NP,z_NP,'*','MarkerSize',10,'Color','m');

% Plot a referece line
plot3([0,0],[0,0],[-L0,L0+L],'b--');

%% Trace several field lines near two nulls
% Define initial positions (randomly sampled within a sphere):
R0 = [0,0,0]; % Position of sphere center: Null 0
R1 = [0,0,1]; % Position of sphere center: Null 1

Radius = 0.3; % Radius of sphere
n_sample = 50; % Number of initial samplings
LineConf.Len = 10000;
LineConf.Color = [0.5,0.5,0.5]; % Field-line color
LineConf.Style = '-'; % Field-line style
LineConf.Width = 0.5; % Field-line width

% The field-line preview function:
hflO = Tool_PreviewFieldLines(B1,B2,B3,x1,x2,x3,'s',...
                              R0,Radius,n_sample,LineConf);
hflO = Tool_PreviewFieldLines(B1,B2,B3,x1,x2,x3,'s',...
                              R1,Radius,n_sample,LineConf);

set(axes1,'XLim',[-L0,L0],'YLim',[-L0,L0],'ZLim',[-L0,L0+L],...
    'DataAspectRatio',[1,1,1],...
    'TickDir','in','layer','top','TickLabelInterpreter','latex',...
    'FontSize',10,'Projection','perspective', 'BoxStyle','back');
view(axes1,[-32,18]);

hold off;