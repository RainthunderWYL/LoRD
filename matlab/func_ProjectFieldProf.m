function out = func_ProjectFieldProf(haxis,B1,B2,B3,x1,x2,x3,r0,e1,e2,e3,MSP_ScaleRatio,NumRefine,MSPDispScale)

r0i = r0;
e1i = e1;
e2i = e2;
e3i = e3;

r1 = r0i+e1i*MSPDispScale;
r2 = r0i+e2i*MSPDispScale;
r3 = r0i+e3i*MSPDispScale;
% Plot Local frame on haxis
hold(haxis,'on');
plot3(haxis,[r0i(1),r3(1)],[r0i(2),r3(2)],[r0i(3),r3(3)],'Color','k','LineStyle','-','LineWidth',2);
plot3(haxis,[r0i(1),r2(1)],[r0i(2),r2(2)],[r0i(3),r2(3)],'Color','k','LineStyle','--','LineWidth',1);
plot3(haxis,[r0i(1),r1(1)],[r0i(2),r1(2)],[r0i(3),r1(3)],'Color','k','LineStyle','--','LineWidth',1);
% hold(haxis,'off');

% Plot Local magnetic field on MSP
[X,Y,Z] = meshgrid(x1,x2,x3);
dL = [mean(diff(x1)),mean(diff(x2)),mean(diff(x3))];
MaxCellProjectScale = norm(dL)*MSP_ScaleRatio;

% Generate local grid
if NumRefine<0
    NumRefine = 0;
end
NumLocalGrids = 2*(NumRefine+1)+1;
x1_s = linspace(-0.5*MaxCellProjectScale,0.5*MaxCellProjectScale,NumLocalGrids);
x2_s = x1_s;
[X1_s,X2_s] = meshgrid(x1_s,x2_s);

%Lab Cor
X1_i = r0i(1) + X1_s*e1i(1) + X2_s*e2i(1);
X2_i = r0i(2) + X1_s*e1i(2) + X2_s*e2i(2);
X3_i = r0i(3) + X1_s*e1i(3) + X2_s*e2i(3);

B1_i = interp3(X,Y,Z,B1,X1_i,X2_i,X3_i,'linear');
B2_i = interp3(X,Y,Z,B2,X1_i,X2_i,X3_i,'linear');
B3_i = interp3(X,Y,Z,B3,X1_i,X2_i,X3_i,'linear');

%B in Local Cor
B1_i_s = B1_i*e1i(1) + B2_i*e1i(2) + B3_i*e1i(3);
B2_i_s = B1_i*e2i(1) + B2_i*e2i(2) + B3_i*e2i(3);

dL_s = mean(diff(x1_s));
step_scale = 0.1;
figure1 = figure('Position',[10,10,800,800]);
ha = axes('Parent',figure1);
hold on;

l_color = [0.5,0.5,0.5];

hfl = streamslice(X1_s,X2_s,B1_i_s,B2_i_s,'linear');
set(hfl,'Color',l_color,'LineWidth',0.5);

hold off;

out.hfl = hfl;
end