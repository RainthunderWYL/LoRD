function NPInfo = ANP( B1,B2,B3,x1,x2,x3,Parameters)
% NPInfo.Data: x1,x2,x3,NPType;
% NPInfo.ExtraData: M11,M21,M31,M12,M22,M32,M13,M23,M33

BThres = Parameters.ANP_NullBThres;
NiterMax = Parameters.ANP_NumMaxIter;
NumSplitBlocks = Parameters.NumRAMBlock;
FixTrace  = Parameters.ANP_FixTrace;

x = x1;
y = x2;
z = x3;
B.x = B1;
B.y = B2;
B.z = B3;
nz = length(z);

fprintf('######## ANP: Analyzing Null Points ...\n');
fprintf('    >>>> Selecting NP Candidates ...\n');
NPCells = [];
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
        NPCells_i = func_ANP_Find_CellWithNullPoints3D(B.x(:,:,idx_z_b),...
                         B.y(:,:,idx_z_b),B.z(:,:,idx_z_b));
        NPCells_i(:,3) = NPCells_i(:,3) + idx_z_b(1) - 1;
        NPCells = [NPCells;...
                   NPCells_i];
    end
else
    NPCells = func_ANP_Find_CellWithNullPoints3D(B.x,B.y,B.z);
end

% Locate Null Point in Cell
fprintf('    >>>> Identifying Null Point Candidates ...\n');
Num_CellWithNP = size(NPCells,1);
NPInfo_All= [];
for i = 1:Num_CellWithNP
    jj = NPCells(i,1);
    ii = NPCells(i,2);
    kk = NPCells(i,3);

    xdomain = [x(ii),x(ii+1)];
    ydomain = [y(jj),y(jj+1)];
    zdomain = [z(kk),z(kk+1)];

    xyz0 = [0.5*sum(xdomain);0.5*sum(ydomain);0.5*sum(zdomain)];  % Set start point of root finding

    % Bx
    Bx_Cell_i(:,:,1) = [B.x(jj,ii,kk), B.x(jj+1,ii,kk);B.x(jj,ii+1,kk),B.x(jj+1,ii+1,kk)];
    Bx_Cell_i(:,:,2) = [B.x(jj,ii,kk+1), B.x(jj+1,ii,kk+1);B.x(jj,ii+1,kk+1),B.x(jj+1,ii+1,kk+1)];
    [cBx,~] = infunc_LinearInterp3(Bx_Cell_i,xdomain,ydomain,zdomain);
    % By
    By_Cell_i(:,:,1) = [B.y(jj,ii,kk), B.y(jj+1,ii,kk);B.y(jj,ii+1,kk),B.y(jj+1,ii+1,kk)];
    By_Cell_i(:,:,2) = [B.y(jj,ii,kk+1), B.y(jj+1,ii,kk+1);B.y(jj,ii+1,kk+1),B.y(jj+1,ii+1,kk+1)];
    [cBy,~] = infunc_LinearInterp3(By_Cell_i,xdomain,ydomain,zdomain);
    % Bz
    Bz_Cell_i(:,:,1) = [B.z(jj,ii,kk), B.z(jj+1,ii,kk);B.z(jj,ii+1,kk),B.z(jj+1,ii+1,kk)];
    Bz_Cell_i(:,:,2) = [B.z(jj,ii,kk+1), B.z(jj+1,ii,kk+1);B.z(jj,ii+1,kk+1),B.z(jj+1,ii+1,kk+1)];
    [cBz,~] = infunc_LinearInterp3(Bz_Cell_i,xdomain,ydomain,zdomain);

    [func_F,func_J] = infunc_GenRootFindingFunctions(cBx,cBy,cBz);

    icount = 0;
    xyz = xyz0;
    while true
        RF_F = func_F(xyz);
        RF_J = func_J(xyz);
        normF = norm(RF_F);

        % Check the RCOND of RF_J
        RCOND_RF_J = rcond(RF_J);
        if RCOND_RF_J<eps || isnan(RCOND_RF_J)
            if normF < BThres % Singular M with small B is Neutral line
                RFLabel = 0;
                break;
            else
                RFLabel = 1; % Abandon singular M with finite B
                break;
            end
        end

        deltaXYZ = RF_J\(-RF_F);

        if normF < BThres
            RFLabel = 0;
            break;
        elseif icount > NiterMax
            RFLabel = 1;
            break;
        else
            xyz = xyz + deltaXYZ;
            icount = icount+1;
        end
    end

    if infunc_indomain(xyz,xdomain,ydomain,zdomain)
        if RFLabel == 0
            M = func_J(xyz);
            % Calculate variables at this null point
            NPInfo_All = [NPInfo_All;...
                xyz(1),xyz(2),xyz(3),nan,...
                reshape(M,[1,9])];
        end
    end

end

fprintf('    >>>> Cleaning Adjacent Points...\n');
NPInfo_Cache = func_ANP_CleanAdjacentPoints(NPInfo_All,x1,x2,x3);

fprintf('    >>>> Classifying Null Points ...\n');
if ~isempty(NPInfo_Cache)
    if Parameters.ANP_NPJacobianMethod == 1
        NPInfo_Cache = func_ANP_UpdateDB(NPInfo_Cache,B1,B2,B3,x1,x2,x3);
    end
    NPInfo_Cache = func_ANP_IdentifyNullPointType(NPInfo_Cache,FixTrace);
end

% Output
NPInfo.Data = [];
NPInfo.ExtraData = [];
if ~isempty(NPInfo_Cache)
    NPInfo.Data = NPInfo_Cache(:,1:4);
    if Parameters.OutputExtraData ~= 0
        NPInfo.ExtraData = NPInfo_Cache(:,5:end);
    end
end

%Output interface
func_ANP_Output(NPInfo,Parameters);

fprintf('    >>>> ANP: Done!\n');
end

function [F,J] = infunc_GenRootFindingFunctions(c1,c2,c3)
F = @(x)  [c1(1)*x(1)*x(2)*x(3) + c1(2)*x(2)*x(3) + c1(3)*x(1)*x(3) + c1(4)*x(1)*x(2) + c1(5)*x(1) + c1(6)*x(2) + c1(7)*x(3) + c1(8);...
    c2(1)*x(1)*x(2)*x(3) + c2(2)*x(2)*x(3) + c2(3)*x(1)*x(3) + c2(4)*x(1)*x(2) + c2(5)*x(1) + c2(6)*x(2) + c2(7)*x(3) + c2(8);...
    c3(1)*x(1)*x(2)*x(3) + c3(2)*x(2)*x(3) + c3(3)*x(1)*x(3) + c3(4)*x(1)*x(2) + c3(5)*x(1) + c3(6)*x(2) + c3(7)*x(3) + c3(8);];

J = @(x) [c1(1)*x(2)*x(3)+c1(3)*x(3)+c1(4)*x(2)+c1(5), c1(1)*x(1)*x(3)+c1(2)*x(3)+c1(4)*x(1)+c1(6), c1(1)*x(1)*x(2)+c1(2)*x(2)+c1(3)*x(1)+c1(7);...
    c2(1)*x(2)*x(3)+c2(3)*x(3)+c2(4)*x(2)+c2(5), c2(1)*x(1)*x(3)+c2(2)*x(3)+c2(4)*x(1)+c2(6), c2(1)*x(1)*x(2)+c2(2)*x(2)+c2(3)*x(1)+c2(7);...
    c3(1)*x(2)*x(3)+c3(3)*x(3)+c3(4)*x(2)+c3(5), c3(1)*x(1)*x(3)+c3(2)*x(3)+c3(4)*x(1)+c3(6), c3(1)*x(1)*x(2)+c3(2)*x(2)+c3(3)*x(1)+c3(7)];
end

function indomain = infunc_indomain(X,xlim,ylim,zlim)
indomain = X(1)<=xlim(2) && X(1)>=xlim(1) && X(2)<=ylim(2) && X(2)>=ylim(1) && X(3)<=zlim(2) && X(3)>=zlim(1); % NP in domain
end

function [c,func] = infunc_LinearInterp3(f,xdomain,ydomain,zdomain)
% f(:,:,1) = [f111,f121; f211,f221];
% f(:,:,2) = [f112,f122; f212,f222];
% f111 = f(iy,ix,iz), f121 = f(iy+1,ix,iz), f211 = f(iy,ix+1,iz), f221 = f(iy+1,ix+1,iz)
% f112 = f(iy,ix,iz+1), f122 = f(iy+1,ix,iz+1), f212 = f(iy,ix+1,iz+1), f222 = f(iy+1,ix+1,iz+1)

% func(x,y,z) = c1*x*y*z + c2*y*z + c3*x*z + c4*x*y + c5*x + c6*y + c7*z + c8
% Notice: this function is approximately equivalent to the interpn (grid is ndgrid)

f111 = f(1,1,1);
f121 = f(1,2,1);
f211 = f(2,1,1);
f221 = f(2,2,1);

f112 = f(1,1,2);
f122 = f(1,2,2);
f212 = f(2,1,2);
f222 = f(2,2,2);

x1 = xdomain(1);
x2 = xdomain(2);
y1 = ydomain(1);
y2 = ydomain(2);
z1 = zdomain(1);
z2 = zdomain(2);

lx = x2 - x1;
ly = y2 - y1;
lz = z2 - z1;
lxlylz = lx*ly*lz;

c = zeros(1,8);

c(1) = (-f111 + f112 + f121 - f122 + f211 - f212 - f221 + f222)/lxlylz;
c(2) = ((-f211 + f212 + f221 - f222)*x1 + (f111 - f112 - f121 + f122)*x2)/lxlylz;
c(3) = ((-f121 + f122 + f221 - f222)*y1 + (f111 - f112 - f211 + f212)*y2)/lxlylz;
c(4) = ((-f112 + f122 + f212 - f222)*z1 + (f111 - f121 - f211 + f221)*z2)/lxlylz;
c(5) = ((-f122 + f222)*y1*z1 + (f112 - f212)*y2*z1 + (f121 - f221)*y1*z2 + (f211 - f111)*y2*z2)/lxlylz;
c(6) = ((-f212 + f222)*x1*z1 + (f112 - f122)*x2*z1 + (f211 - f221)*x1*z2 + (f121 - f111)*x2*z2)/lxlylz;
c(7) = ((-f221 + f222)*x1*y1 + (f121 - f122)*x2*y1 + (f211 - f212)*x1*y2 + (f112 - f111)*x2*y2)/lxlylz;
c(8) = (-(f222*x1*y1*z1) + f122*x2*y1*z1 + f212*x1*y2*z1 - f112*x2*y2*z1 + f221*x1*y1*z2 - f121*x2*y1*z2 - f211*x1*y2*z2 + f111*x2*y2*z2)/lxlylz;

func = @(x)  c(1)*x(1)*x(2)*x(3) + c(2)*x(2)*x(3) + c(3)*x(1)*x(3) + c(4)*x(1)*x(2) + c(5)*x(1) + c(6)*x(2) + c(7)*x(3) + c(8);
end