function SPInfo = func_APNP_ClassifyPNP( x1,x2,B1,B2,BThres,NiterMax,SMALL_NUMBER,Is3D)
% SPInfo: [x, y, SPType, Tr(M), real(lambda1), imag(lambda1), real(lambda2), imag(lambda2)]
% SPType: 9 types as follows:
% 1: 2D_NeutralLine: lambda1 = lambda2 = 0
% 2: 2D_X:           trace-free & lambda1 = -lambda2 are real numbers
% 3: 2D_O:           trace-free & lambda1 = -lambda2 are purely imaginary numbers
% 4: 3D_NeutralLine: lambda1 = 0 & lambda2 ~= 0
% 5: 3D_X:           trace & lambda1 * lambda2 < 0 & Real number
% 6: 3D_Repel:       trace & lambda1 > 0 & lambda2 > 0 & Real number
% 7: 3D_Attract:     trace & lambda1 < 0 & lambda2 < 0 & Real number
% 8: 3D_O_Repel:     trace>0 & lambda1, lambda2 are complex number
% 9: 3D_O_Attract:   trace<0 & lambda1, lambda2 are complex number

NPCells = func_APNP_CellWithNullPoints2D(B1,B2);
% [X1,X2] = meshgrid(x1,x2);

% Locate Singular Points in cells, determine their types
Num_CellWithNP = size(NPCells,1);
SPInfo = [];
for i = 1:Num_CellWithNP
    jj = NPCells(i,1);
    ii = NPCells(i,2);

    x1domain = [x1(ii),x1(ii+1)];
    x2domain = [x2(jj),x2(jj+1)];

    x1x2_0 = [0.5*sum(x1domain);0.5*sum(x2domain)];  % Set start point of root finding

    % Bx
    Bx_Cell_i = [B1(jj,ii), B1(jj+1,ii);B1(jj,ii+1),B1(jj+1,ii+1)];
    [cBx,func_Bx] = infunc_LinearInterp2(Bx_Cell_i,x1domain,x2domain);
    % By
    By_Cell_i = [B2(jj,ii), B2(jj+1,ii);B2(jj,ii+1),B2(jj+1,ii+1)];
    [cBy,func_By] = infunc_LinearInterp2(By_Cell_i,x1domain,x2domain);

    [func_F,func_J] = infunc_GenRootFindingFunctions(cBx,cBy);

    icount = 0;
    x1x2 = x1x2_0;
    while true
        RF_F = func_F(x1x2);
        normF = norm(RF_F);
        RF_J = func_J(x1x2);

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

        deltax1x2 = RF_J\(-RF_F);
        
        if normF < BThres
            RFLabel = 0;
            break;
        elseif icount > NiterMax
            RFLabel = 1;
            break;
        else
            x1x2 = x1x2 + deltax1x2;
            icount = icount+1;
        end
    end

    if infunc_indomain(x1x2,x1domain,x2domain) % Identify SP type
        if RFLabel == 0
            M = func_J(x1x2); % Jacobian of B
%             for i_DB = 1:9
%                 M_intp(i_DB) = interp2(X1,X2,DB{i_DB},x1x2(1),x1x2(2),'linear');                
%             end
%             M_intp = reshape(M_intp,[3,3]);
%             M = M_intp(1:2,1:2);

            lambda = eig(M);
            Tr = trace(M);

            if Is3D  % 3D 
                if 0 == Tr % 2D topology in 3D
                    if sum(lambda==0) == 2
                        SPType = 1;
                    else
                        if sum(infunc_isreal(lambda,SMALL_NUMBER)) == 2
                            SPType = 2;
                        else
                            SPType = 3;
                        end
                    end
                else % 3D topology
                    if sum(infunc_isreal(lambda,SMALL_NUMBER)) == 2
                        if prod(real(lambda)) == 0
                            SPType = 4;
                        elseif prod(real(lambda)) < 0
                            SPType = 5;
                        elseif sum(real(lambda)>0) == 2
                            SPType = 6;
                        else
                            SPType = 7;
                        end
                    else
                        if Tr>0
                            SPType = 8;
                        else
                            SPType = 9;
                        end
                    end

                end
            else % 2D
                if sum(lambda==0) == 2
                    SPType = 1;
                else
                    if sum(infunc_isreal(lambda,SMALL_NUMBER)) == 2
                        SPType = 2;
                    else
                        SPType = 3;
                    end
                end
            end
%             SPInfo = [SPInfo; ...
%                       x1x2(1),x1x2(2),SPType, ...
%                       Tr, real(lambda(1)), imag(lambda(1)), real(lambda(2)), imag(lambda(2))];
            SPInfo = [SPInfo; ...
                      x1x2(1),x1x2(2),SPType, ...
                      M(1,1), M(1,2), M(2,1),M(2,2)];
        end
    end
end
end

function [F,J] = infunc_GenRootFindingFunctions(c1,c2)
F = @(x)  [c1(1)*x(1)*x(2) + c1(2)*x(1) + c1(3)*x(2) + c1(4);...
    c2(1)*x(1)*x(2) + c2(2)*x(1) + c2(3)*x(2) + c2(4)];

J = @(x) [c1(1)*x(2)+c1(2),c1(1)*x(1)+c1(3);...
    c2(1)*x(2)+c2(2),c2(1)*x(1)+c2(3)];
end

function indomain = infunc_indomain(X,xlim,ylim)
indomain = X(1)<=xlim(2) && X(1)>=xlim(1) && X(2)<=ylim(2) && X(2)>=ylim(1); % NP in domain
end

function [c,func] = infunc_LinearInterp2(f,xdomain,ydomain)
% f = [f11,f12; f21,f22]; 
% f11 = f(iy,ix), f12 = f(iy+1,ix), f21 = f(iy,ix+1), f22 = f(iy+1,ix+1)
% f(x,y) = a*x*y + b*x + c*y + d

f11 = f(1,1);
f12 = f(1,2);
f21 = f(2,1);
f22 = f(2,2);

x1 = xdomain(1);
x2 = xdomain(2);
y1 = ydomain(1);
y2 = ydomain(2);

lx = x2 - x1;
ly = y2 - y1;
lxly = lx*ly;

c = zeros(1,4);
c(1) = (f11+f22-f12-f21)/lxly;
c(2) = ((f12-f22)*y1+(f21-f11)*y2)/lxly;
c(3) = ((f21-f22)*x1+(f12-f11)*x2)/lxly;
c(4) = (f22*x1*y1+f11*x2*y2-f12*x2*y1-f21*x1*y2)/lxly;

func = @(x)  c(1)*x(1)*x(2) + c(2)*x(1) + c(3)*x(2) + c(4);
end

function b = infunc_isreal(lambda,SMALL_NUMBER)
    b = imag(lambda) < SMALL_NUMBER;
end