function NPInfo = func_ANP_IdentifyNullPointType(NPInfo_Raw,FixTrace)
% NPType:
% TrM
% DetM: <0, >0, ==0
% J: ~=0, == 0
% DeltaM: >0, <0, ==0
% Jperp: ~=0, ==0
% Jpara

nNPs = size(NPInfo_Raw,1);
NPInfo = NPInfo_Raw;
ExtraInfo = zeros(nNPs,6); % TrM, DetM, J, DeltaM, Jperp, Jpara

for i = 1:nNPs
    M = reshape(NPInfo_Raw(i,5:end),[3,3]);
    
    % Fix trace
    if FixTrace == 1
        TrM0 = trace(M);
        M = M- 1/3*TrM0*diag([1,1,1]);
    end

    A = 0.5*(M-M');
    J = [A(3,2),A(1,3),A(2,1)]'*2;

    Jnorm = norm(J);
    DetM = det(M);
    TrM = trace(M);
    Q = -(M(1,1)^2 + M(2,2)^2 + M(1,1)*M(2,2) + ...
        M(1,2)*M(2,1) + M(1,3)*M(3,1) + M(2,3)*M(3,2));
    R = -DetM;
    DeltaM = -(4*Q^3+27*R^2); % Discriminant

    ExtraInfo(i,1) = TrM;
    ExtraInfo(i,2) = DetM;
    ExtraInfo(i,3) = Jnorm;
    ExtraInfo(i,4) = DeltaM;

    if ~func_iszero(DetM) % 3D cases
        if DetM < 0 % 3D positive null
            if ~func_iszero(Jnorm) % non-potential field
                % Determine direction of spine: e3
                e3 = func_DirectionSpine(M);
                Jpara = sum(J.*e3);
                Jperp = sqrt(Jnorm.^2-Jpara.^2);
                ExtraInfo(i,5) = Jperp;
                ExtraInfo(i,6) = Jpara;

                if ~func_iszero(DeltaM)
                    if DeltaM>0
                        if ~func_iszero(Jperp)
                            NPType = 1;
                        else
                            NPType = 2;
                        end
                    else
                        if ~func_iszero(Jperp)
                            NPType = 3;
                        else
                            NPType = 4;
                        end
                    end
                else
                    if ~func_iszero(Jperp)
                        NPType = 5;
                    else
                        NPType = 6;
                    end
                end
            else % potiential
                ExtraInfo(i,5) = 0;
                ExtraInfo(i,6) = 0;

                if ~func_iszero(DeltaM)
                    if DeltaM>0
                        NPType = 7;
                    else
                        NPType = -2; % Imposible in theory
                    end
                else
                    NPType = 8;
                end
            end
        else % 3D negative null
            if ~func_iszero(Jnorm) % non-potential field
                % Determine direction of spine: e3
                e3 = func_DirectionSpine(M);
                Jpara = sum(J.*e3);
                Jperp = sqrt(Jnorm.^2-Jpara.^2);
                ExtraInfo(i,5) = Jperp;
                ExtraInfo(i,6) = Jpara;

                if ~func_iszero(DeltaM)
                    if DeltaM>0
                        if ~func_iszero(Jperp)
                            NPType = 9;
                        else
                            NPType = 10;
                        end
                    else
                        if ~func_iszero(Jperp)
                            NPType = 11;
                        else
                            NPType = 12;
                        end
                    end
                else
                    if ~func_iszero(Jperp)
                        NPType = 13;
                    else
                        NPType = 14;
                    end
                end
            else % potiential
                ExtraInfo(i,5) = 0;
                ExtraInfo(i,6) = 0;

                if ~func_iszero(DeltaM)
                    if DeltaM>0
                        NPType = 15;
                    else
                        NPType = -2; % Imposible in theory
                    end
                else
                    NPType = 16;
                end
            end
        end
    else % 2D Cases
        if ~func_iszero(Jnorm) % non-potential field
            % Determine direction of spine: e3
            e3 = func_DirectionSpine(M);
            Jpara = sum(J.*e3);
            Jperp = sqrt(Jnorm.^2-Jpara.^2);
            ExtraInfo(i,5) = Jperp;
            ExtraInfo(i,6) = Jpara;

            if ~func_iszero(DeltaM)
                if DeltaM>0
                    if ~func_iszero(Jperp)
                        NPType = 17;
                    else
                        NPType = 18;
                    end
                else
                    if ~func_iszero(Jperp)
                        NPType = 19;
                    else
                        NPType = 20;
                    end
                end
            else
                if ~func_iszero(Jperp)
                    NPType = 21;
                else
                    NPType = 22;
                end
            end
        else % potiential
            ExtraInfo(i,5) = 0;
            ExtraInfo(i,6) = 0;

            if DeltaM>0
                NPType = 23;
            else
                NPType = -3; % Imposible in theory
            end
        end

    end
    NPInfo(i,4) = NPType;
end
NPInfo = [NPInfo ExtraInfo];
end
