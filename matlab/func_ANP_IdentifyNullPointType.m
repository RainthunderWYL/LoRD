function NPInfo = func_ANP_IdentifyNullPointType(NPInfo_Raw,FixTrace)
% NPType:
% TrM
% DetM: <0, >0, ==0
% J: ~=0, == 0
% DeltaM: <0, >0, ==0
% Jperp: ~=0, ==0
% Jpara

nNPs = size(NPInfo_Raw,1);
NPInfo = NPInfo_Raw;
ExtraInfo = zeros(nNPs,5); % DetM, J, Jperp, DM, TrM

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
    DeltaM = 4*Q^3+27*R^2; % Discriminant

    ExtraInfo(i,1) = DetM;
    ExtraInfo(i,2) = Jnorm;
    ExtraInfo(i,3) = zeros(size(DetM));
    ExtraInfo(i,4) = DeltaM;
    ExtraInfo(i,5) = TrM;

    if ~func_iszero(DetM) % 3D cases
        if ~func_iszero(Jnorm) % non-potential
            e3 = func_DirectionSpine(M);
            Jpara = sum(J.*e3);
            Jperp = sqrt(Jnorm.^2-Jpara.^2);
            ExtraInfo(i,3) = Jperp;

            if ~func_iszero(Jperp)
                if ~func_iszero(DeltaM) && DeltaM < 0
                    if DetM>0
                        NPType = 1;
                    else
                        NPType = 2;
                    end
                elseif ~func_iszero(DeltaM) && DeltaM > 0
                    if DetM>0
                        NPType = 3;
                    else
                        NPType = 4;
                    end
                else
                    if DetM>0
                        NPType = 5;
                    else
                        NPType = 6;
                    end
                end
            else
                if ~func_iszero(DeltaM) && DeltaM < 0
                    if DetM>0
                        NPType = 7;
                    else
                        NPType = 8;
                    end
                elseif ~func_iszero(DeltaM) && DeltaM > 0
                    if DetM>0
                        NPType = 9;
                    else
                        NPType = 10;
                    end
                else
                    if DetM>0
                        NPType = 11;
                    else
                        NPType = 12;
                    end
                end
            end
        else % potential
            if ~func_iszero(DeltaM) && DeltaM < 0
                if DetM>0
                    NPType = 13;
                else
                    NPType = 14;
                end
            elseif func_iszero(DeltaM)
                if DetM>0
                    NPType = 15;
                else
                    NPType = 16;
                end
            else
                NPType = -1;
            end
        end
    else % 2D
        if ~func_iszero(Jnorm) % non-potential
            e3 = func_DirectionSpine(M);
            Jpara = sum(J.*e3);
            Jperp = sqrt(Jnorm.^2-Jpara.^2);
            ExtraInfo(i,3) = Jperp;

            if ~func_iszero(Jperp)
                if ~func_iszero(DeltaM) && DeltaM < 0
                    NPType = 17;
                elseif ~func_iszero(DeltaM) && DeltaM > 0
                    NPType = 18;
                else
                    NPType = 19;
                end
            else
                if ~func_iszero(DeltaM) && DeltaM < 0
                    NPType = 20;
                elseif ~func_iszero(DeltaM) && DeltaM > 0
                    NPType = 21;
                else
                    NPType = 22;
                end
            end
        else % potential
            if ~func_iszero(DeltaM) && DeltaM < 0
                NPType = 23;
            else
                NPType = -2;
            end
        end
    end
    NPInfo(i,4) = NPType;
end
NPInfo = [NPInfo ExtraInfo];
end
