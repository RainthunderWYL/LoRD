function [e31,e32,e33] = func_APNP_ProjectPlaneDirection(B1,B2,B3,x1,x2,x3,Parameters)
E3Type = Parameters.APNP_E3Type;

if E3Type == 0  % Default: use B Jacobian
    [e31,e32,e33] =  func_APNP_NormalV_MaxShearPlane_J(B1,B2,B3,x1,x2,x3,Parameters);
elseif E3Type == 1 % Use Symmetric part of B Jacobian
    [e31,e32,e33] =  func_APNP_NormalV_MaxShearPlane(B1,B2,B3,x1,x2,x3,Parameters);
elseif E3Type == 2 % Use local direction of B
    norm_B = sqrt(B1.^2+B2.^2+B3.^2);
    e31 = B1./norm_B;
    e32 = B2./norm_B;
    e33 = B3./norm_B;
elseif E3Type == -1 % Uniform e3 set by user
    norm_ConstE3 = norm(Parameters.APNP_ConstProjectPlaneDirection);
    e31 = ones(size(B1))*Parameters.APNP_ConstProjectPlaneDirection(1)/norm_ConstE3;
    e32 = ones(size(B1))*Parameters.APNP_ConstProjectPlaneDirection(2)/norm_ConstE3;
    e33 = ones(size(B1))*Parameters.APNP_ConstProjectPlaneDirection(3)/norm_ConstE3;
else
    error("Wrong Value of LocalFrameTensor: 0, 1, 2");
end

end