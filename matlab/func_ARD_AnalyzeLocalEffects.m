function [DNparaL,GammaL] = func_ARD_AnalyzeLocalEffects(B,N,DNpara,De1,De2,De3,DNperpOverB,T,Tinv)
Bnorm = func_m_VectorNorm(B);
e_size = size(N.e1);
e1L = func_m_VectorGenUniform([1,0,0],e_size);
e2L = func_m_VectorGenUniform([0,1,0],e_size);
% e3L = func_m_VectorGenUniform([0,0,1],e_size);

% Tranform to local frame
% DNparaL
DNparaL = func_m_VectorTransform(T,DNpara); % Checked

% GammaL
NL = func_m_VectorTransform(T,N);
DNperpOverBL = func_m_MatrixProd(func_m_MatrixProd(T,DNperpOverB),Tinv);
De1L = func_m_MatrixProd(func_m_MatrixProd(T,De1),Tinv);
De2L = func_m_MatrixProd(func_m_MatrixProd(T,De2),Tinv);
De3L = func_m_MatrixProd(func_m_MatrixProd(T,De3),Tinv);

D3e1.e1 = De1L.e13;
D3e1.e2 = De1L.e23;
D3e1.e3 = De1L.e33;

D3e2.e1 = De2L.e13;
D3e2.e2 = De2L.e23;
D3e2.e3 = De2L.e33;

% Calculate GammaPerp
GammaL_T1 = func_m_VectorSum(func_m_VectorScaleProd(e2L,  DNperpOverBL.e13),...
                             func_m_VectorScaleProd(e1L, -DNperpOverBL.e23));
GammaL_T2 = func_m_VectorSum(func_m_VectorScaleProd(D3e2, NL.e1./Bnorm),...
                             func_m_VectorScaleProd(D3e1,-NL.e2./Bnorm));
GammaL_T3 = func_m_VectorTransform(De3L,...
            func_m_VectorSum(func_m_VectorScaleProd(e1L, NL.e2./Bnorm),...
                             func_m_VectorScaleProd(e2L,-NL.e1./Bnorm))); 
GammaL = func_m_VectorSum(GammaL_T1,func_m_VectorSum(GammaL_T2,GammaL_T3));

% GammaPerp_T1 = func_m_VectorCross(e3L,func_m_VectorSum(func_m_VectorScaleProd(e1L,DNperpOverBL.e13),...
%                                 func_m_VectorScaleProd(e2L,DNperpOverBL.e23)));
% GammaPerp_T2 = func_m_VectorCross(e3L,func_m_VectorSum(func_m_VectorScaleProd(D3e1,NL.e1./Bnorm),...
%                                 func_m_VectorScaleProd(D3e2,NL.e2./Bnorm)));
% GammaPerp_T3 = func_m_VectorTransform(De3L,...
%                func_m_VectorSum(func_m_VectorScaleProd(e1L,NL.e2./Bnorm),...
%                func_m_VectorScaleProd(e2L,-NL.e1./Bnorm))); 
% GammaPerp_V = func_m_VectorSum(GammaPerp_T1,func_m_VectorSum(GammaPerp_T2,GammaPerp_T3));
% % Note: Theoretically GammaPerp.e3=0, but numerically this term can be
% % finite resulting from the error in calculating \nabla\hat{b}
end