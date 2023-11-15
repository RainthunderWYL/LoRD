function [Snorm,S] = func_m_MatrixSymmPart3D(DB)
    S.e11 = DB.e11;
    S.e22 = DB.e22;
    S.e33 = DB.e33;

    S.e12 = 0.5*(DB.e12 + DB.e21);
    S.e13 = 0.5*(DB.e13 + DB.e31);
    S.e23 = 0.5*(DB.e23 + DB.e32);

    Snorm = sqrt(S.e11.^2+S.e22.^2+S.e33.^2+...
        2*S.e12.^2+2*S.e13.^2+2*S.e23.^2);
end