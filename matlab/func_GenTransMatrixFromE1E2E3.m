function [T,Tinv] = func_GenTransMatrixFromE1E2E3(e1,e2,e3)
    for i = 1:3
        for j = 1:3
            eval(sprintf('T.e%d%d = e%d.e%d;',i,j,i,j));
        end
    end
    Tinv = func_m_MatrixTranspose(T);
end