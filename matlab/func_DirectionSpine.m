function e3 = func_DirectionSpine(M)
% Suppose trace of M is close to zero

    idx = 1:3;
    [U,lambda] = eig(M);
    D = diag(lambda);
    if isreal(D)
        b0 = func_iszero(D);
        if sum(b0) > 0
            ie3 = idx(~b0);
            ie3 = ie3(1);
        else
            b = D>0;
            if sum(b) == 1
                ie3 = idx(b);
            elseif sum(b) == 2
                ie3 = idx(~b);
            else
                error('Trace of M is far from zero!');
            end
        end
    else
        b = imag(D) == 0;
        ie3 = idx(b);
    end
    e3 = U(:,ie3);
end