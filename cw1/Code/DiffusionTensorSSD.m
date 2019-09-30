function [sumRes, S] = DiffusionTensorSSD(x, Avox, bvals, qhat)
    S0 = x(1);
    Dxx = x(2);
    Dxy = x(3);
    Dxz = x(4);
    Dyy = x(5);
    Dyz = x(6);
    Dzz = x(7);
    D = [[Dxx Dxy Dxz];[Dxy Dyy Dyz];[Dxz Dyz Dzz]];
    % abs is a cheat to add numerical stability
    S = S0*exp(-abs(sum(qhat.*(D*qhat)).*bvals));
    sumRes = sum((Avox - S').^2);
end