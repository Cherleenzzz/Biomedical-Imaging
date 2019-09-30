function [sumRes,S,AIC,BIC] = ZeppelinStickTortSSDIC(x, Avox, bvals, qhat)
    % Extract the parameters
    S0 = x(1)^2;
    diff = x(2)^2;
    f = (1/(1+exp(-x(3))));
    theta = x(4);
    phi = x(5);   
    lambda1 = abs(x(6));
    lambda2 = (1-f)*lambda1;
    % Synthesize the signals
    fibdir = [cos(phi)*sin(theta) sin(phi)*sin(theta) cos(theta)];
    fibdotgrad = sum(qhat.*repmat(fibdir, [length(qhat) 1])');
    S = S0*(f*exp(-bvals*diff.*(fibdotgrad.^2)) + (1-f)*exp(-bvals.*(lambda2 + (lambda1-lambda2).*(fibdotgrad.^2))));
    % Compute the sum of square differences
    sumRes = sum((Avox - S').^2);
    % AIC=2N+Klog(sumRes/K)
    AIC = 12 + (3612*log(sumRes/3612));
%     AIC = 12 + (sumRes/(0.04^2));
    % BIC=NlogK+Klog(sumRes/K)
    BIC = (6*log(3612)) + (3612*log(sumRes/3612));
%     BIC = (6*log(3612)) + (sumRes/(0.04^2));
end