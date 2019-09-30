function [sumRes,S,AIC,BIC] = NewBallStickSSDIC(x, Avox, bvals, qhat)
% Extract the parameters
S0 = x(1)^2;
diff = x(2)^2;
%f = exp(-x(3)^2);
f = (1/(1+exp(-x(3))));
theta = x(4);
phi = x(5);  
% Synthesize the signals
fibdir = [cos(phi)*sin(theta) sin(phi)*sin(theta) cos(theta)];
fibdotgrad = sum(qhat.*repmat(fibdir, [length(qhat) 1])');
S = S0*(f*exp(-bvals*diff.*(fibdotgrad.^2)) + (1-f)*exp(-bvals*diff));
% Compute the sum of square differences
sumRes = sum((Avox - S').^2);
% AIC=2N+Klog(sumRes/K)
AIC = 10 + (3612*log(sumRes/3612));
% AIC = 10 + (sumRes/(0.04^2));
% BIC=NlogK+Klog(sumRes/K)
BIC = (5*log(3612)) + (3612*log(sumRes/3612));
% BIC = (5*log(3612)) + (sumRes/(0.04^2));
end