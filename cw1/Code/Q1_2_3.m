% Load data
load('data_p1/data.mat');
dwis=double(dwis);
dwis=permute(dwis,[4,1,2,3]);

% Load gradient directions
qhat = load('data_p1/bvecs');
% Construct b_value array
bvals = 1000*sum(qhat.*qhat);
% Extract the set of 108 measurements from one single voxel
% Take one near the centre of the image volume
Avox = dwis(:,92,65,72);

% Initialize the minimum RESNORM
min_RESNORM = 10000000000;
% Initislize the minimum parameter_hat
min_parameter = zeros(1,5);
min_hessian = zeros(5,5);

% Iteration times to find minimum resnorm
iterations = 10;

% Define a starting point for the non-linear fit
startx = [3.5e+0 3e-03 2.5e-01 0 0];
% Define various options for the non-linear fitting algorithm
h=optimset( 'MaxFunEvals',  20000, ...
            'Algorithm' , 'levenberg-marquardt', ...  
            'TolX' ,1e-10, 'TolFun' ,1e-10, ...
            'Display', 'off', 'LargeScale', 'off');

% Create normally distributed random numbers
random1 = randn(1,iterations);
random2 = randn(1,iterations);
random3 = randn(1,iterations);
random4 = randn(1,iterations)*2*pi;
random5 = randn(1,iterations)*2*pi;

% Find the minimum RESNORM
for i=1:iterations            
    %Creat a new starting point
    NewStartx(1) = startx(1)*(1+random1(i)); 
    NewStartx(2) = startx(2)*(1+random2(i));
    NewStartx(3) = startx(3)*(1+random3(i));
    NewStartx(4) = startx(4)+random4(i);
    NewStartx(5) = startx(5)+random5(i);          
    % Now run the fitting
    [parameter_hat,RESNORM,EXITFLAG,OUTPUT,grad,hessian]=fminunc(@NewBallStickSSD,NewStartx,h,Avox,bvals,qhat);       
    if RESNORM < min_RESNORM 
        min_RESNORM = RESNORM;
        min_parameter = parameter_hat;
        min_hessian = hessian;
    end       
end

% Recovery transformation
min_parameter(1) = min_parameter(1)^2; 
min_parameter(2) = min_parameter(2)^2;
min_parameter(3) = (1/(1+exp(-min_parameter(3))));

% Get best parameter
best_parameter = min_parameter;

% cov== -inv(min_hessian/(-2*sigmaNoise^2));
cov = -inv(min_hessian/(-2*200^2));
sigma = sqrt(diag(cov));

% Calculate ranges
S0_range2sig_Lap = [best_parameter(1) - 2*sigma(1), best_parameter(1) + 2*sigma(1)];
diff_range2sig_Lap = [best_parameter(2) - 2*sigma(2), best_parameter(2) + 2*sigma(2)];
f_range2sig_Lap = [best_parameter(3) - 2*sigma(3), best_parameter(3) + 2*sigma(3)];

y_bound=[0,1000];

figure,
plot(S0_range2sig_Lap, 0.3*[y_bound(2) ,y_bound(2)],'g-x', 'LineWidth', 0.75);
legend('2-sigma range');
figure,
plot(diff_range2sig_Lap, 0.3*[y_bound(2) ,y_bound(2)],'g-x', 'LineWidth', 0.75);
legend('2-sigma range');
figure,
plot(f_range2sig_Lap, 0.4*[y_bound(2) ,y_bound(2)],'g-x', 'LineWidth', 0.75);
legend('2-sigma range');

% figure,
% plot(S0_range2sig_boot, 0.1*[y_bound(2) ,y_bound(2)],'m-x', 'LineWidth', 0.75);
% hold on
% plot(S0_range2sig_MCMC, 0.2*[y_bound(2) ,y_bound(2)],'b-x', 'LineWidth', 0.75);
% hold on
% plot(S0_range2sig_Lap, 0.3*[y_bound(2) ,y_bound(2)],'g-x', 'LineWidth', 0.75);
% hold on
% plot(S0_range95_boot, 0.4*[y_bound(2) ,y_bound(2)],'r-o', 'LineWidth', 0.75);
% hold on
% plot(S0_range95_MCMC, 0.5*[y_bound(2) ,y_bound(2)],'k-o', 'LineWidth', 0.75);
% legend('2-sigma-boot range','2-sigma-MCMC range','2-sigma-Lap range','95 percent-boot range','95 percent-MCMC range');
% figure,
% plot(diff_range2sig_boot, 0.1*[y_bound(2) ,y_bound(2)],'m-x', 'LineWidth', 0.75);
% hold on
% plot(diff_range2sig_MCMC, 0.2*[y_bound(2) ,y_bound(2)],'b-x', 'LineWidth', 0.75);
% hold on
% plot(diff_range2sig_Lap, 0.3*[y_bound(2) ,y_bound(2)],'g-x', 'LineWidth', 0.75);
% hold on
% plot(diff_range95_boot, 0.4*[y_bound(2) ,y_bound(2)],'r-o', 'LineWidth', 0.75);
% hold on
% plot(diff_range95_MCMC, 0.5*[y_bound(2) ,y_bound(2)],'k-o', 'LineWidth', 0.75);
% legend('2-sigma-boot range','2-sigma-MCMC range','2-sigma-Lap range','95 percent-boot range','95 percent-MCMC range');
% figure,
% plot(f_range2sig_boot, 0.1*[y_bound(2) ,y_bound(2)],'m-x', 'LineWidth', 0.75);
% hold on
% plot(f_range2sig_MCMC, 0.2*[y_bound(2) ,y_bound(2)],'b-x', 'LineWidth', 0.75);
% hold on
% plot(f_range2sig_Lap, 0.3*[y_bound(2) ,y_bound(2)],'g-x', 'LineWidth', 0.75);
% hold on
% plot(f_range95_boot, 0.4*[y_bound(2) ,y_bound(2)],'r-o', 'LineWidth', 0.75);
% hold on
% plot(f_range95_MCMC, 0.5*[y_bound(2) ,y_bound(2)],'k-o', 'LineWidth', 0.75);
% legend('2-sigma-boot range','2-sigma-MCMC range','2-sigma-Lap range','95 percent-boot range','95 percent-MCMC range');
