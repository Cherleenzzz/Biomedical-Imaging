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
% Avox = dwis(:,56,39,48);

% Bootstrap iterations
iterations_boot = 5000;

% Initialize the minimum RESNORM
min_RESNORM = 10000000000;
% Initislize the minimum parameter_hat
min_parameter = zeros(1,5);

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
    [parameter_hat,RESNORM,EXITFLAG,OUTPUT]=fminunc(@NewBallStickSSD,NewStartx,h,Avox,bvals,qhat);       
    if RESNORM < min_RESNORM 
        min_RESNORM = RESNORM;
        min_parameter = parameter_hat;
    end       
end
% Recovery transformation
min_parameter(1) = min_parameter(1)^2; 
min_parameter(2) = min_parameter(2)^2;
min_parameter(3) = (1/(1+exp(-min_parameter(3))));

% Get best parameter
best_parameter = min_parameter;

% Use ballstick function to get predicted signal
S_pred = BallStick(best_parameter, bvals, qhat)';
% Calculate initial RESNORM
init_RESNORM = sum((Avox - S_pred).^2);
% constant=sqrt(init_RESNORM/K-N)
sigma = sqrt(init_RESNORM/103);

startx = best_parameter;
startx(1) = startx(1)/1000;
S0_array= zeros(1,iterations_boot);
diff_array = zeros(1,iterations_boot);
f_array= zeros(1, iterations_boot);

for i=1:iterations_boot  
    % Now run the fitting                                          
    S_new = S_pred + sigma.*randn(size(S_pred,1), 1);
    [parameter_hat,RESNORM,EXITFLAG,OUTPUT]=fminunc(@NewBallStickSSD,startx,h, S_new,bvals,qhat);    
    % Rocovery transformation
    parameter_hat(1) = parameter_hat(1)^2; 
    parameter_hat(2) = parameter_hat(2)^2;
    parameter_hat(3) = (1/(1+exp(-parameter_hat(3))));
    
    S0_array(i) = parameter_hat(1);
    diff_array(i) = parameter_hat(2);
    f_array(i) = parameter_hat(3);
end

% Calculate means and standard deviation
mean_S0 = mean(S0_array);
std_S0 = std(S0_array);
mean_diff = mean(diff_array);
std_diff = std(diff_array);
mean_f = mean(f_array);
std_f = std(f_array);

% Calculate ranges
S0_range2sig_boot = [mean_S0 - 2*std_S0, mean_S0 + 2*std_S0];
S0_range95_boot = [prctile(S0_array,2.5), prctile(S0_array,97.5)];
diff_range2sig_boot = [mean_diff - 2*std_diff, mean_diff + 2*std_diff];
diff_range95_boot = [prctile(diff_array,2.5), prctile(diff_array,97.5)];
f_range2sig_boot = [mean_f- 2*std_f, mean_f + 2*std_f];
f_range95_boot = [prctile(f_array,2.5), prctile(f_array,97.5)];

y_bound=[0,1000];

figure,
hist(S0_array, 50);
hold on
plot(S0_range2sig_boot, 0.2*[y_bound(2) ,y_bound(2)],'g-x', 'LineWidth', 0.75);
hold on
plot(S0_range95_boot, 0.3*[y_bound(2) ,y_bound(2)],'b-x', 'LineWidth', 0.75);
legend('S0 likelihood distribution','2-sigma range', '95 percent range');
figure,
hist(diff_array, 50);
hold on
plot(diff_range2sig_boot, 0.3*[y_bound(2) ,y_bound(2)],'g-x', 'LineWidth', 0.75);
hold on
plot(diff_range95_boot, 0.4*[y_bound(2) ,y_bound(2)],'b-x', 'LineWidth', 0.75);
legend('diff likelihood distribution','2-sigma range', '95 percent range');
figure,
hist(f_array, 50);
hold on
plot(f_range2sig_boot, 0.3*[y_bound(2) ,y_bound(2)],'g-x', 'LineWidth', 0.75);
hold on
plot(f_range95_boot, 0.4*[y_bound(2) ,y_bound(2)],'b-x', 'LineWidth', 0.75);
legend('f likelihood distribution','2-sigma range', '95 percent range');

% % Plot CDF
% figure,        
% cdfplot(S0_array);
% title('S0 CDF');
% figure,
% cdfplot(diff_array);
% title('diff CDF');
% figure,
% cdfplot(f_array);
% title('f CDF');
