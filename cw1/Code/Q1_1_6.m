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

% Define a starting point for the non-linear fit
startx = [3.5e+0 3e-03 2.5e-01 0 0];
% Define various options for the non-linear fitting algorithm.
h=optimset( 'MaxFunEvals',  20000, ...
            'Algorithm' , 'active-set',   ...  
            'TolX' ,1e-10, 'TolFun' ,1e-10, 'Display', 'off',...
            'LargeScale', 'off');

t1 = cputime;

% Now run the fitting
lb = [0; 0; 0; -inf; -inf];
ub = [inf; inf; 1; inf; inf];
[parameter_hat,RESNORM,EXITFLAG,OUTPUT]=fmincon(@BallStickSSD,startx,[],[],[],[],lb, ub, [], h, Avox,bvals,qhat);

t2 = cputime - t1

% Use NewBallStickSSD function to computes predicted signal and the sum of square differences
[sumRes, S_pred] =NewBallStickSSD(parameter_hat, Avox, bvals, qhat);

% Mean and standard deviation of predicted signal
mean_S_pred = mean(S_pred);
std_S_pred = std(S_pred);

% Mean and standard deviation of actual signal
mean_S_actu = mean(Avox);
std_S_actu = std(Avox);

figure;
% Plot actual measurements
plot(Avox, 'bs');
hold on;
% Add predictions to the plot
plot(S_pred, 'rx');
xlabel('g index');
ylabel('S');
legend('actual S', 'predicted S');
