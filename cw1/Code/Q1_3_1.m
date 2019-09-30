% Load the diffusion signal
fid = fopen('data_p1/isbi2015_data_normalised.txt', 'r', 'b');
fgetl(fid); % Read in the header
D = fscanf(fid, '%f', [6, inf])'; % Read in the data
fclose(fid);

% Select the first of the 6 voxels
meas = D(:,1);

% Load the protocol
fid = fopen('data_p1/isbi2015_protocol.txt', 'r', 'b');
fgetl(fid);
A = fscanf(fid, '%f', [7, inf]);
fclose(fid);

% Create the protocol
grad_dirs = A(1:3,:);
G = A(4,:);
delta = A(5,:);
smalldel = A(6,:);
TE = A(7,:);
GAMMA = 2.675987E8;
bvals = ((GAMMA*smalldel.*G).^2).*(delta-smalldel/3);

% Convert bvals units from s/m^2 to s/mm^2
bvals = bvals/10^6;

% Define a starting point for the non-linear fit
startx = [3.5e+0 3e-03 2.5e-01 0 0];

% Set random seed for reproducibility
rng(3);

% Define various options for the non-linear fitting algorithm.
h=optimset('MaxFunEvals',20000,...
 'Algorithm','levenberg-marquardt',...
 'TolX',1e-10,...
 'TolFun',1e-10);

t1 = cputime;

% Use inverse transformation to	maintain the same starting point
NewStartx(1) = sqrt(startx(1));
NewStartx(2) = sqrt(startx(2));
%NewStartx(3) = sqrt(-log(startx(3)));
NewStartx(3) = -log((1/startx(3))-1);
NewStartx(4) = startx(4);
NewStartx(5) = startx(5);

% Now run the fitting using NewBallStickSSD function
[parameter_hat,RESNORM,EXITFLAG,OUTPUT]=fminunc(@NewBallStickSSD,NewStartx,h,meas,bvals,grad_dirs);

t2 = cputime - t1;

% Use NewBallStickSSD function to computes predicted signal and the sum of square differences
[sumRes, S_pred] =NewBallStickSSD(parameter_hat, meas, bvals, grad_dirs);

% Mean and standard deviation of predicted signal
mean_S_pred = mean(S_pred);
std_S_pred = std(S_pred);

% Mean and standard deviation of actual signal
mean_S_actu = mean(meas);
std_S_actu = std(meas);

% Recover parameter_hat2 through inverse transformation
parameter_hat(1) = parameter_hat(1)^2; 
parameter_hat(2) = parameter_hat(2)^2;
parameter_hat(3) = (1/(1+exp(-parameter_hat(3))));

figure;
% Plot actual measurements
plot(meas, 'bs');
hold on;
% Add predictions to the plot
plot(S_pred, 'rx');
xlabel('g index');
ylabel('S');
legend('actual S', 'predicted S');

% Iteration times to find the global minimum of our objective function	
iterations = 1000;

% Create normally distributed random numbers
random1 = randn(1,iterations);
random2 = randn(1,iterations);
random3 = randn(1,iterations);
random4 = randn(1,iterations)*2*pi;
random5 = randn(1,iterations)*2*pi;

% Initalize the RESNORM array
RESNORM_array = zeros(1,iterations);

% Initialise the minimum RESNORM
min_RESNORM = 100000000000;
% Initislize the minimum parameter_hat
min_parameter_hat = zeros(1,5);

% Find the minimum RESNORM
for i=1:iterations
    %Creat a new starting point
    NewStartx(1) = startx(1)*(1+random1(i)); 
    NewStartx(2) = startx(2)*(1+random2(i));
    NewStartx(3) = startx(3)*(1+random3(i));
    NewStartx(4) = startx(4)+random4(i);
    NewStartx(5) = startx(5)+random5(i);
    % Define various options for the non-linear fitting algorithm.
    h=optimset( 'MaxFunEvals',  20000, ...
                'Algorithm' , 'levenberg-marquardt',   ...  
                'TolX' ,1e-10, 'TolFun' ,1e-10, ...
                'Display', 'off', 'LargeScale', 'off');
    % Now run the fitting
    [parameter_hat,RESNORM,EXITFLAG,OUTPUT]=fminunc(@NewBallStickSSD,NewStartx,h,meas,bvals,grad_dirs);
    RESNORM_array(i)=RESNORM;
    if RESNORM < min_RESNORM 
        min_RESNORM = RESNORM;
        min_parameter_hat=parameter_hat;
    end
end
propotion=sum(RESNORM_array < min_RESNORM + 0.0001)/iterations