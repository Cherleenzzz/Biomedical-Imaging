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
%Avox = dwis(:,92,65,72);
Avox = dwis(:,83,56,92);

% Define a starting point for the non-linear fit
startx = [3.5e+0 3e-03 2.5e-01 0 0];

% Set random seed for reproducibility
rng(3);

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
    [parameter_hat,RESNORM,EXITFLAG,OUTPUT]=fminunc('NewBallStickSSD',NewStartx, h,Avox,bvals,qhat);
    RESNORM_array(i)=RESNORM;
    if RESNORM < min_RESNORM 
        min_RESNORM = RESNORM;
        min_parameter_hat=parameter_hat;
    end
end

propotion=sum(RESNORM_array < min_RESNORM + 1)/iterations

