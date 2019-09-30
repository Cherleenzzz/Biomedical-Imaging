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

% convert bvals units from s/m^2 to s/mm^2
bvals = bvals/10^6;

% Define a starting point for the non-linear fit
startx = [3.5e+0 3e-03 2.5e-01 0 0];

% Set random seed for reproducibility
rng(3);

% Iteration times to find the global minimum of our objective function	
iterations = 100;

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
min_parameter = zeros(1,5);

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

F = Fisher(min_parameter_hat, meas, grad_dirs, bvals);

function [fisher] = Fisher(parameter_hat, Avox, qhat, bvals)
    length_params = length(parameter_hat);
    length_Axox = length(Avox);
    fisher = zeros(length_params,length_params);
    derivatives = getDerivatives(parameter_hat,bvals,qhat);
    for i=1:length_Axox
        fisher = fisher+derivatives(i,:)'*derivatives(i,:);
    end
    fisher = fisher./(0.04^2);
    fisher = fisher./(parameter_hat'*parameter_hat);
end

function [derivatives] = getDerivatives(x, bvals, qhat)
    % Extract the parameters
    S0 = x(1);
    diff = x(2);
    f = x(3);
    theta = x(4);
    phi = x(5);
    
    fibdir = [cos(phi)*sin(theta) sin(phi)*sin(theta) cos(theta)];
    fibdotgrad = sum(qhat.*repmat(fibdir, [length(qhat) 1])');    
    S_i = exp(-bvals * diff .* (fibdotgrad.^2));
    S_e = exp(-bvals*diff);
    
    deriv_F = S0*(S_i-S_e);
    
    intermediate = S0*f*S_i .* (-(2*bvals*diff) .* fibdotgrad);
    ddTheta = repmat([cos(phi)*cos(theta),sin(phi)*cos(theta),-sin(theta)], [length(qhat) 1]);
    ddPhi = repmat([-sin(phi)*sin(theta),cos(phi)*sin(theta),0], [length(qhat) 1]);
    
    deriv_Phi = intermediate.*sum(ddPhi'.*qhat);
    deriv_Theta = intermediate.*sum(ddTheta'.*qhat);
    
    deriv_Diff = S0*(f*S_i.*(-bvals*diff.*(fibdotgrad.^2))) + ((1-f)*S_e.*(-bvals));
    deriv_S0 = f*S_i + (1-f)*S_e;
    
    derivatives = [deriv_S0;deriv_Diff;deriv_F;deriv_Theta;deriv_Phi]';
end