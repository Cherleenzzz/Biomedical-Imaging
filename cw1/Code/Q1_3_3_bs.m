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
rng(5);

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
[parameter_hat,RESNORM,EXITFLAG,OUTPUT]=fminunc(@NewBallStickSSDIC,NewStartx,h,meas,bvals,grad_dirs);

t2 = cputime - t1;

% Use NewBallStickSSD function to computes predicted signal and the sum of square differences
[sumRes, S_pred, AIC, BIC] =NewBallStickSSDIC(parameter_hat, meas, bvals, grad_dirs);
AIC
BIC