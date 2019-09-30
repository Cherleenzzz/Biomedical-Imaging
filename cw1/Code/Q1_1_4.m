% Load data
load('data_p1/data.mat');
dwis=double(dwis);
dwis=permute(dwis,[4,1,2,3]);

% Load gradient directions
qhat = load('data_p1/bvecs');
% Construct b_value array
bvals = 1000*sum(qhat.*qhat);

% Iteration times to find global minimum of each voxel
iterations=5;

xvoxels = 174;
yvoxels = 145;

% Initialize S0 matrix
S0_mat = zeros(xvoxels,yvoxels);
% Initialize diff matrix
d_mat = zeros(xvoxels,yvoxels);
% Initialize f matrix
f_mat = zeros(xvoxels,yvoxels);
% Initialize theta matrix
theta_mat = zeros(xvoxels,yvoxels);
% Initialize phi matrix
phi_mat = zeros(xvoxels,yvoxels);
% Initialize RESNORM matrix
RESNORM_mat = zeros(xvoxels,yvoxels);

% Define various options for the non-linear fitting algorithm.
h=optimset( 'MaxFunEvals',  20000, ...
            'Algorithm' , 'levenberg-marquardt',   ...  
            'TolX' ,1e-10, 'TolFun' ,1e-10, ...
            'Display', 'off', 'LargeScale', 'off');

% Loop over each voxel
for y=1:xvoxels
    for z=1:yvoxels
        % Extract set of 108 measurements from one single voxel
        Avox = dwis(:,z,y,72);
        % Define a starting point for the non-linear fit
        startx = [3.5e+0 3e-03 2.5e-01 0 0];
        % Initialize the minimum RESNORM
        min_RESNORM = 10000000000;
        % Initislize the minimum parameter_hat
        min_parameter = zeros(1,5);
        % Create normally distributed random numbers
        random1 = randn(1,iterations);
        random2 = randn(1,iterations);
        random3 = randn(1,iterations);
        random4 = randn(1,iterations)*2*pi;
        random5 = randn(1,iterations)*2*pi;
        for i=1:iterations    
            NewStartx(1) = startx(1)*(1+random1(i)); 
            NewStartx(2) = startx(2)*(1+random2(i));
            NewStartx(3) = startx(3)*(1+random3(i));
            NewStartx(4) = startx(4)+random4(i);
            NewStartx(5) = startx(5)+random5(i);                                  
            try
            % Now run the fitting
            [parameter_hat,RESNORM,EXITFLAG,OUTPUT]=fminunc(@NewBallStickSSD,NewStartx, h, Avox,bvals,qhat);       
            catch
            NewStartx = startx;
            end
            if RESNORM < min_RESNORM 
                min_RESNORM = RESNORM;
                min_parameter = parameter_hat;
            end   
        end 
        
        % Recovery through transformation
        S0_mat(y,z) = min_parameter(1)^2;
        d_mat(y,z) = min_parameter(2)^2;
        f_mat(y,z) = (1/(1+exp(-min_parameter(3))));
        theta_mat(y,z) = min_parameter(4); 
        phi_mat(y,z) = min_parameter(5);
        RESNORM_mat(y,z) = min_RESNORM;
        z
    end
    y
end   

% % Plot paremeters
% figure,
% surf(S0_mat);
% figure,
% surf(d_mat);
% figure,
% surf(f_mat);
% % Plot RESNORM
% figure,
% surf(RESNORM_mat);

figure,
imshow(S0_mat, [min(S0_mat(:)), max(S0_mat(:))]);
title('S0 map');
figure,
imshow(d_mat, [min(d_mat(:)), max(d_mat(:))]);
title('d map');
figure,
imshow(f_mat, [min(f_mat(:)), max(f_mat(:))]);
title('f map');
figure,
imshow(RESNORM_mat, [min(RESNORM_mat(:)), max(RESNORM_mat(:))]);
title('RESNORM map');


[xfib,yfib,zfib]=sph2cart(theta_mat, phi_mat, ones(xvoxels,yvoxels));
% Plot fibre direction map
figure,
quiver(xfib,yfib,1.5);
title('Fibre direction map');
figure,
quiver(xfib.*f_mat, yfib.*f_mat, 1.5);
title('Coherent fibre direction map');
