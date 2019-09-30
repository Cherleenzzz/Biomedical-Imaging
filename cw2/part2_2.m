% load files
CPA_I = [4,5,6,7,8,9,10,11]; % identifier for CPA file names
PPA_I = [3,6,9,10,13,14,15,16]; % identifier for PPA file names
number_subjects = 8;
% group1 had 8 subjects, 64000=40*40*40
CPA = zeros(64000,number_subjects);
j = 0;
for i=drange(CPA_I)   
    filename = sprintf('CPA%d_diffeo_fa.img',i);
    fid = fopen(filename, 'r', 'l'); % little-endian
    CPA_sub = fread(fid, 'float'); % 16-bit floating point
    j=j+1;
    CPA(:,j) = CPA_sub;    
end
% group2 had 8 subjects, 64000=40*40*40
PPA = zeros(64000, number_subjects);
j = 0;
for i=drange(PPA_I)   
    filename = sprintf('PPA%d_diffeo_fa.img',i);
    fid = fopen(filename, 'r', 'l'); % little-endian
    PPA_sub = fread(fid, 'float'); % 16-bit floating point
    j=j+1;
    PPA(:,j) = PPA_sub;    
end
% wm_mask.img is an additional binary volume defining the ROI for 
% statistical analysis 
fid = fopen('wm_mask.img', 'r', 'l');
mask = fread(fid, 'float');
voxels = [CPA PPA];

%% 2_(a) 
% we use the GLM in part1_1c Y = X1b1 + X2b2 + e.
% generate design matrix 
oness = ones(number_subjects,1);
zeross = zeros(number_subjects,1);
design_mat_X = [oness,zeross;zeross,oness];
% contrast vector
lambda = [1 -1]';
% determine betahat
betahat = voxels*(inv(design_mat_X'*design_mat_X)*design_mat_X')';
% determine ehat
ehat = voxels-(betahat*design_mat_X');
% covariance matrix of estimated model parameters
variance = sum((ehat.*ehat), 2)./ 14;
t_statiatic = (betahat*lambda)./sqrt(variance.*(lambda'*(inv(design_mat_X'*design_mat_X))*lambda));
% compute t-statisti for ROI
ROItstatistic = t_statiatic.*mask;
% compute maximum t-statistic. 
maxtstatistic = max(ROItstatistic);

%% 2_(b)
s1_size = 8;
s2_size = 8;
indices = 1:(s1_size+s2_size);
% construct all the valid permutations of indices
C1 = combnk(indices, s1_size);
% initialise C2
C2 = zeros(length(C1),s2_size); 
for i=1:length(C1)
    % make sure there is no repeat for valid permutations
    C2(i,:) = setdiff(indices, C1(i,:));
end
bothgroups = [C1 C2];
% initialise array to store all max t-statistic
max_tstatistic_array = zeros((length(C1)),1);
% just keep the values in ROI
voxels = voxels.*mask;
for i=1:length(C1)
    voxels_now = voxels(:,bothgroups(i,:));
    betahat_now = voxels_now*(inv(design_mat_X'*design_mat_X)*design_mat_X')';
    ehat_now = voxels_now-betahat_now*design_mat_X';
    % find std of error
    variance = sum((ehat_now.*ehat_now), 2)./ 14;
    tstatistic = (betahat_now*lambda)./sqrt(variance.*(lambda'*(inv(design_mat_X'*design_mat_X))*lambda));
    max_tstatistic_array(i,1) = max(tstatistic); 
end
figure,
hist(max_tstatistic_array,50);
title('The empirical distribution of max t-statistic');

%% (c)
p_value_d = nnz(max_tstatistic_array >=maxtstatistic)/numel(max_tstatistic_array);

%% (d)
% determine maximum t-statistic threshold corresponding to p-value of 5% 
threeshold=prctile(max_tstatistic_array, 95);

