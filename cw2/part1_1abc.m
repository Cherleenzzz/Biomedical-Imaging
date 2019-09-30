% set random seed
rng(3)
% size of sample 1
s1_size = 25;
% size of sample 2
s2_size = 25;
% set the means to g1 and g2
mean_s1 = 1;
mean_s2 = 1.5;
% stochastic error
e1 = 0.25*randn(s1_size, 1);
e2 = 0.25*randn(s2_size, 1);
% generate datasets
sample1 = mean_s1+e1;
sample2 = mean_s2+e2;

%% 1_(a)
% computer the mean for new datasets
new_mean_s1 = mean(sample1);
new_mean_s2 = mean(sample2);
% computer the std for new datasets
new_std_s1 = std(sample1);
new_std_s2 = std(sample2);

%% 1_(b)
% perform two-sample t-test
[hyp,p,ci,stats] = ttest2(sample1, sample2);

%% 1_(c)_(i)
% generate design matrix 
oness = ones(s1_size,1);
zeross = zeros(s2_size,1);
design_mat_X = [oness,zeross;zeross,oness];
% dimension of col space of design matrix
dim_X = rank(design_mat_X);

%% 1_(c)_(ii)
% perpendicular projection operator Px corresponding to any C(X)
projection_X = design_mat_X*inv(design_mat_X'*design_mat_X)*design_mat_X';
% check idempotence
proj_id = projection_X^2; 
% check symmetric 
proj_sy = projection_X';

%% 1_(c)_(iii)
% both groups
bothgroups = [sample1; sample2];
% use projection_X to determine projection_Y
projection_Yhat = projection_X*bothgroups;

%% 1_(c)_(iv)
% compute Rx = I - Px 
Rx = eye(size(bothgroups,1))- projection_X;
% check idempotence
Rxsq = Rx^2;
% check symmetric 
Rxt = Rx';

%% 1_(c)_(v)
% determine ehat, the projection of Y into the error space
ehat = Rx*bothgroups;

%% 1_(c)_(vi)
% determine the angle between ehat and Yhat.
angle = acos(dot(ehat,projection_Yhat));

%% 1_(c)_(vii)
% use derived formula to determine betahat
betahat = inv(design_mat_X'*design_mat_X)*design_mat_X'*bothgroups;

%% 1_(c)_(vii)
% estimate the variance of the stochastic component ehat
variance = (ehat'*ehat)/48;

%% 1_(c)_(ix)
% covariance matrix of estimated model parameters
covariance_mat = variance*inv(design_mat_X'*design_mat_X);
% use covariance_mat to determine std of model parameters
std_parameters = sqrt(covariance_mat(1,1));

%% 1_(c)_(x)
% derive contrast vector for comparing group differences in means 
lambda = [1 -1]';
% reduced model X0 corresponding to the null hypothesis 
X0 = ones(length(bothgroups),1);

%% 1_(c)_(xi)
% recalculate new betahat
new_betahat = inv(X0'*X0)*X0'*bothgroups;
% determine error via e = Y-X0B0
new_errorhat = bothgroups-X0*new_betahat;
% additional error as a result of placing the constraint
additional_error = new_errorhat - ehat;
% projection matrix for new design mat X0
projection_X0 = X0*inv(X0'*X0)*X0';
% calculate m and n-k for f-statistic
m = trace(projection_X-projection_X0);
n_k = trace(eye(length(bothgroups))-projection_X);
% estimate f-statistic of comparing X0 to X
SSE1 = new_errorhat'*new_errorhat;
SSE2 = ehat'*ehat;
F = ((SSE1 - SSE2)/m)/(SSE2/n_k);

%% 1_(c)_(xii)
% determine t-statistic to test whether one group has higher mean
t= (lambda'*betahat)/(sqrt(lambda'*covariance_mat*lambda));

%% 1_(c)_(xiii)
% ground truth beta composed of ground truth means
beta_gt = [1 1.5]';

%% 1_(c)_(xiv)
% compute the projection of the ground truth deviation e into C(X)
egtd = bothgroups - design_mat_X*beta_gt;