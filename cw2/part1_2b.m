% set random seed
rng(3);
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

% perform two-sample t-test
[hyp,p,ci,stats] = ttest(sample1, sample2);

% generate design matrix 
oness = ones(s1_size,1);
zeross = zeros(s2_size,1);
identity = eye(s1_size); 
design_mat_X = [oness,oness,identity;oness,zeross,identity];
% dimension of design matrix
dim_X = rank(design_mat_X);

% perpendicular projection operator Px corresponding to any C(X)
projection_X = design_mat_X*pinv(design_mat_X'*design_mat_X)*design_mat_X';

% both groups
bothgroups = [sample1; sample2];
% use projection_X to determine projection_Y
projection_Yhat = projection_X*bothgroups;

% compute Rx = I - Px 
Rx = eye(size(bothgroups,1))- projection_X;

% determine ehat, the projection of Y into the error space
ehat = Rx*bothgroups;


% use derived formula to determine betahat
betahat = pinv(design_mat_X'*design_mat_X)*design_mat_X'*bothgroups;

% estimate the variance of the stochastic component ehat
variance = (ehat'*ehat)/24;

% covariance matrix of estimated model parameters
covariance_mat = variance*pinv(design_mat_X'*design_mat_X);
% use covariance_mat to determine std of model parameters
std_parameters = sqrt(covariance_mat(1,1));

% derive contrast vector for comparing group differences in means 
lambda = [0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]';
% reduced model X0 corresponding to the null hypothesis 
X0 = ones(length(bothgroups),1);

% determine t-statistic to test whether one group has higher mean
t= (lambda'*betahat)/(sqrt(lambda'*covariance_mat*lambda));