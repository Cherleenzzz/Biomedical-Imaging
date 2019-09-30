% set random seed
rng(10)
% size of sample 1
s1_size = 6;
% size of sample 2
s2_size = 8;
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
% perform two-sample t-test
[hyp,p1,ci,stats] = ttest2(sample1, sample2);
original_tstatistic = stats.tstat;
original_meandiff = mean(sample1)-mean(sample2);


%% 1_(b)_(i)
% both groups
bothgroups = [sample1; sample2];
% construct an one-dimensional array D
D=bothgroups';

%% 1_(b)_(ii)&(iii)
% construct all the valid permutations of D
C1 = combnk(D, s1_size); 
% initialise C2
C2 = zeros(length(C1),s2_size); 

% initialise array to store all t-statistics
tstatistic_array = zeros(1, length(C1)); % initialise array to store all t-statistics
for i=1:length(C1)
    % make sure there is no repeat for valid permutations
    C2(i,:) = setdiff(D, C1(i,:));
    % compute the t-statistics for all the membership permutations
    [hyp,p2,ci,stats] = ttest2(C1(i,:), C2(i,:));
    tstatistic_array(i) = stats.tstat;    
end
% construct empirical distribution of t-statistic
figure,
hist(tstatistic_array, 50);
title('The empirical distribution of t-statistic');

%% 1_(b)_(iv)
% determine the p-value by finding the percentage of the permutations with a t-statistic greater
% than and equal to that of the original labeling.
p_value_b = nnz(tstatistic_array <= original_tstatistic)/numel(tstatistic_array);

%% 1_(c)
% repeat (b) but rather than using the t-statistic, use the difference between
% the means as the test statistic
mean_diff_array = zeros(1, length(C1));
for i=1:length(C1)
   mean_diff = mean(C1(i,:)) - mean(C2(i,:));
   mean_diff_array(i) = mean_diff;
end
figure,
hist(mean_diff_array, 50);
title('The empirical distribution of the difference between means');
% p-value computed with mean difference
p_value_c = nnz(mean_diff_array <= original_meandiff)/numel(mean_diff_array);