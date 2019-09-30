% set random seed
rng(10);
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

% perform two-sample t-test
[hyp,p1,ci,stats] = ttest2(sample1, sample2);
original_tstatistic = stats.tstat;
original_meandiff = mean(sample1)-mean(sample2);

% both groups
bothgroups = [sample1; sample2];
% construct an one-dimensional array D
D=bothgroups';

%% d(i)
total_permutations = 1000;
tstatistic_array_d = zeros(1,total_permutations);
all_permutations = zeros(total_permutations, 6);
for i=1:total_permutations
    % use randperm to generate a random set of permutations of the
    % intergers
    random_permutations = randperm(14);
    all_permutations(i,:) = random_permutations(1:6);
    bothgroups = D(random_permutations);
    group1 = bothgroups(1:6);
    group2 = bothgroups(7:14);
    [hyp,p2,ci,stats] = ttest2(group1, group2);
    tstatistic_array_d(i) = stats.tstat;        
end                                  
figure,
hist(tstatistic_array_d,50);
title('The empirical distribution of t-statistic');
p_value_d = nnz(tstatistic_array_d <=original_tstatistic)/numel(tstatistic_array_d);

%% (diii)
% sort the elements of each row of permutations
sort_permutations = sort(all_permutations(:,1:6),2);
% get the position of elments without repetition 
[~, id1, ~]=unique(sort_permutations, 'rows');

extended_permutations = 1500;
tstatistic_array_diii = zeros(1,total_permutations);
all_permutations_iii= zeros(extended_permutations, 14);
for i = 1:extended_permutations
    % randperm generates random set of integer permutations
    % integers are indices for D
    random_permutations_iii = randperm(14);
    all_permutations_iii(i,:) = random_permutations_iii;
    bothgroups_iii = D(random_permutations_iii);
    group1_iii = bothgroups_iii(1:6);
    group2_iii = bothgroups_iii(7:14);
    [hyp,p3,ci,stats] = ttest2(group1_iii, group2_iii);
    tstatistic_array_diii(i) = stats.tstat;            
end                                  
% sort the elements of each row of permutations
sort_permutations_iii = sort(all_permutations_iii(:,1:6),2);
% get the position of elments without repetition 
[~, id2, ~]=unique(sort_permutations_iii, 'rows');

new_permutations = sort_permutations_iii(id2,:);
new_permutations = new_permutations(1:total_permutations, :);
tstatistic_array_diii = tstatistic_array_diii(id2);
tstatistic_array_diii = tstatistic_array_diii(10:1010);

% p-value without repetition
p_value_diii = nnz(tstatistic_array_diii  <= original_tstatistic)/numel(tstatistic_array_diii);
