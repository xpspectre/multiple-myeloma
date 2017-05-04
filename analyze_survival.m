% Similar to analyze_enpoints.py, but using better validated Matlab code
% Plus survival regression
% https://www.mathworks.com/help/stats/kaplan-meier-methods.html
% https://www.mathworks.com/help/stats/coxphfit.html
clear; close all; clc
rng('default');

baseline_clinical_data_file = 'data/processed/baseline_clinical_data_imputed_knn.csv';
endpoints_file = 'data/processed/baseline_clinical_endp.csv';
treatments_file = 'data/processed/patient_treat.csv';

% M = csvread(baseline_clinical_data_file, 1, 1); % skip 1st row (col headers) and 1st col: Patient IDs
data = readtable(baseline_clinical_data_file);
endp = readtable(endpoints_file);
treat = readtable(treatments_file);
n_data_col = width(data) - 1; % minus the PUBLIC_ID col
n_endp_col = width(endp) - 1;

% Join data by common col PUBLIC_ID and resplit to ensure consistency
all_data = join(data, endp);
all_data = join(all_data, treat);
all_data.Properties.RowNames = data{:,'PUBLIC_ID'};
all_data.('PUBLIC_ID') = [];

data = all_data(:,1:n_data_col);
endp = all_data(:,n_data_col+1:n_data_col+n_endp_col);
treat = all_data(:,n_data_col+n_endp_col+1:end);

% Keep treatments as data covariates
data = join(data, treat, 'Keys', 'RowNames');

col_names = data.Properties.VariableNames;

data_orig = data;

%% Drop some not-so-good cols
% Note/TODO: Lots of cols have overwhelming majority of patients of 1 class, a few (~dozen or less) of the other class
drop_cols_names = {'CMMC_REASONFORPROC', 'CMMC_REASONCODE', 'CMMC_RECDY', 'CMMC'}; % Lots of data missing - these may be dropped by Ege so check if they're present first
drop_cols_names = [drop_cols_names, {'D_CM_ANEUPLOIDYCAT', 'D_IM_CD38_PC_PERCENT', 'D_IM_CD138_PC_PERCENT', 'D_IM_CD45_PC_PERCENT', 'D_IM_CD56_PC_PERCENT', 'D_IM_CD117_PC_PERCENT'}]; % Slightly weird and not-as-informative D_IM_CD*_PERCENT measurements
drop_cols_names = [drop_cols_names, {'D_IM_CD138_DETECTED'}]; % everyone except 1 person has ~the same val; causes matrix to be ~singular
drop_cols_names = [drop_cols_names, {'D_IM_CD38_DETECTED', 'D_IM_CD45_PC_TYPICAL_DETECTED', 'D_IM_CD56_DETECTED', 'D_IM_CD13_DETECTED', 'D_IM_CD20_DETECTED', 'D_IM_CD33_DETECTED', 'D_IM_CD52_DETECTED', 'D_IM_CD117_DETECTED'}]; % similar to above
drop_cols_names = [drop_cols_names, {'sct_bresp'}]; % less treatment treatment col
for col = drop_cols_names;
    if ismember(col, col_names)
        col_name = col{1};
        fprintf('Dropping col %s because of some reason\n', col_name)
        data.(col_name) = [];
    end
end
n_data_col = width(data);
col_names = data.Properties.VariableNames;

% Preprocess data for regression
% 1st throw out any cols that are constant - single val for all rows
drop_cols = false(1,n_data_col);
for i = 1:n_data_col
    n = length(unique(data{:,i}));
    if n == 1
        drop_cols(i) = true;
    end
end
drop_cols = find(drop_cols);
for col = drop_cols
    fprintf('Dropping col %s because all vals are the same\n', col_names{col})
end
data(:,drop_cols) = [];
n_data_col = width(data);
col_names = data.Properties.VariableNames;

% Toss cols that are identical to another col - keep the 1st col
%   Else the matrix will be rank deficient
%   This is messy and done as a matrix
X = table2array(data);
[C, ia, ic] = unique(X', 'rows', 'stable');
redundant_cols = setdiff(1:n_data_col, ia);
for col = redundant_cols
    fprintf('Dropping col %s because it is identical to another col\n', col_names{col})
end
data(:,redundant_cols) = [];
n_data_col = width(data);
col_names = data.Properties.VariableNames;

% Toss outliers - anything beyond n std devs +/- set to n std dev
n_dev = 3;
n_outliers = zeros(1,n_data_col);
for i = 1:n_data_col
    x = data{:,i};
    col_mean = mean(x);
    col_std = std(x);
    lo_lim = col_mean - n_dev*col_std;
    hi_lim = col_mean + n_dev*col_std;
    lo = x < lo_lim;
    hi = x > hi_lim;
    x(lo) = lo_lim;
    x(hi) = hi_lim;
    n_outliers(i) = sum(lo) + sum(hi);
    data{:,i} = x;
end

% Standardize all cols to between [0,1]
col_bounds = zeros(2,n_data_col);
for i = 1:n_data_col
    x = data{:,i};
    col_min = min(x);
    col_max = max(x);
    col_bounds(:,i) = [col_min, col_max];
    data{:,i} = (x - col_min) / (col_max - col_min);
end

%% Turn data into matrix for regression
X = table2array(data);
n_data_col = size(X,2);
col_names = data.Properties.VariableNames;

died = endp{:,'D_PT_deathdy'};
last_observed = max(died, endp{:,'D_PT_lstalive'});
censored = isnan(died); % 1 = censored, 0 = death observed

%% Sanity check - some overall survival plots
figure
ecdf(last_observed, 'censoring', censored, 'function', 'survivor');
title('Overall Survival')

%% Cox proportional hazards regression on main feature matrix data
coxphopt = statset('coxphfit');
coxphopt.Display = 'iter'; % doesn't work?
coxphopt.MaxIter = 100; % default

[b, logL, H, stats] = coxphfit(X, last_observed, 'Censoring', censored, 'Options', coxphopt);
pvals = stats.p;
% positive value of beta means hi val = more hazard
% Harder to converge optimization w/ all features

% Sort features by most important (magnitude, pos or neg) that are significant
%   Lots of features w/ large mags also have large p-vals
p_val_thres = 0.05;
keep = pvals <= p_val_thres;

beta = b(keep);
pvals = pvals(keep);
col_names = col_names(keep);

[~, sort_ind] = sort(abs(beta), 'descend');
n_show = 20;
n_show = min(n_show, length(beta)); % if there are not enough sig features
fprintf('The %i most important significant features are:\n', n_show)
fprintf('beta\t\tp-val\tName\n')
for i = 1:n_show
    ind = sort_ind(i);
    fprintf('%8.3f\t%5.4f\t%s\n', beta(ind), pvals(ind), col_names{ind})
end

%% Look at differential survival due to some predictors
% Find 1st important significant feature where more than 50 people are in
%   each group - prevents 1 good/bad person from biasing it too much
% Note: D_LAB_chem_creatinine and D_CM_OTHER give the wrong trend? Is this
%   just what happens when you fit the full model (and you take into account everything else)?
%   TODO: Validate w/ logistic regression - survival at 1 yr
col = 'D_IM_igl';
hi = max(data_orig{:,col});
lo = min(data_orig{:,col});
% group1 = data_orig{:,col} > (hi+lo)/2;
group1 = data_orig{:,col} > mean(data_orig{:,col});
% group1 = data_orig{:,col} > median(data_orig{:,col});

[hf1, hx1, hlo1, hhi1] = ecdf(last_observed(group1), 'censoring', censored(group1), 'function', 'cumulative hazard');
[hf2, hx2, hlo2, hhi2] = ecdf(last_observed(~group1), 'censoring', censored(~group1), 'function', 'cumulative hazard');
[f1, x1, lo1, hi1] = ecdf(last_observed(group1), 'censoring', censored(group1), 'function', 'survivor');
[f2, x2, lo2, hi2] = ecdf(last_observed(~group1), 'censoring', censored(~group1), 'function', 'survivor');

% Histogram of feature distribution in the study
figure('Position', [100, 100, 1200, 400]);
subplot(1,3,1)
edges = linspace(lo, hi, 20);
histogram(data_orig{group1,col}, edges)
hold on
histogram(data_orig{~group1,col}, edges)
hold off
legend({'Hi','Lo'}, 'Location', 'best')
title('Feature Distribution')

% Cumulative hazard
subplot(1,3,2)
% Group 1 
h1 = stairs(hx1, hf1);
c = get(h1, 'Color');
hold on
stairs(hx1, hlo1, ':', 'Color', c);
stairs(hx1, hhi1, ':', 'Color', c);
% Group 2
h2 = stairs(hx2, hf2);
c = get(h2, 'Color');
stairs(hx2, hlo2, ':', 'Color', c);
stairs(hx2, hhi2, ':', 'Color', c);
hold off
xlabel('Time (days)')
ylabel('Cumulative Hazard')
legend([h1, h2], {'Hi','Lo'}, 'Location', 'northwest')
title('Cumulative Hazard')

% Survival
subplot(1,3,3)
% Group 1 
h1 = stairs(x1, f1);
c = get(h1, 'Color');
hold on
stairs(x1, lo1, ':', 'Color', c);
stairs(x1, hi1, ':', 'Color', c);
% Group 2
h2 = stairs(x2, f2);
c = get(h2, 'Color');
stairs(x2, lo2, ':', 'Color', c);
stairs(x2, hi2, ':', 'Color', c);
hold off
xlabel('Time (days)')
ylabel('Survival')
title('Survival')

th = suptitle(col);
set(th, 'Interpreter', 'none')

% Individual cox regression on col only
%   This seems to show the right result - the above must just end up really
%   strange when taking everything else into account
[b, logL, H, stats] = coxphfit(data{:,col}, last_observed, 'Censoring', censored, 'Options', coxphopt);

%%
figure
group1 = data_orig{:,'line1sct'} == 1;
[f1, x1] = ecdf(last_observed(group1), 'censoring', censored(group1), 'function', 'survivor');
[f2, x2] = ecdf(last_observed(~group1), 'censoring', censored(~group1), 'function', 'survivor');
stairs(x1, f1)
hold on
stairs(x2, f2)
hold off
xlabel('Time (days)')
ylabel('Survival')
legend('SCT', 'No SCT')


