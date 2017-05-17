% Similar to analyze_enpoints.py, but using better validated Matlab code
% Plus survival regression
% https://www.mathworks.com/help/stats/kaplan-meier-methods.html
% https://www.mathworks.com/help/stats/coxphfit.html
function analyze_survival
clear; close all; clc
rng('default');

data = load_data();

[data, outlier_bounds, col_bounds] = preprocess_data(data);
data_orig = data;

% plot_overall_survival(data)

export_file = 'data/processed/baseline_data_all.csv';
export_data(data, export_file);

% cox_reg(data);

% cox_reg_weighted(data);

% plot_feature_results(data_orig, 'D_IM_FLOWCYT_PCT_PC_IN_BM_LOW')

log_reg(data)

end

function [data, endp, treat] = load_data()
% Load the various baseline datasets (clinical, rnaseq, mutation), basic
%   treatments, endpoints, and train/test split
use_muts = false;
% use_muts = true;
use_sig_genes = false;
% use_sig_genes = true;
use_gsea = false;
% use_gsea = true;

baseline_clinical_data_file = 'data/processed/imputed_updated_May_15.csv';
muts_file = 'data/processed/ns_mut_baseline_features.mat';
sig_gene_file = 'data/processed/rnaseq_cufflinks_fpkm_baseline_sig_gene_expr.mat';
gsea_file = 'data/processed/gsea_features_enrichment.mat';
endpoints_file = 'data/processed/patient_endp.csv';
treatments_file = 'data/processed/patient_treat.csv';
test_train_split_file = 'data/processed/test_train_split.mat';

%% Load main baseline clinical data file
data = readtable(baseline_clinical_data_file);
data.Var1 = []; % drop extraneous 1st col

%% Load NS mutations subset
if use_muts
    loaded = load(muts_file);
    X = loaded.X;
    X = double(X > 0); % treat multiple mutations as 1
    muts = array2table(X, 'VariableNames', loaded.genes);
    muts.PUBLIC_ID = loaded.patients;
    % Not all patients have this, so remove missing patients - yeah, left join
    %   doesn't automatically fill in with NaNs
    keep = ismember(data.PUBLIC_ID, muts.PUBLIC_ID);
    data = data(keep,:);
    data = join(data, muts, 'Keys', 'PUBLIC_ID');
end

%% Load significantly expressed genes
if use_sig_genes
    loaded = load(sig_gene_file);
    % Filter by q-val again - strict to keep non-regularized regression tractable
    q_cutoff = 0.005;
    keep = loaded.qs < q_cutoff;
    ids = loaded.gene_ids(keep);
    nids = length(ids);
    expr = loaded.T(:,'PUBLIC_ID');
    for i = 1:nids
        id = [ids{i} '_expr'];
        expr.(id) = loaded.T.(id);
    end
    % Not all patients have this, so remove missing patients
    keep = ismember(data.PUBLIC_ID, expr.PUBLIC_ID);
    data = data(keep,:);
    data = join(data, expr, 'Keys', 'PUBLIC_ID');
end

%% Load GSEA results
if use_gsea
    loaded = load(gsea_file);
    enriched = loaded.enriched;
    
    % Filter gene sets that have some threshold present
    gene_sets = enriched.Properties.VariableNames(2:end);
    ng = length(gene_sets);
    n = height(enriched);
    cutoff = 0.01;
    for ig = 1:ng
        col = gene_sets{ig};
        x = enriched.(col);
        if sum(x)/n < cutoff
            fprintf('Deleting gene set %s due to not enough hits\n', col)
            enriched.(col) = [];
        end
    end
    
    % Not all patients have this, so remove missing patients
    keep = ismember(data.PUBLIC_ID, enriched.PUBLIC_ID);
    data = data(keep,:);
    data = join(data, enriched, 'Keys', 'PUBLIC_ID');
end

%% Load endpoints
endp = readtable(endpoints_file);
endp = sortrows(endp, 'PUBLIC_ID');
died = endp{:,'D_PT_deathdy'};
last_observed = max(died, endp{:,'D_PT_lstalive'}); % main output
censored = isnan(died); % 1 = censored, 0 = death observed
survival = table(endp.PUBLIC_ID, last_observed, censored, endp.D_PT_iss, 'VariableNames', {'PUBLIC_ID','LAST_OBSERVED','CENSORED', 'ISS'});
data = join(data, survival, 'Keys', 'PUBLIC_ID');

% For stratified analysis, drop ISS = NaN
data = data(~isnan(data.ISS),:);

%% Load basic treatments
% 3 coarse-grained 1st line therapies + 1st line SCT
treat = readtable(treatments_file);
treats = table(treat.PUBLIC_ID, treat.TREAT_BOR, treat.TREAT_CAR, treat.TREAT_IMI, treat.line1sct, 'VariableNames', {'PUBLIC_ID','TREAT_BOR','TREAT_CAR','TREAT_IMI','TREAT_SCT'});
data = join(data, treats, 'Keys', 'PUBLIC_ID');

%% Load train/test split
loaded = load(test_train_split_file);
train = zeros(length(loaded.patients),1);
train(loaded.train_inds) = 1;
ttsplit = table(loaded.patients, train, 'VariableNames', {'PUBLIC_ID','TRAIN_SET'});
data = join(data, ttsplit, 'Keys', 'PUBLIC_ID');
end

function export_data(data, export_file)
% Export entire processed baseline dataset for analysis in other software
writetable(data, export_file);
end

function [data, outlier_bounds, col_bounds] = preprocess_data(data)
col_names = data.Properties.VariableNames;

%% First pull out significant cols that aren't features
meta_cols = {'PUBLIC_ID', 'LAST_OBSERVED', 'CENSORED', 'ISS', 'TRAIN_SET'};
meta = data(:,meta_cols);
data(:,meta_cols) = [];

%% Drop some not-so-good cols
% Note: Most of these cols should have been dropped by Ege in the imputation code
% Note/TODO: Lots of cols have overwhelming majority of patients of 1 class, a few (~dozen or less) of the other class
drop_cols_names = {'CMMC_REASONFORPROC', 'CMMC_REASONCODE', 'CMMC_RECDY', 'CMMC'}; % Lots of data missing - these may be dropped by Ege so check if they're present first
drop_cols_names = [drop_cols_names, {'D_CM_ANEUPLOIDYCAT', 'D_IM_CD38_PC_PERCENT', 'D_IM_CD138_PC_PERCENT', 'D_IM_CD45_PC_PERCENT', 'D_IM_CD56_PC_PERCENT', 'D_IM_CD117_PC_PERCENT'}]; % Slightly weird and not-as-informative D_IM_CD*_PERCENT measurements
drop_cols_names = [drop_cols_names, {'D_IM_CD138_DETECTED'}]; % everyone except 1 person has ~the same val; causes matrix to be ~singular
drop_cols_names = [drop_cols_names, {'D_IM_CD38_DETECTED', 'D_IM_CD45_PC_TYPICAL_DETECTED', 'D_IM_CD56_DETECTED', 'D_IM_CD13_DETECTED', 'D_IM_CD20_DETECTED', 'D_IM_CD33_DETECTED', 'D_IM_CD52_DETECTED', 'D_IM_CD117_DETECTED'}]; % similar to above
% drop_cols_names = [drop_cols_names, {'DEMOG_AMERICANINDIA', 'IC_PATIENTHADANO', 'D_IM_reason', 'D_TRI_CF_TRISOMIES21'}]; % badly behaved
n_drop_cols = length(drop_cols_names);
for id = 1:n_drop_cols;
    col = drop_cols_names{id};
    if ismember(col, col_names)
        fprintf('Dropping col %s because of some reason\n', col)
        data.(col) = [];
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
[~, ia, ~] = unique(X', 'rows', 'stable');
redundant_cols = setdiff(1:n_data_col, ia);
for col = redundant_cols
    fprintf('Dropping col %s because it is identical to another col\n', col_names{col})
end
data(:,redundant_cols) = [];
n_data_col = width(data);
col_names = data.Properties.VariableNames;

% Take the log of all the lab results
% Also the expression data
% May give -Inf for 0 vals. These will be ignored and tossed in the next step
%   after conversion to NaN and then thresholded to the min val.
lab_inds = find(strncmp(col_names, 'D_LAB_', 6));
% expr_inds = find(strncmp(col_names, 'ENSG00', 6));
% lab_inds = [lab_inds, expr_inds];
n_lab = length(lab_inds);
for il = 1:n_lab
    col = col_names{lab_inds(il)};
    x = data.(col);
    x(x < 0) = 0; % fix negative vals - none should be negative in regular space
    x = log10(x);
    data.(col) = x;
end

% Toss outliers - anything beyond n std devs +/- set to n std dev
% Only calculate mean+std on training data
train = logical(meta.TRAIN_SET);
n_dev = 3;
n_outliers = zeros(1,n_data_col);
outlier_bounds = zeros(2,n_data_col);
for i = 1:n_data_col
    x = data{:,i};
    x_ = x(train);
    x_(isinf(x_)) = [];
    col_mean = mean(x_);
    col_std = std(x_);
    lo_lim = col_mean - n_dev*col_std;
    hi_lim = col_mean + n_dev*col_std;
    lo = x < lo_lim;
    hi = x > hi_lim;
    x(lo) = lo_lim;
    x(hi) = hi_lim;
    n_outliers(i) = sum(lo) + sum(hi);
    outlier_bounds(:,i) = [lo_lim, hi_lim];
    data{:,i} = x;
    
    % Sanity checks - dataset must be complete and well-behaved
    assert(~any(isinf(x)))
    assert(~any(isnan(x)))
end

% Standardize all cols to between [0,1]
% Only calculate min+max on training data
col_bounds = zeros(2,n_data_col);
for i = 1:n_data_col
    x = data{:,i};
    col_min = min(x(train));
    col_max = max(x(train));
    col_bounds(:,i) = [col_min, col_max];
    data{:,i} = (x - col_min) / (col_max - col_min);
end

% Another filter to toss cols that ended up being all the same or are
%   otherwise invalid
drop_cols = false(1,n_data_col);
for i = 1:n_data_col
    x = data{:,i};
    if any(isnan(x))
        drop_cols(i) = true;
    end
    n = length(unique(x));
    if n == 1
        drop_cols(i) = true;
    end
end
drop_cols = find(drop_cols);
for col = drop_cols
    fprintf('Dropping col %s because all vals are the same or invalid vals found\n', col_names{col})
end
data(:,drop_cols) = [];

% Rejoin with metadata cols
data(:,meta_cols) = meta;
end

function [X, last_observed, censored, col_names] = data_to_mat(data)
% Turn data into matrices for regression and easy plotting
last_observed = data.LAST_OBSERVED;
censored = data.CENSORED;
meta_cols = {'PUBLIC_ID', 'LAST_OBSERVED', 'CENSORED', 'ISS', 'TRAIN_SET'};
data(:,meta_cols) = [];
X = table2array(data);
col_names = data.Properties.VariableNames;
end

function plot_overall_survival(data)
% Sanity check - some overall survival plots
[~, last_observed, censored] = data_to_mat(data);

[f, x, lo, hi] = ecdf(last_observed, 'censoring', censored, 'function', 'survivor');

figure
h = stairs(x, f, 'LineWidth', 2);
hold on
c = get(h, 'Color');
stairs(x, lo, ':', 'Color', c, 'LineWidth', 2);
stairs(x, hi, ':', 'Color', c, 'LineWidth', 2);
hold off
xlabel('Time (days)')
ylabel('Survival')
title('Overall Survival')

save2pdf('fig_overall_survival', 1);
end

function cox_reg_weighted(data)
% Reweight pop so that SCT vs non-SCT groups have the same "severity", which is
%   given by ISS stage
% This is for inference - don't train/test split
[X, last_observed, censored, col_names] = data_to_mat(data);
ISS = data.ISS;
np = height(data);

% Output
sct = data.TREAT_SCT;

has_sct = sct == 1;
% X_sct = X(has_sct);
% X_nosct = X(~has_sct);
ISS_sct = ISS(has_sct);
ISS_nosct = ISS(~has_sct);

% Make both SCT and non-SCT groups equal pop weight
counts_pop = zeros(1,3);
for i = 1:3
    counts_pop(i) = sum(ISS == i);
end
counts_sct = zeros(1,3); % ISS = [1, 2, 3]
counts_nosct = zeros(1,3);
for i = 1:3
    counts_sct(i) = sum(ISS_sct == i);
    counts_nosct(i) = sum(ISS_nosct == i);
end
% Major difference: No SCT has much more ISS 3
weights_sct = counts_pop ./ counts_sct;
weights_nosct = counts_pop ./ counts_nosct;
% So give the most extra weight to SCT ISS 3

weights = zeros(np,1);
for i = 1:3
    weights(has_sct & ISS == i) = weights_sct(i);
    weights(~has_sct & ISS == i) = weights_nosct(i);
end
% Normalize - total sum = np
weights = weights .* np ./ sum(weights);
% weights = ones(np,1); % control test - gives a bigger difference

coxphopt = statset('coxphfit');
coxphopt.Display = 'iter'; % doesn't work?
coxphopt.MaxIter = 250; % default

[b, logL, H, stats] = coxphfit(sct, last_observed, 'Censoring', censored, 'Frequency', weights, 'Options', coxphopt);
pvals = stats.p;
fprintf('SCT beta = %g, p-val = %g\n', b, pvals);

%% SCT vs noSCT, ISS reweighting
[f1, x1, lo1, hi1] = ecdf(last_observed(has_sct), 'censoring', censored(has_sct), 'frequency', weights(has_sct), 'function', 'survivor');
[f2, x2, lo2, hi2] = ecdf(last_observed(~has_sct), 'censoring', censored(~has_sct), 'frequency', weights(~has_sct), 'function', 'survivor');
figure
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
ylim([0.4, 1])
xlabel('Time (days)')
ylabel('Survival')
h = legend([h1, h2], {'SCT','No SCT'}, 'Location', 'southwest');
title(h, sprintf('Hazard Ratio = %.3f', exp(b)));
title('ISS-Reweighted Stem Cell Transplant Treatment Effect')

%% All drugs, ISS reweighting
pop_bor = [sum(data.TREAT_BOR == 0), sum(data.TREAT_BOR == 1)]./np; % 0, 1
pop_car = [sum(data.TREAT_CAR == 0), sum(data.TREAT_CAR == 1)]./np; % 0, 1 - very few get this
pop_imi = [sum(data.TREAT_IMI == 0), sum(data.TREAT_IMI == 1)]./np; % 0, 1
pop_sct = [sum(data.TREAT_SCT == 0), sum(data.TREAT_SCT == 1)]./np; % 0, 1
pop_iss = [sum(data.ISS == 1), sum(data.ISS == 2), sum(data.ISS == 3)]./np; % 1, 2, 3
all_weights = setprod(pop_bor, pop_car, pop_imi, pop_sct, pop_iss);
weights = zeros(48,5);
for i = 1:5 % each treatment weighted by the others
    inds = setdiff(1:5, i);
    weights(:,i) = prod(all_weights(:,inds), 2);
end
mask = setprod([0,1], [0,1], [0,1], [0,1], [1,2,3]);
data_mask = data{:,{'TREAT_BOR','TREAT_CAR','TREAT_IMI','TREAT_SCT','ISS'}};
data_weights = zeros(np,5);
for j = 1:np
    % Get corresponding mask row
    ind = ismember(mask, data_mask(j,:), 'rows');
    for i = 1:5
        data_weights(j,i) = weights(ind,i);
    end
end
for i = 1:5
    data_weights(:,i) = data_weights(:,i) .* np ./ sum(data_weights(:,i));
end

lw = 2;
treatments = {'Bortezomib','Carfilzomib','IMiDs','SCT'};
for i = 1:4
    [b, logL, H, stats] = coxphfit(data_mask(:,i), last_observed, 'Censoring', censored, 'Frequency', data_weights(:,i), 'Options', coxphopt);
    [f1, x1, lo1, hi1] = ecdf(last_observed(data_mask(:,i) == 1), 'censoring', censored(data_mask(:,i) == 1), 'frequency', data_weights(data_mask(:,i) == 1,i), 'function', 'survivor');
    [f2, x2, lo2, hi2] = ecdf(last_observed(data_mask(:,i) == 0), 'censoring', censored(data_mask(:,i) == 0), 'frequency', data_weights(data_mask(:,i) == 0,i), 'function', 'survivor');
    figure
    h1 = stairs(x1, f1, 'LineWidth', lw);
    c = get(h1, 'Color');
    hold on
    stairs(x1, lo1, ':', 'Color', c, 'LineWidth', lw);
    stairs(x1, hi1, ':', 'Color', c, 'LineWidth', lw);
    % Group 2
    h2 = stairs(x2, f2, 'LineWidth', lw);
    c = get(h2, 'Color');
    stairs(x2, lo2, ':', 'Color', c, 'LineWidth', lw);
    stairs(x2, hi2, ':', 'Color', c, 'LineWidth', lw);
    hold off
    ylim([0.4, 1])
    xlabel('Time (days)')
    ylabel('Survival')
    h = legend([h1, h2], {'Treated','Untreated'}, 'Location', 'southwest');
    title(h, sprintf('Hazard Ratio = %.3f', exp(b)));
    title(sprintf('ISS+Other Treatments Reweighted %s Treatment Effect', treatments{i}))
end

i = 5; % ISS
[b, logL, H, stats] = coxphfit(data_mask(:,i), last_observed, 'Censoring', censored, 'Frequency', data_weights(:,i), 'Options', coxphopt);
[f1, x1, lo1, hi1] = ecdf(last_observed(data_mask(:,i) == 1), 'censoring', censored(data_mask(:,i) == 1), 'frequency', data_weights(data_mask(:,i) == 1,i), 'function', 'survivor');
[f2, x2, lo2, hi2] = ecdf(last_observed(data_mask(:,i) == 2), 'censoring', censored(data_mask(:,i) == 2), 'frequency', data_weights(data_mask(:,i) == 2,i), 'function', 'survivor');
[f3, x3, lo3, hi3] = ecdf(last_observed(data_mask(:,i) == 3), 'censoring', censored(data_mask(:,i) == 3), 'frequency', data_weights(data_mask(:,i) == 3,i), 'function', 'survivor');
figure
h1 = stairs(x1, f1, 'LineWidth', lw);
c = get(h1, 'Color');
hold on
stairs(x1, lo1, ':', 'Color', c, 'LineWidth', lw);
stairs(x1, hi1, ':', 'Color', c, 'LineWidth', lw);
% Group 2
h2 = stairs(x2, f2, 'LineWidth', lw);
c = get(h2, 'Color');
stairs(x2, lo2, ':', 'Color', c, 'LineWidth', lw);
stairs(x2, hi2, ':', 'Color', c, 'LineWidth', lw);
% Group 3
h3 = stairs(x3, f3, 'LineWidth', lw);
c = get(h3, 'Color');
stairs(x3, lo3, ':', 'Color', c, 'LineWidth', lw);
stairs(x3, hi3, ':', 'Color', c, 'LineWidth', lw);
hold off
xlabel('Time (days)')
ylabel('Survival')
ylim([0.4, 1])
h = legend([h1, h2, h3], {'1','2','3'}, 'Location', 'southwest');
title(h, sprintf('Hazard Ratio = %.3f', exp(b)));
title('Treatments Reweighted Survival by ISS')

for i = 1:6
    save2pdf(sprintf('fig_weighted_survival%i', i), i)
end

end


function cox_reg(data)
run_fit = true; % the fit is expensive so cache it during development
fit_cache_file = 'data/processed/coxph_fit_cache.mat';

% Cox proportional hazards regression on main feature matrix data
[X, last_observed, censored, col_names] = data_to_mat(data);
train = logical(data.TRAIN_SET);

X_train = X(train,:);
y_train = last_observed(train);
c_train = censored(train);

X_test = X(~train,:);
y_test = last_observed(~train);
c_test = censored(~train);

coxphopt = statset('coxphfit');
coxphopt.Display = 'iter'; % doesn't work?
coxphopt.MaxIter = 250; % default

if run_fit
    [b, logL, H, stats] = coxphfit(X_train, y_train, 'Censoring', c_train, 'Options', coxphopt);
    save(fit_cache_file, 'b', 'logL', 'H', 'stats');
else
    load(fit_cache_file);
end
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
end

function plot_feature_results(data, feature, weights)
%% Look at differential survival due to some predictors
% Find 1st important significant feature where more than 50 people are in
%   each group - prevents 1 good/bad person from biasing it too much
% Note: D_LAB_chem_creatinine and D_CM_OTHER give the wrong trend? Is this
%   just what happens when you fit the full model (and you take into account everything else)?
%   TODO: Validate w/ logistic regression - survival at 1 yr
[~, last_observed, censored] = data_to_mat(data);

col = feature;
hi = max(data{:,col});
lo = min(data{:,col});
% group1 = data_orig{:,col} > (hi+lo)/2;
group1 = data{:,col} > mean(data{:,col});
% group1 = data_orig{:,col} > median(data_orig{:,col});

[hf1, hx1, hlo1, hhi1] = ecdf(last_observed(group1), 'censoring', censored(group1), 'function', 'cumulative hazard');
[hf2, hx2, hlo2, hhi2] = ecdf(last_observed(~group1), 'censoring', censored(~group1), 'function', 'cumulative hazard');
[f1, x1, lo1, hi1] = ecdf(last_observed(group1), 'censoring', censored(group1), 'function', 'survivor');
[f2, x2, lo2, hi2] = ecdf(last_observed(~group1), 'censoring', censored(~group1), 'function', 'survivor');

% Histogram of feature distribution in the study
figure('Position', [100, 100, 1200, 400]);
subplot(1,3,1)
edges = linspace(lo, hi, 20);
histogram(data{group1,col}, edges)
hold on
histogram(data{~group1,col}, edges)
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
[b, logL, H, stats] = coxphfit(data{:,col}, last_observed, 'Censoring', censored);
fprintf('Individual feature %s beta: %8.3f\n', col, b)
end

function log_reg(data)
% Logistic regression for comparison
%   No censoring - throw out patients whose results are unknown (censored before n-th yr)
%   Use death = 1, live = 0 so sign of coefficients have the same interpretation
%       as Cox
%   1-yr is hard because there's less dead than features; later yrs are even
%       harder because the vast majority are censored (study hasn't gone on long enough)
[X, last_observed, censored, col_names] = data_to_mat(data);

% Convert last_observed and censored results into n-th yr survival
n = 1;
n_day = n * 365;
dead = last_observed < n_day & ~censored; % not censored means death recorded
keep = (last_observed >= n_day) | dead; % patients who don't have the full measurement timeframe are only kept if dead

% true/false output
y = dead;

% Toss patients who we don't know the result (last observed < nday & censored)
%   => may have survived past n-th yr but don't know
X = X(keep,:);
y = y(keep);

% Split train/test
train = logical(data.TRAIN_SET);
train = train(keep);
X_train = X(train,:);
y_train = y(train);
X_test = X(~train,:);
y_test = y(~train);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L1-regularized logistic regression (binomial defaults to logit link)
opts = statset('UseParallel', true);
[B, FitInfo] = lassoglm(X_train, y_train, 'binomial', 'NumLambda', 10, 'CV', 5, 'PredictorNames', col_names, 'Options', opts);

lassoPlot(B, FitInfo, 'PlotType', 'CV');

lassoPlot(B, FitInfo, 'PlotType', 'Lambda', 'XScale', 'log');

1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Logistic regression
% mnrfit needs integers of 1:k
y_ = zeros(size(y));
y_(y) = 1; % dead
y_(~y) = 2; % alive - reference
[b, dev, stats] = mnrfit(X, y_);
col_names = ['Intercept', col_names]; % Add single intercept for dead term to front for mnrfit order
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stepwise logistic reg
% m = stepwiseglm(X, y, 'linear', 'Link', 'logit', 'VarNames', col_names);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pvals = stats.p;
% positive value of beta means hi val = more hazard

% Sort features by most important (magnitude, pos or neg) that are significant
%   Lots of features w/ large mags also have large p-vals
p_val_thres = 0.1;
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
end

