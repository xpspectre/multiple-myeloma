%% Disease progression prediction using timeseries data
clear; close all; clc
rng('default'); % for consistency

assembled_pred_prob_file = 'data/processed/timeseries_assembled_pred_prob.mat';
baseline_file = 'data/processed/imputed_h_and_w_update.csv';
treatments_file = 'data/processed/patient_treat.csv'; % baseline treatments includes SCT

% Progression mapping strategy
prog_map = 'orig'; % original progressive disease = 1, sCR = 6
% prog_map = 'orig-centered'; % original, but stable disease = 0, positive outcomes shifted down, negative outcome made more negative
% prog_map = 'simple'; % stable disease = 0, good = +1, worse = -1

% Test-validate-train split
% tt_split = 'by-obs';
tt_split = 'by-patient'; % fairer than obs?

% Load timeseries data
loaded = load(assembled_pred_prob_file);
data = loaded.t_data;
present = loaded.t_present;

% Clean up timeseries data
%   Map progression as desired
%   Note that interpolated data won't be the clean original 1-6
% In all cases, make neutral = 0, bad = negative, good = positive
n = height(data);
prog_cols = {'PROGRESSION', 'PROGRESSION_CURR'};
for i = 1:2
    prog_col = prog_cols{i};
    switch prog_map
        case 'orig' % mak
            prog = round(data.(prog_col));
            prog(prog < 1) = 1;
            prog(prog > 6) = 6;
            prog = prog - 2;
            prog_val = -1:4;
            prog_names = {'Progressive Disease', 'Stable Disease', ...
                'Partial Response', 'Very Good Partial Response (VGPR)', ...
                'Complete Response', 'Stringent Complete Response (sCR)'}; % too long for yticks
            prog_names = {'Progressive', 'Stable', ...
                'Partial Resp', 'VGPR', ...
                'CR', 'sCR'};
        case 'orig-centered'
            prog = round(data.(prog_col));
            prog(prog < 1) = -1;
            prog(prog > 6) = 6;
            prog = prog - 2; % stable disease = 2, shift to 0
            prog(prog < 0) = -3; % extra weight to progressive disease, tunable
        case 'simple'
            prog = round(data.(prog_col));
            prog(prog < 2) = -1;
            prog(prog == 2) = 0;
            prog(prog > 2) = 1; % this is same to ordering of overwriting ops
    end
    data.(prog_col) = prog;
end

% Quick and dirty validate prog vals
% unique(prog)
% figure
% histogram(prog)

%% Load baseline data and keep age, demographics, that are unchanging but important
base = readtable(baseline_file);
base.Var1 = [];
base_keep = {'PUBLIC_ID', 'D_PT_age', 'demog_height', 'demog_weight', 'DEMOG_GENDER', ...
    'DEMOG_AMERICANINDIA', 'DEMOG_BLACKORAFRICA', 'DEMOG_WHITE', 'DEMOG_ASIAN', ...
    'DEMOG_ETHNICITY'};
base = base(:,base_keep);

% Join baseline to time data, repeating baseline data for each patient
data = join(data, base, 'Keys', 'PUBLIC_ID');

%% Load baseline treatment data and keep just the 1st line stem cell therapy status
treat = readtable(treatments_file);
treat_keep = {'PUBLIC_ID', 'line1sct'};
treat = treat(:,treat_keep);

data = join(data, treat, 'Keys', 'PUBLIC_ID'); % datapoint present for everyone

%% Display some patient profiles
% Collect data
patient = 'MMRF_1011';
patient_rows = strcmp(patient, data.PUBLIC_ID);
x = data(patient_rows,:);
x_present = present(patient_rows,:);
start = min(x.INTERVAL);
stop = max(x.INTERVAL);
intervals = x.INTERVAL;

lab_cols = {'D_LAB_cbc_abs_neut', 'D_LAB_chem_albumin', 'D_LAB_chem_bun', ...
    'D_LAB_chem_calcium', 'D_LAB_chem_creatinine', 'D_LAB_chem_glucose', ...
    'D_LAB_cbc_hemoglobin', 'D_LAB_chem_ldh', ...
    'D_LAB_cbc_platelet', 'D_LAB_chem_totprot', 'D_LAB_cbc_wbc'};
lab_names = cellfun(@(x) x(7:end), lab_cols, 'UniformOutput', false);
lab_vals = table2array(x(:,lab_cols));
lab_vals_present = table2array(x_present(:,lab_cols));
lab_vals_marks = lab_vals;
lab_vals_marks(~lab_vals_present) = NaN;

mm_lab_cols = {'D_LAB_serum_kappa', 'D_LAB_serum_m_protein', 'D_LAB_serum_iga', ...
    'D_LAB_serum_igg', 'D_LAB_serum_igm', 'D_LAB_serum_lambda'};
mm_lab_names = cellfun(@(x) x(7:end), mm_lab_cols, 'UniformOutput', false);
mm_lab_vals = table2array(x(:,mm_lab_cols));
mm_lab_vals_present = table2array(x_present(:,mm_lab_cols));
mm_lab_vals_marks = mm_lab_vals;
mm_lab_vals_marks(~mm_lab_vals_present) = NaN;

treatment_cols = {'Bortezomib', 'Carfilzomib', 'Cyclophosphamide', 'Dexamethasone', ...
    'Lenalidomide', 'Melphalan', 'Other'};
n_treatments = length(treatment_cols);
treatments = table2array(x(:,treatment_cols));
tr = cell(1,n_treatments); % line segments of intervals
for i = 1:n_treatments
    tri = zeros(0,4); % [x1 x2 y1 y2]
    tri_inds = find(treatments(:,i));
    intervals_i = intervals(tri_inds);
    n_tri_inds = length(tri_inds);
    for j = 1:n_tri_inds
        interval_start = intervals_i(j);
        interval_stop = interval_start + 1;
        tri = [tri; interval_start, interval_stop, i, i];
    end
    tr{i} = tri;
end

prog = x.PROGRESSION;


% Make plots
% Subplots 1,2, and 3 are offset by half to indicate middle of interval
figure('Position', [300,300,900,900])
subplot(4,1,1)
plot(intervals - 0.5, log10(lab_vals))
hold on
ax = gca;
ax.ColorOrderIndex = 1;
plot(intervals - 0.5, log10(lab_vals_marks), '+')
hold off
xlim([start - 1, stop + 1])
ylabel('log_{10} Value')
legend(lab_names, 'Location', 'east', 'Interpreter', 'none')
legend('boxoff')
title('General Lab Measurements')

subplot(4,1,2)
plot(intervals - 0.5, log10(mm_lab_vals))
hold on
ax = gca;
ax.ColorOrderIndex = 1;
plot(intervals - 0.5, log10(mm_lab_vals_marks), '+')
hold off
xlim([start - 1, stop + 1])
ylabel('log_{10} Value')
legend(mm_lab_names, 'Location', 'east', 'Interpreter', 'none')
legend('boxoff')
title('Multiple Myeloma Lab Measurements')

subplot(4,1,3)
hold on
ax = gca;
for i = 1:n_treatments
    tri = tr{i};
    n_tri_inds = size(tri,1);
    for j = 1:n_tri_inds
        ax.ColorOrderIndex = i;
        trij = tri(j,:);
        plot(trij(1:2), trij(3:4))
    end
end
hold off
set(gca, 'YTick', 1:n_treatments)
set(gca, 'YTickLabel', treatment_cols)
xlim([start - 1, stop + 1])
ylim([0, n_treatments+1])
title('Treatments')

subplot(4,1,4)
ax = gca;
plot(intervals - 0.5, prog, '+-')
set(gca, 'YTick', prog_val)
set(gca, 'YTickLabel', prog_names)
ylim([min(prog_val) - 1, max(prog_val)+1])
ylabel('Progression')
xlabel('Time Interval (90 days each)')
title('Disease Progression')
h = suptitle(sprintf('Patient %s Profile', patient));
set(h, 'Interpreter', 'none');

%% Preprocess data for regression
% Note/TODO: It's fairer to do the preprocessing on the
%   train data only, then apply the same transformations to the test data
% Hold onto the ID and outcome
% Keep INTERVAL as a feature - maybe there's a general timing factor
hold_out_cols = {'PUBLIC_ID','PROGRESSION'};
data_ = data(:,hold_out_cols);
data(:,hold_out_cols) = [];

% Throw out some cols that behave badly
%   The AT_* cols that are differences from last time aren't coded very well
drop_cols = {'AT_INCREASEOF25F', 'AT_SERUMMCOMPONE', 'AT_URINEMCOMPONE', ...
    'AT_ONLYINPATIENT', 'AT_ONLYINPATIENT2', 'AT_DEVELOPMENTOF'};
drop_delta_cols = strcat('delta_', drop_cols);
drop_cols = [drop_cols, drop_delta_cols, 'delta_SS_SYMPTOMATICHY']; % last one has issues
data(:,drop_cols) = [];

% Throw out any cols that are constant - single val for all rows
cols = data.Properties.VariableNames;
n_col = length(cols);
drop_cols = false(1,n_col);
for i = 1:n_col
    n = length(unique(data{:,i}));
    if n == 1
        drop_cols(i) = true;
    end
end
drop_cols = find(drop_cols);
for col = drop_cols
    fprintf('Dropping col %s because all vals are the same\n', cols{col})
end
data(:,drop_cols) = [];
cols = data.Properties.VariableNames;
n_col = length(cols);

% Toss cols that are identical to another col - keep the 1st col
%   Else the matrix will be rank deficient
%   This is messy and done as a matrix
X = table2array(data);
[~, ia, ~] = unique(X', 'rows', 'stable');
redundant_cols = setdiff(1:n_col, ia);
for col = redundant_cols
    fprintf('Dropping col %s because it is identical to another col\n', cols{col})
end
data(:,redundant_cols) = [];
cols = data.Properties.VariableNames;
n_col = length(cols);

% Toss outliers - anything beyond n std devs +/- set to n std dev
n_dev = 3;
n_outliers = zeros(1,n_col);
outlier_bounds = zeros(2,n_col);
for i = 1:n_col
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
    outlier_bounds(:,i) = [lo_lim, hi_lim];
    data{:,i} = x;
end

% Standardize all cols to between [0,1]
col_bounds = zeros(2,n_col);
for i = 1:n_col
    x = data{:,i};
    col_min = min(x);
    col_max = max(x);
    col_bounds(:,i) = [col_min, col_max];
    data{:,i} = (x - col_min) / (col_max - col_min);
end

% Re-add back in held out cols
data(:,hold_out_cols) = data_;

%% Split test-train datasets
switch tt_split
    case 'by-obs'
        n = height(data);
        [train_ind, val_ind, test_ind] = dividerand(n, 0.7 ,0.1, 0.1); % [train, val, test]
        TTSET = zeros(n,1);
        TTSET(train_ind) = 1;
        TTSET(val_ind) = 2;
        TTSET(test_ind) = 3;
        data.TTSET = TTSET;
    case 'by-patient'
        PUBLIC_ID = unique(base.PUBLIC_ID);
        n = length(PUBLIC_ID);
        [train_ind, val_ind, test_ind] = dividerand(n, 0.7 ,0.1, 0.1);
        TTSET = zeros(n,1);
        TTSET(train_ind) = 1;
        TTSET(val_ind) = 2;
        TTSET(test_ind) = 3;
        T = table(PUBLIC_ID, TTSET);
        data = join(data, T, 'Keys', 'PUBLIC_ID');
end

data_train = data(data.TTSET == 1,:);
data_test = data(data.TTSET == 3,:) ;

% Turn into features X and outputs y
drop_cols = {'PUBLIC_ID', 'TTSET'};
X_train = data_train;
X_train(:,drop_cols) = [];
% X_train = table2array(X_train);
% y_train = data_train.PROGRESSION;

X_test = data_test;
X_test(:,drop_cols) = [];
% X_test = table2array(X_test);
% y_test = data_test.PROGRESSION;

col_names = data_train.Properties.VariableNames;
% col_names = ['intercept', col_names];

% Some diagnostic plots
% Distribution of progression - not normalize
% figure
% histogram(y_train)
% hold on
% histogram(y_test)
% hold off
% legend('Train','Test')

% Run some regressions
%   glmfit defaults to adding a y-intercept term as the 1st feature automatically
% [beta, dev, stats] = glmfit(X_train, y_train);
mdl = fitglm(X_train, 'ResponseVar', 'PROGRESSION');
% mdl = stepwiseglm(X_train, y_train, 'interactions', 'VarNames', col_names); % takes forever, unprincipled


% Display sorted significant params
sortrows(mdl.Coefficients, 'pValue')

% pvals = stats.p;
% pval_cutoff = 0.1;
% keep = pvals <= pval_cutoff;
% beta = beta(keep);
% pvals = pvals(keep);
% col_names = col_names(keep);
% [~, sort_ind] = sort(abs(beta), 'descend');
% fprintf('Beta\tpval\tName\n')
% for i = 1:length(beta)
%     ind = sort_ind(i);
%     fprintf('%8.3f\t%5.4f\t%s\n', beta(ind), pvals(ind), col_names{ind})
% end

% Analysis
% Higher PROGRESSION = better => positive beta = predicts good result
% Higher ECOG_PERFORMANCEST means worse results (0 = normal, 5 = dead)
% !!! Naive analysis suggests SCT helps, but Cyclophosphamide and Dexamethasone hurt !!!
%   Possible Simpson's paradox: these drugs go to sicker patients
%   Many common drugs show this - drugs are given to sick patients (which
%   have a higher chance of dying soon) and OK patients (a while since last
%   episode) are more likely to be OK.
%   -> Need to regress on severity, then treatments

