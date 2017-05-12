%% Disease progression prediction using timeseries data
clear; close all; clc
rng('default'); % for consistency

assembled_pred_prob_file = 'data/processed/timeseries_assembled_pred_prob.mat';
baseline_file = 'data/processed/imputed_update.csv';

% Progression mapping strategy
prog_map = 'orig'; % original progressive disease = 1, sCR = 6
% prog_map = 'orig-centered'; % original, but stable disease = 0, positive outcomes shifted down, negative outcome made more negative
% prog_map = 'simple'; % stable disease = 0, good = +1, worse = -1

% Test-validate-train split
tt_split = 'by-obs';
% tt_split = 'by-patient'; % fairer than obs?

% Load timeseries data
loaded = load(assembled_pred_prob_file);
data = loaded.t_data;

% Clean up timeseries data
%   Map progression as desired
%   Note that interpolated data won't be the clean original 1-6
n = height(data);
switch prog_map
    case 'orig'
        prog = round(data.PROGRESSION);
        prog(prog < 1) = 1;
        prog(prog > 6) = 6;
    case 'orig-centered'
        prog = round(data.PROGRESSION);
        prog(prog < 1) = 1;
        prog(prog > 6) = 6;
        prog = prog - 2; % stable disease = 2, shift to 0
        prog(prog < 0) = -3; % extra weight to progressive disease, tunable
    case 'simple'
        prog = round(data.PROGRESSION);
        prog(prog < 2) = -1;
        prog(prog == 2) = 0;
        prog(prog > 2) = 1; % this is same to ordering of overwriting ops
end
data.PROGRESSION = prog;
% Quick and dirty validate prog vals
% unique(prog)
% figure
% histogram(prog)

% Load baseline data and keep age, demographics, that are unchanging but important
base = readtable(baseline_file);
base.Var1 = [];
base_keep = {'PUBLIC_ID', 'D_PT_age', 'demog_height', 'demog_weight', 'DEMOG_GENDER', ...
    'DEMOG_AMERICANINDIA', 'DEMOG_BLACKORAFRICA', 'DEMOG_WHITE', 'DEMOG_ASIAN', ...
    'DEMOG_ETHNICITY'};
base = base(:,base_keep);

% Join baseline to time data, repeating baseline data for each patient
data = join(data, base, 'Keys', 'PUBLIC_ID');

% Split test-train datasets
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
X_train = data_train;
X_train.PUBLIC_ID = [];
X_train.TTSET = [];
X_train.PROGRESSION = [];
X_train = table2array(X_train);
y_train = data_train.PROGRESSION;

X_test = data_test;
X_test.PUBLIC_ID = [];
X_test.TTSET = [];
X_test.PROGRESSION = [];
X_test = table2array(X_test);
y_test = data_test.PROGRESSION;

% Some diagnostic plots
% Distribution of progression - not normalize
figure
histogram(y_train)
hold on
histogram(y_test)
hold off
legend('Train','Test')

% Run some regressions



