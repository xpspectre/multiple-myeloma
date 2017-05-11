%% Analyze per-visit data over time
% No genomic data for now (or ever: we don't have much of that during follow-up visits)
% Setup: fill in missing time data, break into 3 month chunks for analysis
clear; close all; clc
rng('default');

%% Load data
base_data_file = 'data/processed/imputed_update.csv'; % imputed dataset for baseline data from Ege, similar to baseline data from prep_baseline_clinical_data.py
time_data_file = 'data/processed/clinical_data.csv'; % main processed per-visit data from prep_clinical_data.py
time_endp_file = 'data/processed/clinical_endp.csv'; % main processed per-visit disease progression data from prep_clinical_data.py
data_file = 'data/processed/patient_data.csv'; % demographic data from prep_patient_data.py
endp_file = 'data/processed/patient_endp.csv'; % overall survival data from prep_patient_data.py
treat_file = 'data/processed/patient_treatment_matrix.mat'; % detailed treatment over time matrix from analyze_treatments.py

base_data = readtable(base_data_file);
time_data = readtable(time_data_file);
time_endp = readtable(time_endp_file);
data = readtable(data_file);
endp = readtable(endp_file);
treat = load(treat_file);

% Baseline and per-visit cols to drop
%   List of cols to drop from prev analysis
drop_cols = {'CMMC_REASONFORPROC', 'CMMC_REASONCODE', 'CMMC_RECDY', 'CMMC', ...
    'D_CM_ANEUPLOIDYCAT', 'D_IM_CD38_PC_PERCENT', 'D_IM_CD138_PC_PERCENT', ...
    'D_IM_CD45_PC_PERCENT', 'D_IM_CD56_PC_PERCENT', 'D_IM_CD117_PC_PERCENT', ...
    'D_IM_CD138_DETECTED', 'D_IM_CD38_DETECTED', 'D_IM_CD45_PC_TYPICAL_DETECTED', ...
    'D_IM_CD56_DETECTED', 'D_IM_CD13_DETECTED', 'D_IM_CD20_DETECTED', ...
    'D_IM_CD33_DETECTED', 'D_IM_CD52_DETECTED', 'D_IM_CD117_DETECTED', ...
    'sct_bresp'};
n_drop_cols = length(drop_cols);
base_data_cols = base_data.Properties.VariableNames;
time_data_cols = time_data.Properties.VariableNames;
for i = 1:n_drop_cols
    col = drop_cols{i};
    if ismember(col, base_data_cols)
        base_data.(col) = [];
    end
    if ismember(col, time_data_cols)
        time_data.(col) = [];
    end
end

%% Combine and process data
% Each patient has a set of unchanging features ("demographic data")

% Each patient has changing features ("lab data") at irregular intervals
%   Augment t = 0 with imputed baseline data

% Each patient has treatments ongoing during specified intervals

