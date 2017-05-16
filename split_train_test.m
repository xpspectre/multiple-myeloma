%% Main script for assigning train-test to ... everything
clear; close all; clc
rng('default');

% Load all clinical datasets for patient ids
baseline_data_file = 'data/processed/patient_data.csv';
data_file = 'data/processed/clinical_data.csv';
baseline_data = readtable(baseline_data_file);
data = readtable(data_file);

patients = union(baseline_data.PUBLIC_ID, data.PUBLIC_ID);
patients = sort(patients); % just in case, union should sort
n = length(patients);

train_frac = 0.7;
test_frac = 1 - train_frac;

[train_inds, ~, test_inds] = dividerand(n, train_frac, 0, test_frac);
test_train_split_file = 'data/processed/test_train_split.mat';
save(test_train_split_file, 'patients', 'train_frac', 'test_frac', 'train_inds', 'test_inds')
