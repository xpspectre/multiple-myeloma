%% Prep RNAseq expression data
% Use MMRF_CoMMpass_IA9_E74GTF_Cufflinks_Gene_FPKM.txt file
%   - Gene-based is easier to interpret/smaller
%   - FPKM is normalized, proper to interpret
% Replaces the prep_rnaseq_data.py file
% Run this after prep_mut.m, which extracts a sub-table from the per-visit
%   csv data to match patient PUBLIC_ID with sequencing IDs
% Notes for classifiers: 
%   - Can z-score each gene if using them as raw features
%   - Run 2-sample t-tests for each gene vs outcome of interest and keep to
%   top ones only
clear; close all; clc
rng('default');

load_raw = 1;

rnaseq_file = 'data/raw/MMRF_CoMMpass_IA9_E74GTF_Cufflinks_Gene_FPKM.txt';
per_visit_lim_out = 'data/processed/per_visit_lim.mat';
raw_file = 'data/processed/rnaseq_cufflinks_fpkm_raw.mat';
baseline_feature_file = 'data/processed/rnaseq_cufflinks_fpkm_baseline_features.mat';

%% Load main RNAseq file
switch load_raw
    case 2
        % Read data
        r1 = 1; % 1st row is sample name
        c1 = 2; % 1st col is gene id, 2nd col is description
        M = dlmread(rnaseq_file, '\t', r1, c1);
        [ng,np] = size(M);
        M = M'; % transpose to get proper patients x genes format
        
        % Read sample names
        f = fopen(rnaseq_file);
        line1 = fgets(f);
        headers = strsplit(line1, '\t');
        patients = headers(3:end);
        assert(length(patients) == np)
        
        % Read gene ids and descriptions
        C = textscan(f, '%s\t%s\t%*[^\n]');
        gene_ids = C{1};
        gene_descs = C{2};
        assert(length(gene_ids) == ng);
        
        fclose(f);
        
        save(raw_file, 'M', 'patients', 'gene_ids', 'gene_descs')
    case 1
        loaded = load(raw_file);
        M = loaded.M;
        patients = loaded.patients;
        gene_ids = loaded.gene_ids;
        gene_descs = loaded.gene_descs;
end

%% Basic dataset processing
% Delete genes that all patients have the same val (0 probably) or close;
eps = 1e-6;
[np,ng] = size(M);
toss = false(1,ng);
for ig = 1:ng
    Mi = M(:,ig);
    exp_uniq = unique(Mi);
    dev_var = std(Mi)/mean(Mi);
    if isnan(dev_var) || dev_var < eps || length(exp_uniq) == 1
        toss(ig) = true;
    end
end
M(:,toss) = [];
gene_ids(toss) = [];
gene_descs(toss) = [];
ng = length(gene_ids);
fprintf('Threw out %i genes where all patients had the same val\n', sum(toss))

%% Look at distribution of variances
Xmean = mean(M,1);
devs = std(M,0,1)./Xmean;
figure
histogram(log10(devs));
xlabel('log_{10} Std Dev/Mean')
ylabel('Count')
title('Distribution of Expression Variance by Gene')

%% Load map from per-visit patients to RNAseq patients
loaded = load(per_visit_lim_out);
vdata = loaded.vdata;
vdata = vdata(vdata{:,'VISIT'}<=0,:);
seq_ids = unique(vdata{:,'SPECTRUM_SEQ'});
seq_ids(strcmp('',seq_ids)) = [];
patient_ids = strcat(seq_ids, '_BM');
[keep_ids, ia, ib] = intersect(patient_ids, patients); % ib is the X rows to keep

M = M(ib,:);
Mmean = M(ib); % keep this for comparing patients against the pop
patients = patients(ib);
patients = regexprep(patients,'.....$',''); % delete the last 5 chars to get the PUBLIC_ID to join on (probably)
save(baseline_feature_file, 'M', 'Mmean', 'patients', 'gene_ids', 'gene_descs')

