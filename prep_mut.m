%% Load mutation data
% Get just the non-synonomous somatic mutations from baseline samples
% Only use results from patients in the training set
clear; close all; clc
rng('default');

% Flags for different types of processing - because some of the steps are
%   expensive
load_file = 1;
slice_file = 1;
get_counts = 2;
get_baseline_features = 2;
load_per_visit_file = 2;

test_train_split_file = 'data/processed/test_train_split.mat';
mut_file = 'data/raw/MMRF_CoMMpass_IA9_NS.mut';
raw_out = 'data/processed/ns_mut_raw.mat';
lim_out = 'data/processed/ns_mut_lim.mat';
counts_file = 'data/processed/ns_mut_counts.mat';
baseline_feature_file = 'data/processed/ns_mut_baseline_features.mat'; % feature matrix ready to join to other baseline feature matrices

%% Load train-test split
loaded = load(test_train_split_file);
train_patients = loaded.patients(loaded.train_inds);

%% Load text file and convert to convenient form
switch load_file
    case 2
        T = readtable(mut_file, 'FileType', 'text', 'Delimiter', '\t');
        % Note: some col names are changed to be valid Matlab identifiers
        %   Col prepended with 'x'. Then 1st letter of name is Uppercased if
        %   it's a letter.
        % Note: empty fields in the mut file are filled in witha dot '.' This
        %   will make some cols strings. Handle this later
        save(raw_out, 'T');
    case 1
        loaded = load(raw_out);
        T = loaded.T;
    case 0
        % do nothing
end

%% Extract cols of interest
switch slice_file
    case 2
        data = T(:,{'sample', 'GENEID', 'gene', 'type', 'BIOTYPE'});
        save(lim_out, 'data');
    case 1
        loaded = load(lim_out);
        data = loaded.data;
end
nm = height(data);

%% Get counts
% Number of mutations per patient
% Number of incidences of a mutation
%   Note: a patient can have multiple mutations in a single gene
switch get_counts
    case 2
        % Per-row counts
        mut_counts = containers.Map('KeyType', 'char', 'ValueType', 'double');
        gene_id_symbol_map = containers.Map('KeyType', 'char', 'ValueType', 'char'); % for translating between gene id and symbol
        for im = 1:nm
            sample = data{im,'sample'}{1};
            
            gene_id = data{im,'GENEID'}{1};
            gene_symbol = data{im,'gene'}{1};
            
            % Add all mutations to id-symbol map
            gene_id_symbol_map(gene_id) = gene_symbol;
            
            % Only take baseline training patient data
            tokens = regexp(sample, '^(MMRF_\d{4})_(\d)_.*+', 'tokens');
            sample_time = str2double(tokens{1}{2});
            if sample_time > 1
                continue
            end
            patient = tokens{1}{1};
            if ~ismember(patient, train_patients)
                continue
            end
            
            if mut_counts.isKey(gene_id)
                mut_counts(gene_id) = mut_counts(gene_id) + 1;
            else
                mut_counts(gene_id) = 1;
            end
            
            if mod(im, 1000) == 0
                fprintf('Loaded %i rows\n', im)
            end
        end
        
        % Per-patient (actually per-patient, per-visit) counts
        %   Treats multiple mutations from a patient as 1
        patient_mut_counts = containers.Map('KeyType', 'char', 'ValueType', 'double');
        mut_counts_unique = containers.Map('KeyType', 'char', 'ValueType', 'double');
        samples = unique(data{:, 'sample'});
        ns = length(samples);
        patient_muts = cell(ns,1); % list of mutations found for each patient
        for is = 1:ns
            sample = samples{is};
            
            % Keep only the baseline data
            % Get counts for train and test patients; calc stats for only train patients
            tokens = regexp(sample, '^(MMRF_\d{4})_(\d)_.*+', 'tokens');
            sample_time = str2double(tokens{1}{2});
            if sample_time > 1
                continue
            end
            patient = tokens{1}{1};
            
            data_ip = data(strcmp(data{:,'sample'}, sample),:);
            patient_mut_counts(patient) = height(data_ip);
            muts_unique = unique(data_ip{:,'GENEID'});
            patient_muts{is} = muts_unique;
            
            % Count unique mutations only for train patients
            if ~ismember(patient, train_patients)
                continue
            end
            nmu = length(muts_unique);
            for imu = 1:nmu
                mut = muts_unique{imu};
                if mut_counts_unique.isKey(mut)
                    mut_counts_unique(mut) = mut_counts_unique(mut) + 1;
                else
                    mut_counts_unique(mut) = 1;
                end
            end
        end
        patient_muts(cellfun('isempty', patient_muts)) = []; % remove entries corresponding to non-baseline times
        
        save(counts_file, 'patient_mut_counts', 'mut_counts', 'mut_counts_unique', 'patient_muts', 'gene_id_symbol_map');
    case 1
        loaded = load(counts_file);
        patient_mut_counts = loaded.patient_mut_counts;
        mut_counts = loaded.mut_counts;
        mut_counts_unique = loaded.mut_counts_unique;
        patient_muts = loaded.patient_muts;
        gene_id_symbol_map = loaded.gene_id_symbol_map;
end

%% View counts
figure
histogram(cell2mat(values(patient_mut_counts)))
xlabel('Number of Mutations')
ylabel('Patient Count')
title('Distribution of Number of Mutations per Patient')

mut_counts_log = log10(cell2mat(values(mut_counts)));
figure
histogram(mut_counts_log)
xlabel('log_{10} Number of Incidences')
ylabel('Count')
title('Distribution of Number of Patients with Mutation')

% only count each mutation once per patient
mut_counts_unique_log = log10(cell2mat(values(mut_counts_unique)));
figure
histogram(mut_counts_unique_log)
xlabel('log_{10} Number of Incidences')
ylabel('Count')
title('Distribution of Number of Patients with Mutation Present')

% Many mutations are only seen once
fprintf('%i unique genes mutated\n', length(mut_counts_unique_log))

%% Make a feature matrix from the top N mutated genes (or genes where > M patients have it)
% Top mutated genes are based only on the training patients
genes = keys(mut_counts_unique)';
counts = cell2mat(values(mut_counts_unique))';
[counts, sort_inds] = sort(counts, 'descend');
genes = genes(sort_inds);

% Threshold by mutation incidence
M = 70; % was 100, but multiplied by train_frac
keep = counts >= M;
counts = counts(keep);
genes = genes(keep);

% Make gene-pos map for faster lookup in next step
%   Note that map returns keys in alphabetical order - but we're going to
%   keep them in mutation incidence order in the feature matrix
N = length(genes);
genes_map = containers.Map(genes, 1:N);

% Build a feature matrix of the kept mutations
% Keep multiple mutation counts for now. Condense to presence/absence later
%   if desired.
samples = sort(unique(data{:, 'sample'}));
ns = length(samples);

patients = cell(ns,1);
keep_patient = true(ns,1); % mask of baseline patients to keep

switch get_baseline_features
    case 2
        X = zeros(ns,N);
        for is = 1:ns
            sample = samples{is};
            
            % Keep only the baseline data
            % Get counts for train and test patients; calc stats for only train patients
            tokens = regexp(sample, '^(MMRF_\d{4})_(\d)_.*+', 'tokens');
            sample_time = str2double(tokens{1}{2});
            if sample_time > 1
                keep_patient(is) = false;
                continue
            end
            patient = tokens{1}{1};
            patients{is} = patient;
            
            data_ip = data(strcmp(data{:,'sample'}, sample),:);
            nm = height(data_ip);
            for im = 1:nm
                gene = data_ip{im,'GENEID'}{1};
                if genes_map.isKey(gene)
                    X(is,genes_map(gene)) = X(is,genes_map(gene)) + 1;
                end
            end
            if mod(is, 100) == 0
                fprintf('Saw %i samples\n', is)
            end
        end
        patients(~keep_patient) = []; % only keep baseline data
        X(~keep_patient,:) = []; 
        
        save(baseline_feature_file, 'X', 'patients', 'genes')
    case 1
        loaded = load(baseline_feature_file);
        X = loaded.X;
        patients = loaded.patients;
        genes = loaded.genes;
end

%% Make a list of genes mutated in each patient for GO enrichment analysis
np = length(patients);
for ip = 1:np
    patient_muts_file = sprintf('data/processed/go/ns_mut_baseline_mut_lists_%s.txt', patients{ip});
    f = fopen(patient_muts_file, 'w');
    patient_mut = patient_muts{ip};
    n = length(patient_mut);
    genes = cell(n,1);
    for i = 1:n
        genes{i} = gene_id_symbol_map(patient_mut{i});
    end
    fprintf(f, strjoin(genes, '\n'));
    fclose(f);
end

