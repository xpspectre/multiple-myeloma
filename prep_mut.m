%% Load mutation data
% Start with the non-synonomous somatic mutations
clear; close all; clc
rng('default');

% Flags for different types of processing - because some of the steps are
%   expensive
load_file = 0;
slice_file = 1;
get_counts = 1;
get_raw_features = 1;
load_per_visit_file = 1;

mut_file = 'data/raw/MMRF_CoMMpass_IA9_NS.mut';
raw_out = 'data/processed/ns_mut_raw.mat';
lim_out = 'data/processed/ns_mut_lim.mat';
counts_file = 'data/processed/ns_mut_counts.mat';
raw_features_file = 'data/processed/ns_mut_raw_features.mat'; % all patients (including multiple times), features are num of most common mutations
per_visit_lim_out = 'data/processed/per_visit_lim.mat';
baseline_feature_file = 'data/processed/ns_mut_baseline_features.mat'; % feature matrix ready to join to other baseline feature matrices

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
        for im = 1:nm
            patient = data{im,'sample'}{1};
            gene_id = data{im,'GENEID'}{1};
            
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
        patients = unique(data{:, 'sample'});
        np = length(patients);
        for ip = 1:np
            patient = patients{ip};
            data_ip = data(strcmp(data{:,'sample'},patient),:);
            patient_mut_counts(patient) = height(data_ip);
            muts_unique = unique(data_ip{:,'GENEID'});
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
        
        save(counts_file, 'patient_mut_counts', 'mut_counts', 'mut_counts_unique')
    case 1
        loaded = load(counts_file);
        patient_mut_counts = loaded.patient_mut_counts;
        mut_counts = loaded.mut_counts;
        mut_counts_unique = loaded.mut_counts_unique;
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
genes = keys(mut_counts_unique)';
counts = cell2mat(values(mut_counts_unique))';
[counts, sort_inds] = sort(counts, 'descend');
genes = genes(sort_inds);

% Threshold by mutation incidence
M = 100;
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
patients = sort(unique(data{:, 'sample'}));
np = length(patients);

switch get_raw_features
    case 2
        X = zeros(np,N);
        for ip = 1:np
            patient = patients{ip};
            data_ip = data(strcmp(data{:,'sample'},patient),:);
            nm = height(data_ip);
            for im = 1:nm
                gene = data_ip{im,'GENEID'}{1};
                if genes_map.isKey(gene)
                    X(ip,genes_map(gene)) = X(ip,genes_map(gene)) + 1;
                end
            end
            if mod(ip, 100) == 0
                fprintf('Saw %i patients\n', ip)
            end
        end
        save(raw_features_file, 'X', 'patients', 'genes')
    case 1
        loaded = load(raw_features_file);
        X = loaded.X;
end

%% Make a baseline feature vector to join to patients
switch load_per_visit_file
    case 2
        per_visit_file = 'data/raw/clinical_data_tables/CoMMpass_IA9_FlatFiles/PER_PATIENT_VISIT.csv';
        per_visit = readtable(per_visit_file);
        vdata = per_visit(:,{'PUBLIC_ID', 'VISIT', 'SPECTRUM_SEQ'});
        save(per_visit_lim_out, 'vdata');
    case 1
        loaded = load(per_visit_lim_out);
        vdata = loaded.vdata;
end
vdata = vdata(vdata{:,'VISIT'}<=0,:);
seq_ids = unique(vdata{:,'SPECTRUM_SEQ'});
seq_ids(strcmp('',seq_ids)) = [];
patient_ids = strcat(seq_ids, '_BM');
[keep_ids, ia, ib] = intersect(patient_ids, patients); % ib is the X rows to keep

X = X(ib,:);
patients = patients(ib);
patients = regexprep(patients,'.....$',''); % delete the last 5 chars to get the PUBLIC_ID to join on (probably)
save(baseline_feature_file, 'X', 'patients', 'genes')
