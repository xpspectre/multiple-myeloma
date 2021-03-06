%% Prep RNAseq expression data
% Use MMRF_CoMMpass_IA9_E74GTF_Cufflinks_Gene_FPKM.txt file
%   - Gene-based is easier to interpret/smaller
%   - FPKM is normalized, proper to interpret
% Replaces the prep_rnaseq_data.py file
% Run this after prep_mut.m, which extracts a sub-table from the per-visit
%   csv data to match patient PUBLIC_ID with sequencing IDs
% Make sure the gene set file h.all.v6.0.symbols.gmt and the gene chip file
%   GENE_SYMBOL.chip are present in the data/processed/ dir.
% Notes for classifiers:
%   - Can z-score each gene if using them as raw features
%   - Run 2-sample t-tests for each gene vs outcome of interest and keep to
%   top ones only
% Only use results from patients in the training set
clear; close all; clc
rng('default');

load_raw = 2;
build_gene_table = 2;
build_dataset = 2;

test_train_split_file = 'data/processed/test_train_split.mat';
rnaseq_file = 'data/raw/MMRF_CoMMpass_IA9_E74GTF_Cufflinks_Gene_FPKM.txt';
raw_file = 'data/processed/rnaseq_cufflinks_fpkm_raw.mat';
baseline_feature_file = 'data/processed/rnaseq_cufflinks_fpkm_baseline_features.mat';
gene_ids_file = 'data/processed/rnaseq_cufflinks_fpkm_gsea_gene_ids.txt';
gene_table_in = 'gene_table.txt';
gene_table_out = 'data/processed/gene_table.mat';
missing_genes_file = 'data/processed/genes_to_lookup.txt';

% Paths for large-scale GSEA
%   Use absolute paths
base_dir = '/data/kshi/multiple_myeloma/processed/';
gsea_inputs_dir = '/data/kshi/multiple_myeloma/processed/gsea2/';
gsea_outputs_dir = '/data/kshi/multiple_myeloma/processed/gsea_outputs2/';

if ~exist(gsea_inputs_dir, 'dir')
    mkdir(gsea_inputs_dir);
end
if ~exist(gsea_outputs_dir, 'dir')
    mkdir(gsea_outputs_dir);
end

%% Load train-test split
loaded = load(test_train_split_file);
train_patients = loaded.patients(loaded.train_inds);

%% Load main RNAseq file
switch load_raw
    case 2
        % Read data
        r1 = 1; % 1st row is sample name
        c1 = 2; % 1st col is gene id, 2nd col is description
        M = dlmread(rnaseq_file, '\t', r1, c1);
        [ng,ns] = size(M);
        M = M'; % transpose to get proper patients x genes format
        
        % Read sample names
        f = fopen(rnaseq_file);
        line1 = fgets(f);
        headers = strsplit(line1, '\t');
        samples = headers(3:end);
        assert(length(samples) == ns)
        
        % Read gene ids and descriptions
        C = textscan(f, '%s\t%s\t%*[^\n]');
        gene_ids = C{1};
        gene_descs = C{2};
        assert(length(gene_ids) == ng);
        fclose(f);
        
        % Keep only the baseline BM data
        patients = cell(ns,1);
        keep_patient = true(ns,1); % mask of baseline patients to keep
        for is = 1:ns
            sample = samples{is};
            tokens = regexp(sample, '^(MMRF_\d{4})_(\d)_BM', 'tokens');
            if isempty(tokens)
                keep_patient(is) = false;
                continue
            end
            sample_time = str2double(tokens{1}{2});
            if sample_time > 1
                keep_patient(is) = false;
                continue
            end
            patients{is} = tokens{1}{1};
        end
        
        % Toss non-baseline samples
        M(~keep_patient,:) = [];
        patients(~keep_patient) = [];
        
        % Sort rows for convenience
        [patients, sort_inds] = sort(patients);
        M = M(sort_inds,:);
        
        save(raw_file, 'M', 'patients', 'gene_ids', 'gene_descs')
    case 1
        loaded = load(raw_file);
        M = loaded.M;
        patients = loaded.patients;
        gene_ids = loaded.gene_ids;
        gene_descs = loaded.gene_descs;
end

%% Build map of Ensembl Gene IDs to Gene Symbols
% The gene_table.txt file is from:
% General tool URL: http://www.ensembl.org/biomart/martview/8bb3d198e9ec47bed36e5bf1dc6c9019
%   Dataset: Human genes (GRCh38.p10)
%   Filters: Gene stable ID(s) [e.g. ENSG00000000003]: [ID-list specified]
%   Attributes: Gene stable ID, Gene name, Gene description
%   No Other Datasets
% Note: Copy and paste the entire list of Ensembl IDs into the box; don't upload the file
switch build_gene_table
    case 2
        gene_table = readtable(gene_table_in, 'FileType', 'text', 'Delimiter', '\t');
        gene_map = containers.Map(gene_table.GeneStableID, gene_table.GeneName);
        gene_name_map = containers.Map(gene_table.GeneStableID, gene_table.GeneDescription);
        save(gene_table_out, 'gene_table', 'gene_map', 'gene_name_map');
    case 1
        loaded = load(gene_table_out);
        gene_table = loaded.gene_table;
        gene_map = loaded.gene_map;
        gene_name_map = loaded.gene_name_map;
end

%% Build basic dataset
% Decide to keep genes and calculate mean expression from training set only
train_mask = ismember(patients, train_patients);
switch build_dataset
    case 2
        % Delete genes that all patients have the same val (0 probably) or close
        eps = 1e-6;
        [np,ng] = size(M);
        toss = false(1,ng);
        for ig = 1:ng
            Mi = M(train_mask,ig);
            exp_uniq = unique(Mi);
            dev_var = std(Mi)/mean(Mi);
            if isnan(dev_var) || dev_var < eps || length(exp_uniq) == 1
                toss(ig) = true;
            end
        end
        M(:,toss) = [];
        gene_ids(toss) = [];
        gene_descs(toss) = [];
        fprintf('Threw out %i genes where all patients had the same val\n', sum(toss))
        
        % Look at distribution of variances
        devs = std(M,0,1)./mean(M,1);
        figure
        histogram(log10(devs));
        xlabel('log_{10} Std Dev/Mean')
        ylabel('Count')
        title('Distribution of Expression Variance by Gene')
        
        % Calculate mean expression
        Mmean = mean(M(train_mask,:), 1); % keep this for comparing patients against the pop; for fairness, only do this for the kept patients
        
        % Save "complete" expression dataset
        save(baseline_feature_file, 'M', 'Mmean', 'patients', 'gene_ids', 'gene_descs')
    case 1
        load(baseline_feature_file);
end

%% Generate GSEA format files
% http://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats
% N = 10; % Testing: just take a subset of patients
N = size(M,1); % Production: all patients
ng = size(M,2);

% Expression data
%   Basically the original input file format
%   Make a separate file for each analysis
prec = 16;
expr_files = cell(N,1);
labels_files = cell(N,1);

% Preallocate kept gene id<->name mapping
keep = false(ng,1);
gene_symbols = cell(ng,1);
for ig = 1:ng
    gene_id = gene_ids{ig};
    if gene_map.isKey(gene_id)
        keep(ig) = true;
        gene_symbols{ig} = gene_map(gene_id);
    end
end
% Only keep the genes that we have symbols for
gene_symbols(~keep) = [];
Mmean(~keep) = [];
M(:,~keep) = [];
ng = length(Mmean);

% Copy needed extra files to gsea run dir
gsea_jar = 'gsea2-2.2.4.jar';
gene_set_file = 'h.all.v6.0.symbols.gmt';
gene_chip_file = 'GENE_SYMBOL.chip';
copyfile([base_dir gsea_jar], [gsea_inputs_dir gsea_jar])
copyfile([base_dir gene_set_file], [gsea_inputs_dir gene_set_file])
copyfile([base_dir gene_chip_file], [gsea_inputs_dir gene_chip_file])

% Write input files
parfor iN = 1:N
    patient = patients{iN};
    n_classes = 2; % MEAN and patient
    n_samples = 2;
    
    fprintf('Writing files for patient %s\n', patient)
    
    % Expression data file
    % Note: This is the rate-limiting step
    expr_file = [gsea_inputs_dir sprintf('%s_expr.txt', patient)];
    f = fopen(expr_file, 'w');
    header = {'NAME', 'DESCRIPTION', 'MEAN', patient};
    fprintf(f, [strjoin(header, '\t'), '\n']);
    expr_strs = cell(ng,1);
    for ig = 1:ng
        gene_symbol = gene_symbols{ig};
        gene_name = 'na';
        expr_strs{ig} = [gene_symbol '\t' gene_name '\t' num2str(Mmean(ig), prec) '\t' num2str(M(iN,ig), prec)];
    end
    fprintf(f, strjoin(expr_strs, '\n'));
    fclose(f);
    
    % Phenotype labels file
    %   Each class has 1 sample, user-visible names are the same as internal names
    labels_file = [gsea_inputs_dir sprintf('%s_labels.cls', patient)];
    f = fopen(labels_file, 'w');
    fprintf(f, '%i\t%i\t1\n', n_samples, n_classes);
    names = {'MEAN', patient};
    fprintf(f, ['#\t' strjoin(names, '\t'), '\n']);
    fprintf(f, [strjoin(names, '\t'), '\n']);
    fclose(f);
    
    % Run script
    run_file = [gsea_inputs_dir sprintf('run_%i.sh', iN)];
    f = fopen(run_file, 'w');
    fprintf(f, '#!/usr/bin/env bash\n');
    fprintf(f, '# Run GSEA for patient %s\n', patient);
    fprintf(f, 'java -cp gsea2-2.2.4.jar -Xmx1024m xtools.gsea.Gsea -res %s -cls %s#%s_versus_MEAN -gmx %s -collapse false -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set -rnd_type no_balance -scoring_scheme weighted -metric Signal2Noise -sort real -order descending -chip %s -include_only_symbols true -make_sets true -median false -num 100 -plot_top_x 20 -rnd_seed timestamp -save_rnd_lists false -set_max 500 -set_min 15 -zip_report false -out %s -gui false -rpt_label %s\n', expr_file, labels_file, patient, gene_set_file, gene_chip_file, gsea_outputs_dir, patient);
    fclose(f);
    system(['chmod 744 ' run_file]);
end

