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

load_raw = 0;
build_gene_table = 1;
build_dataset = 1;

rnaseq_file = 'data/raw/MMRF_CoMMpass_IA9_E74GTF_Cufflinks_Gene_FPKM.txt';
per_visit_lim_out = 'data/processed/per_visit_lim.mat';
raw_file = 'data/processed/rnaseq_cufflinks_fpkm_raw.mat';
baseline_feature_file = 'data/processed/rnaseq_cufflinks_fpkm_baseline_features.mat';
gene_ids_file = 'data/processed/rnaseq_cufflinks_fpkm_gsea_gene_ids.txt';
gene_table_in = 'gene_table.txt';
gene_table_out = 'data/processed/gene_table.mat';
missing_genes_file = 'data/processed/genes_to_lookup.txt';

gsea_inputs_dir = 'data/processed/gsea/';
gsea_outputs_dir = 'data/processed/gsea_outputs/';
gsea_cmd_file = 'rnaseq_cufflinks_fpkm_gsea_cmds.sh'; % Generate script to run GSEA in an automated fashion

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
switch build_dataset
    case 2
        
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
        devs = std(M,0,1)./mean(M,1);
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
        Mmean = mean(M,1); % keep this for comparing patients against the pop; for fairness, only do this for the kept patients
        patients = patients(ib);
        patients = regexprep(patients,'.....$',''); % delete the last 5 chars to get the PUBLIC_ID to join on (probably)
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
prec = 12;
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

% Write actual files
%   This is inefficient, but good enough
for iN = 1:N
    patient = patients{iN};
    n_classes = 2; % MEAN and patient
    n_samples = 2;
    
    fprintf('Writing files for patient %s\n', patient)
    
    % Expression data file
    expr_file = [gsea_inputs_dir sprintf('rnaseq_cufflinks_fpkm_gsea_expr_%s.txt', patient)];
    f = fopen(expr_file, 'w');
    header = {'NAME', 'DESCRIPTION', 'MEAN', patient};
    fprintf(f, [strjoin(header, '\t'), '\n']);
    for ig = 1:ng
        vals = [Mmean(ig); M(iN,ig)]';
        vals_strs = strread(num2str(vals, prec), '%s');
        gene_symbol = gene_symbols{ig};
        gene_name = 'na';
        entries = [gene_symbol, gene_name, vals_strs']; % blank description
        fprintf(f, [strjoin(entries, '\t'), '\n']);
    end
    fclose(f);
    expr_files{iN} = expr_file;
    
    % Phenotype labels file
    %   Each class has 1 sample, user-visible names are the same as internal names
    labels_file = [gsea_inputs_dir sprintf('rnaseq_cufflinks_fpkm_gsea_phenotype_labels_%s.cls', patient)];
    f = fopen(labels_file, 'w');
    fprintf(f, '%i\t%i\t1\n', n_samples, n_classes);
    names = {'MEAN', patient};
    fprintf(f, ['#\t' strjoin(names, '\t'), '\n']);
    fprintf(f, [strjoin(names, '\t'), '\n']);
    fclose(f);
    labels_files{iN} = labels_file;
    
%     if iN == 3
%         break
%     end
end

%% Command to run GSEA
% Make sure the gsea jar is in this dir and run this script from this dir.
gene_set_file = 'data/processed/h.all.v6.0.symbols.gmt';
gene_chip_file = 'data/processed/GENE_SYMBOL.chip';
f = fopen(gsea_cmd_file, 'w');
fprintf(f, '#!/usr/bin/env bash\n');
fprintf(f, '# Script to run GSEA on all the patients vs MEAN\n');
for iN = 1:N
    patient = patients{iN};
    expr_file = expr_files{iN};
    label_file = labels_files{iN};
    fprintf(f, 'java -cp gsea2-2.2.4.jar -Xmx1024m xtools.gsea.Gsea -res %s -cls %s#%s_versus_MEAN -gmx %s -collapse false -mode Max_probe -norm meandiv -nperm 1000 -permute gene_set -rnd_type no_balance -scoring_scheme weighted -metric Signal2Noise -sort real -order descending -chip %s -include_only_symbols true -make_sets true -median false -num 100 -plot_top_x 20 -rnd_seed timestamp -save_rnd_lists false -set_max 500 -set_min 15 -zip_report false -out %s -gui false -rpt_label %s\n', expr_file, label_file, patient, gene_set_file, gene_chip_file, gsea_outputs_dir, patient);
    
%     if iN == 3
%         break
%     end
end
fclose(f);


