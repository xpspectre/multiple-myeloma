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
build_gene_table = 1;
build_dataset = 1;

rnaseq_file = 'data/raw/MMRF_CoMMpass_IA9_E74GTF_Cufflinks_Gene_FPKM.txt';
per_visit_lim_out = 'data/processed/per_visit_lim.mat';
raw_file = 'data/processed/rnaseq_cufflinks_fpkm_raw.mat';
baseline_feature_file = 'data/processed/rnaseq_cufflinks_fpkm_baseline_features.mat';
gsea_expr_file = 'data/processed/rnaseq_cufflinks_fpkm_gsea_expr.txt'; % main expression file for GSEA
gsea_phenotype_labels_file = 'data/processed/rnaseq_cufflinks_fpkm_gsea_phenotype_labels.cls'; % main phenotype labels file for GSEA
gene_ids_file = 'data/processed/rnaseq_cufflinks_fpkm_gsea_gene_ids.txt';
gene_table_in = 'gene_table.txt';
gene_table_out = 'data/processed/gene_table.mat';
missing_genes_file = 'data/processed/genes_to_lookup.txt';

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
% Testing: just take a subset of patients
N = 10;
ng = size(M,2);

% Expression data
%   Basically the original input file format
missing_genes = cell(ng,1); % keep track of missing genes to lookup. These will be skipped in the expr file.
prec = 12;
f = fopen(gsea_expr_file, 'w');
header = [{'NAME', 'DESCRIPTION', 'MEAN'}, patients(1:N)];
fprintf(f, [strjoin(header, '\t'), '\n']);
for ig = 1:ng
    vals = [Mmean(ig); M(1:N,ig)]';
    vals_strs = strread(num2str(vals, prec), '%s');
    gene_id = gene_ids{ig};
    if gene_map.isKey(gene_id)
        gene_symbol = gene_map(gene_id);
        gene_name = gene_name_map(gene_id);
        if isempty(gene_name) % yeah, some of the vals are blank
            gene_name = 'na';
        end
    else
        fprintf('Gene %s symbol not found\n', gene_id)
        missing_genes{ig} = gene_id;
        continue
    end
    entries = [gene_symbol, gene_name, vals_strs']; % blank description
    fprintf(f, [strjoin(entries, '\t'), '\n']);
end
fclose(f);
missing_genes = missing_genes(~cellfun(@isempty, missing_genes));
nmg = length(missing_genes);

% Make list of missing genes
%   There are a couple thousand - probably no big deal
f = fopen(missing_genes_file, 'w');
for img = 1:nmg
    fprintf(f, '%s\n', missing_genes{img});
end
fclose(f);

% Phenotype labels
%   Each class has 1 sample, user-visible names are the same as internal names
n_classes = N + 1;
n_samples = n_classes;
f = fopen(gsea_phenotype_labels_file, 'w');
fprintf(f, '%i\t%i\t1\n', n_samples, n_classes);
names = [{'MEAN'}, patients(1:N)];
fprintf(f, ['#\t' strjoin(names, '\t'), '\n']);
fprintf(f, [strjoin(names, '\t'), '\n']);
fclose(f);
