%% Find genes from RNAseq expression vals that are individually correlated to outcomes
% Outcomes: overall survival, survival at x yrs
% Obviously run this after prep_rnaseq.m
% Only use results from patients in the training set
clear; close all; clc
rng('default');

run_survival = 2;

test_train_split_file = 'data/processed/test_train_split.mat';
baseline_feature_file = 'data/processed/rnaseq_cufflinks_fpkm_baseline_features.mat';
endpoints_file = 'data/processed/baseline_clinical_endp.csv';
gene_table_file = 'data/processed/gene_table.mat';
sig_gene_output_file = 'data/processed/rnaseq_cufflinks_fpkm_baseline_sig_gene_expr.mat';

survival_results_file = 'data/processed/rnaseq_analysis_survival.mat';

%% Load train-test split
loaded = load(test_train_split_file);
train_patients = loaded.patients(loaded.train_inds);

%% Load expression data
% The patients are sorted here
loaded = load(baseline_feature_file);
M = loaded.M;
gene_ids = loaded.gene_ids;
patients = loaded.patients(:); % these are the proper patient IDs
[N, ng] = size(M);

%% Load survival data
endp = readtable(endpoints_file);
endp_patients = endp.PUBLIC_ID;
[endp_patients, sort_ind] = sort(endp_patients);
endp = endp(sort_ind,:);

%% Match up endpoints to RNAseq data
combined_patients = intersect(patients, endp_patients);
keep_expr = ismember(patients, combined_patients);
keep_endp = ismember(endp_patients, combined_patients);
M = M(keep_expr,:);
endp = endp(keep_endp,:);
patients = patients(keep_expr);
train_mask = ismember(patients, train_patients);

died = endp{:,'D_PT_deathdy'};
last_observed = max(died, endp{:,'D_PT_lstalive'});
censored = isnan(died); % 1 = censored, 0 = death observed

% Preprocessing: z-score on training set
[Mz, mus, sigmas] = zscore(M(train_mask,:));

% Same training set observations for each var
y = last_observed(train_mask);
censored = censored(train_mask);

% Run individual gene survivals
% ng = 3; % test subset of genes
switch run_survival
    case 2
        bs = zeros(ng,1);
        ps = zeros(ng,1);
        errs = false(ng,1); % report if an error occurred
        parfor ig = 1:ng
            X = Mz(:,ig);
            if length(unique(X)) < 10 % handles case where everyone is the same; everyone's the same except a couple of outliers, etc
                errs(ig) =  true
                continue
            end
            [b, logL, H, stats] = coxphfit(X, y, 'Censoring', censored);
            bs(ig) = b;
            ps(ig) = stats.p;
        end
        
        save(survival_results_file, 'bs', 'ps', 'errs');
    case 1
        load(survival_results_file);
end

all_gene_ids = gene_ids; % keep for reference

% First toss errored results
ps = ps(~errs);
bs = bs(~errs);
gene_ids = gene_ids(~errs);

% Use Bonferroni corrected p-vals
p_cutoff = 0.05;
% p_cutoff_bon = p_cutoff/ng;
% Use false discovery rate q-vals
[fdr, qs] = mafdr(ps);

% Only keep sig results
% keep = ps < p_cutoff_bon;
keep = qs < p_cutoff;
ps = ps(keep);
qs = qs(keep);
bs = bs(keep);
gene_ids = gene_ids(keep);

% Output sig results in order of b magnitude
[~, sort_ind] = sort(abs(bs), 'descend');
ps = ps(sort_ind);
qs = qs(sort_ind);
bs = bs(sort_ind);
gene_ids = gene_ids(sort_ind);

% Lookup gene IDs
nb = length(bs);
loaded = load(gene_table_file);
gene_map = loaded.gene_map;
gene_name_map = loaded.gene_name_map;
gene_symbols = cell(nb,1);
gene_descs = cell(nb,1);
for ib = 1:nb
    gene_id = gene_ids{ib};
    if gene_map.isKey(gene_id)
        gene_symbols{ib} = gene_map(gene_id);
        gene_descs{ib} = gene_name_map(gene_id);
        if isempty(gene_descs{ib})
            gene_descs{ib} = '??';
        end
    else
        gene_symbols{ib} = '?';
        gene_descs{ib} = '??';
    end
end

% Print results
fprintf('beta\t\tq-val\t\tGene ID\t\tGene Symbol\tDesc\n')
for ib = 1:nb
    fprintf('%8.3f\t%e\t%s\t%s\t%s\n', bs(ib), qs(ib), gene_ids{ib}, gene_symbols{ib}, gene_descs{ib});
end
% Result: The only sig genes all increase hazard

% Create table of these values ready to join to other features
T = table(patients, 'VariableNames', {'PUBLIC_ID'});
for ib = 1:nb
    gene_id = gene_ids{ib};
    col_name = [gene_id '_expr'];
    col_pos = find(ismember(all_gene_ids, gene_id));
    T.(col_name) = M(:,col_pos);
end

save(sig_gene_output_file, 'T', 'bs', 'ps', 'qs', 'gene_ids', 'gene_symbols', 'gene_descs')

