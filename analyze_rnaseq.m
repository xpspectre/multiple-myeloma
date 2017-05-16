%% Analyze RNAseq GSEA results
clear; close all; clc
rng('default');

% Flags for running expensive ops
raw_enrichment = 2;

gene_set_file = 'data/processed/h.all.v6.0.symbols.gmt';
gsea_outputs_dir = 'data/processed/gsea_outputs2/';

raw_enrichment_file = 'data/processed/gsea_features_raw_enrichment.mat';
enrichment_file = 'data/processed/gsea_features_enrichment.mat';

%% Get gene sets
f = fopen(gene_set_file);
C = textscan(f, '%s\t%*[^\n]'); % just get 1st col
fclose(f);
gene_sets = C{1};
ng = length(gene_sets);

%% Build map of gene set to col index in feature matrix
gene_set_map = containers.Map('KeyType', 'char', 'ValueType', 'double');
for ig = 1:ng
    gene_set = gene_sets{ig};
    gene_set_map(gene_set) = ig + 1; % offset: 1st col is PUBLIC_ID
end

%% Get GSEA outputs
%   Assume there's only 1 entry for each ID
listing = dir(gsea_outputs_dir);
nl = length(listing);
result_dirs = cell(nl,3); % [dir, SUBJECT_ID, timestamp]
for il = 1:nl
    isdir = listing(il).isdir;
    if ~isdir
        continue
    end
    name = listing(il).name;
    tokens = regexp(name, '^(MMRF_\d{4}).Gsea.(\d+)$', 'tokens');
    if isempty(tokens)
        continue
    end
    full_name = [gsea_outputs_dir name];
    result_dirs(il,:) = {[full_name '/'], tokens{1}{1}, tokens{1}{2}};
end

% Cleanup patient dirs: only take real patients and sort by index
result_dirs(all(cellfun(@isempty, result_dirs), 2), :) = []; % remove empty rows
[patients, sort_inds] = sort(result_dirs(:,2));
result_dirs = result_dirs(sort_inds,:);
nl = size(result_dirs,1);

%% Build main table of significantly enriched gene sets
switch raw_enrichment
    case 2
        % Feature table
        enriched = table(patients, 'VariableNames', {'PUBLIC_ID'});
        for ig = 1:ng
            enriched.(gene_sets{ig}) = zeros(nl,1);
        end
        
        % Raw table
        raw = cell(nl,2);
        
        % Read tables of significantly enriched genes in subject id (vs mean)
        % Col headers, with underscore suffix denoting don't interpret
        headers = {'GeneSet', 'GS_', 'Details_', 'Size', 'ES', 'NES', 'NOMpval', 'FDRqval', 'FWERpval', 'RankAtMax', 'LeadingEdge', 'Whatever_'};
        % By default, GSEA outputs gene sets significant to FDR < 0.25 ("strong") and
        %   nominal p-val < 0.05 ("weak"). Code these as 2 and 1, respectively. Later
        %   analysis can do the cutoff solely by the FDR criteria, as desired.
        qval_cutoff = 0.25;
        pval_cutoff = 0.05;
        for il = 1:nl
            result_dir = result_dirs{il,1};
            subject_id = result_dirs{il,2};
            raw{il,1} = subject_id;
            
            tstamp = result_dirs{il,3};
            result_file = [result_dir 'gsea_report_for_' subject_id '_' tstamp '.xls'];
            T = readtable(result_file, 'FileType', 'text', 'Delimiter', '\t', 'HeaderLines', 1, 'ReadVariableNames', false);
            if isempty(T) % no significant hits
                continue
            end
            T.Properties.VariableNames = headers;
            
            raw{il,2} = T; % this will be empty for people who have no sig enriched gene sets at all
            
            ne = height(T);
            for ie = 1:ne
                gene_set = T{ie,'GeneSet'}{1};
                gene_set_ind = gene_set_map(gene_set);
                pval = T{ie,'NOMpval'};
                qval = T{ie,'FDRqval'};
                if qval <= qval_cutoff
                    enriched{il,gene_set_ind} = 2;
                elseif pval <= pval_cutoff % not sig at q but sig at p
                    enriched{il,gene_set_ind} = 1;
                end
            end
        end
        
        save(raw_enrichment_file, 'enriched', 'raw');
    case 1
        loaded = load(raw_enrichment_file);
        enriched = loaded.enriched;
        raw = loaded.raw;
end

%% Postprocess enrichment hits
% Look at distribution of enrichment stats
nsigq = zeros(ng,1); % strong
nsigp = zeros(ng,1); % weak
for ig = 1:ng
    gene_set = gene_sets{ig};
    x = enriched.(gene_set);
    nsigq(ig) = sum(x == 2);
    nsigp(ig) = sum(x == 1);
end

figure
subplot(1,2,1)
histogram(nsigq)
xlabel('# Patients')
ylabel('# Gene Sets')
title('Sig at FDR q-val <= 0.25')
subplot(1,2,2)
histogram(nsigp) % not sig at q but sig at p, ignore these
xlabel('# Patients')
title('Sig at p-val <= 0.05')
suptitle('Number of Significantly Enriched Patients')

% Throw out weak hits - keep only sig at q
for ig = 1:ng
    gene_set = gene_sets{ig};
    enriched.(gene_set) = double(enriched.(gene_set) > 1);
end

% Throw out gene sets that aren't enriched for anyone
for ig = 1:ng
    gene_set = gene_sets{ig};
    x = enriched.(gene_set);
    if all(x == 0)
        fprintf('Dropping gene set %s because no one was enriched\n', gene_set);
        enriched.(gene_set) = [];
    end
end
gene_sets = enriched.Properties.VariableNames(2:end);
fprintf('Kept %i gene sets\n', length(gene_sets));

save(enrichment_file, 'enriched');

