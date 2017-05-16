%% Analyze per-visit data over time
% No genomic data for now (or ever: we don't have much of that during follow-up visits)
% Setup: fill in missing time data, break into 3 month chunks for analysis
% Save current TREATMENTRESP as PROGRESSION_CURR - this is a feature that
%   indicates current severity
% Assesses death as a possible progression
%   TODO: Separate MM-caused and non-MM (another form of right censoring) death
%       But this can be hard to judge
clear; close all; clc
rng('default');

% Switches for expensive steps
interval_avg = 1;
impute = 1;
add_deltas = 1;
assemble_pred_prob = 2;

interval_avg_file = 'data/processed/timeseries_interval_avg.mat';
impute_file = 'data/processed/timeseries_impute.mat';
add_deltas_file = 'data/processed/timeseries_add_deltas.mat';
assembled_pred_prob_file = 'data/processed/timeseries_assembled_pred_prob.mat';

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
treat.patients = cellstr(treat.patients);
treat.therapies = cellstr(treat.therapies);
treat.intervals = double(treat.intervals);
treat.data = double(treat.data);

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

% Prep survival data
endp = sortrows(endp, 'PUBLIC_ID');
died = endp{:,'D_PT_deathdy'};
last_observed = max(died, endp{:,'D_PT_lstalive'});
censored = isnan(died); % 1 = censored, 0 = death observed
survival = table(endp.PUBLIC_ID, last_observed, censored, 'VariableNames', {'PUBLIC_ID','LAST_OBSERVED','CENSORED'});

%% Combine and process data
% - Each patient has a set of unchanging features ("demographic data")
% - Each patient has changing features ("lab data") at irregular intervals
%   Augment t = 0 with imputed baseline data
% - Each patient has treatments ongoing during specified intervals
%   Treatments lag 1 interval behind the measurements - treatment given in the 90 days right
%       before the prediction occurs are counted. Treatment given in the 90
%       days before the

% Joint treatment response to per-visit data - just append them - they came
%   from the same place
time_data.TREATMENTRESP = time_endp.AT_TREATMENTRESP;

% Markov assumption:
%   - Baseline data + 1st treatment interval (day 1-90) used to predict
%   survival/progression at 90 days
%   - 1st 90 days data + 2nd treatment interval (day 91-180) used to
%   predict survival/progression at 180 days
%   - 2nd 90 days data (day 91-180) + 3rd treatment interval (day 181-270)
%   used to predict survival/progression at 270 days
%   - ...

% Convert VISITDY day to intervals
%   For now, use VISIT as the time
%       More rigorously, we could go back and get the precise collection of
%       the particular measurement by type, but the coarse-grained
%       intervals here don't need that.
%   VISITDY <= 0 assigned to interval 0
%   Interval = 90 is the 1st measurement after baseline, encompassing the
%       last 90 days
% Also convert last observed day to last observed interval
intervals = treat.intervals;
% Drop VISIT col
time_data.VISIT = [];
% Drop all rows whose VISITDY is NaN - these are the VISIT = 20, filled with NaNs
time_data(isnan(time_data.VISITDY),:) = [];
% Assign intervals - this is dumb, but works. At least vectorizing v and
%   looping thru the intervals is fast
interval_lo = intervals(1:end-1);
interval_hi = intervals(2:end);
ni = length(interval_lo);

v = time_data.VISITDY;
nv = length(v);
x = zeros(nv,1); % baseline (v <= 0 implicitly given interval 0)

last_obs = survival.LAST_OBSERVED;
nobs = length(last_obs);
y = zeros(nobs,1);

for i = 1:ni
    lo = interval_lo(i);
    hi = interval_hi(i);
    
    in = lo < v & v <= hi;
    x(in) = i;
    
    in = lo < last_obs & last_obs <= hi;
    y(in) = i;
end
time_data.VISITDY = [];
time_data.INTERVAL = x;

survival.LAST_OBS_INTERVAL = y;

%% Average all values for the same patient, same interval, ignoring NaNs
%   (except when all values are NaN/missing)
% It's acceptable for TREATMENTRESP to be averaged as well here for now
switch interval_avg
    case 2
        [g, t_data] = findgroups(time_data(:,{'PUBLIC_ID', 'INTERVAL'}));
        x = time_data;
        x.PUBLIC_ID = [];
        x.INTERVAL = [];
        % Splitapply each column individually
        time_data_cols = x.Properties.VariableNames;
        n_cols = length(time_data_cols);
        for i = 1:n_cols
            col = time_data_cols{i};
            t_data.(col) = splitapply(@nanmean, x.(col), g);
        end
        save(interval_avg_file, 't_data');
    case 1
        loaded = load(interval_avg_file);
        t_data = loaded.t_data;
end

%% Only keep cols that have > x% data present
present_cutoff = 0.5;
always_keep = {'PUBLIC_ID','INTERVAL','TREATMENTRESP'};
cols = t_data.Properties.VariableNames;
cols = setdiff(cols, always_keep);
n_cols = length(cols);
n = height(t_data);
for i = 1:n_cols
    col = cols{i};
    x = t_data.(col);
    present = ~isnan(x);
    if sum(present)/n < present_cutoff
        fprintf('Tossing col %s due to not enough data\n', col)
        t_data.(col) = [];
    end
end
n_cols = width(t_data) - 3;
fprintf('%i cols kept with > %g frac data present\n', n_cols, present_cutoff)

%% Impute missing data
% For each patient, fill in missing data, including TREATMENTRESP
% Record imputed vals so plots later can note them
% TODO: Figure out blocks that are all NaNs - patient didn't have ANY (or
%   only had 1) measurements at any visit
% TODO: Alternative is Gaussian process or something
switch impute
    case 2
        always_keep = {'PUBLIC_ID','INTERVAL'};
        cols = t_data.Properties.VariableNames;
        cols = setdiff(cols, always_keep);
        n_cols = length(cols);
        
        % Keep a companion table that tracks present/absent data
        t_present = t_data;
        for i = 1:n_cols
            col = cols{i};
            t_present.(col) = ~isnan(t_present.(col));
        end
        
        patients = unique(data.PUBLIC_ID);
        np = length(patients);
        for ip = 1:np
            patient = patients{ip};
            patient_pos = strcmp(t_data.PUBLIC_ID, patient);
            t_data_i = t_data(patient_pos,:);
            x_ = t_data_i.INTERVAL; % sample points
            for i = 1:n_cols
                col = cols{i};
                t_data_ii = t_data_i.(col);
                missing = isnan(t_data_ii);
                x = x_(~missing);
                v = t_data_ii(~missing);
                xq = x_(missing);
                try
                    vq = interp1(x, v, xq, 'linear', 'extrap');
                catch
                    % do nothing - leave as all NaN (or all but 1 NaN)
                    vq = nan(sum(missing),1);
                end
                t_data_ii(missing) = vq;
                t_data_i.(col) = t_data_ii;
            end
            t_data(patient_pos,:) = t_data_i;
        end
        save(impute_file, 't_data', 't_present');
    case 1
        loaded = load(impute_file);
        t_data = loaded.t_data;
        t_present = loaded.t_present;
end

%% Postprocess imputed data
% The imputation method is dumb and doesn't respect limits of col values
cols = t_data.Properties.VariableNames(3:end-1); % except PUBLIC_ID, INTERVAL, and TREATMENTRESP
n_cols = length(cols);
% All selected cols need to be non-negative
for i = 1:n_cols
    col = cols{i};
    t_data{:,col} = max(t_data{:,col}, 0);
end

% ECOG is only from 0 to 5
t_data{:,'ECOG_PERFORMANCEST'} = min(t_data{:,'ECOG_PERFORMANCEST'}, 5);

% Treatment response is from 1 to 6
t_data{:,'TREATMENTRESP'} = max(t_data{:,'TREATMENTRESP'}, 1);
t_data{:,'TREATMENTRESP'} = min(t_data{:,'TREATMENTRESP'}, 6);

% Threshold binary variables
bin_cols = {'PN_PERIPHERALSEN', 'PN_PERIPHERALMOT', 'RU_DIDTHEPATIENT2', 'RU_DIDTHEPATIENT', ...
    'SS_DOESTHEPATIEN', 'SS_SPINALCORDCOM', 'SS_BONEPAIN', 'SS_FATIGUE', ...
    'SS_HYPERCALCEMIA', 'SS_RENALINSUFFIC', 'SS_ANEMIAHEMOGLO', 'SS_BONELESIONSLY', ...
    'SS_BONELESIONSOS', 'SS_SOFTTISSUEPLA', 'SS_OFTTISSUEPLAS', 'SS_RECURRENTBACT', ...
    'SS_SYMPTOMATICHY', 'SS_AMYLOIDOSIS', 'SS_OTHER' };
n_bin_cols = length(bin_cols);
for i = 1:n_bin_cols
    col = bin_cols{i};
    t_data{:,col} = double(t_data{:,col} > 0.5);
end

%% Augment with change in parameter since last measurement
% To try to better satisfy the Markov property
switch add_deltas
    case 2
        cols = t_data.Properties.VariableNames(3:end); % all cols except PUBLIC_ID and INTERVAL
        n_cols = length(cols);
        delta_cols = strcat('delta_', cols);
        nr = height(t_data);
        delta_vals = zeros(nr,n_cols);
        for ir = 1:nr
            patient = t_data{ir,'PUBLIC_ID'}{1};
            interval = t_data{ir,'INTERVAL'};
            
            curr_vals = t_data{ir,3:end};
            
            % Go back until a prev interval is found or notice you're a baseline measurement
            prev_interval = interval - 1;
            prev_vals = zeros(1,n_cols);
            while prev_interval >= 0 % going negative indicates that was the 1st interval (probably baseline)
                row = find(strcmp(t_data.PUBLIC_ID, patient) & t_data.INTERVAL == prev_interval);
                if ~isempty(row) % found previous interval
                    prev_vals = t_data{row,3:end};
                    break
                end
                prev_interval = prev_interval - 1; % previous interval missing, go back another
            end
            delta_vals(ir,:) = curr_vals - prev_vals;
        end
        
        delta_table = array2table(delta_vals, 'VariableNames', delta_cols);
        t_data = [t_data, delta_table];
        
        save(add_deltas_file, 't_data');
    case 1
        loaded = load(add_deltas_file);
        t_data = loaded.t_data;
    case 0
        % don't add deltas
end

%% Build prediction problem
switch assemble_pred_prob
    case 2
        patients = unique(data.PUBLIC_ID);
        np = length(patients);
        
        % Placeholders for therapies
        nr = height(t_data);
        therapies = treat.therapies;
        n_th = length(therapies);
        for i = 1:n_th;
            therapy = therapies{i};
            t_data.(therapy) = nan(nr,1);
        end
        
        % Placeholder for progressions
        t_data.PROGRESSION = nan(nr,1);
        t_data.PROGRESSION_CURR = nan(nr,1); % this is just TREATMENTRESP
        t_data.DIED = zeros(nr,1); % whether the patient died (1) in this interval or not (0)
        
        % Modify t_data so that each row is an observation, where
        %   - the last col PROGRESSION is the output
        %   - the other cols are features
        for ir = 1:height(t_data);
            % Get measurements at interval (e.g., 0, 1)
            patient = t_data{ir,'PUBLIC_ID'}{1};
            interval = t_data{ir,'INTERVAL'};
            
            % Get treatment at next interval (e.g., 0, 1) (treatment data offset by 1 interval)
            treat_ind = find(ismember(treat.patients, patient));
            treat_data = treat.data(treat_ind,:,interval+1); % 1-based indexing
            for i = 1:n_th;
                therapy = therapies{i};
                t_data{ir,therapy} = treat_data(i);
            end
            
            % Get progression at current interval (e.g., 0, 1)
            t_data{ir,'PROGRESSION_CURR'} = t_data{ir,'TREATMENTRESP'};
            
            % Get death status at this interval
            survival_ind = find(strcmp(survival.PUBLIC_ID, patient));
            patient_last_interval = survival{survival_ind,'LAST_OBS_INTERVAL'};
            patient_censored = survival{survival_ind,'CENSORED'};
            death_status = false;
            if interval == patient_last_interval && patient_censored == 0
                t_data{ir,'DIED'} = 1;
                death_status = true;
            end
            
            % Get progression at next interval (e.g., 1, 2)
            next_interval = interval + 1;
            row = find(strcmp(t_data.PUBLIC_ID, patient) & t_data.INTERVAL == next_interval);
            if isempty(row) % Next result is unknown
                if death_status % Unknown because died
                    prog = -Inf; % should be unique to denote very bad
                else % Unknown because missing
                    continue
                end
            else % Next result is known
                prog = t_data{row,'TREATMENTRESP'};
            end
            t_data{ir,'PROGRESSION'} = prog;
        end
        t_data.TREATMENTRESP = []; % Remove missing cols
        
        % Drop rows where there's no output - these can't be predicted
        drop_rows = isnan(t_data.PROGRESSION);
        t_data(drop_rows,:) = [];
        t_present(drop_rows,:) = [];
        save(assembled_pred_prob_file, 't_data', 't_present');
    case 1
        loaded = load(assembled_pred_prob_file);
        t_data = loaded.t_data;
        t_present = loaded.t_present;
end
