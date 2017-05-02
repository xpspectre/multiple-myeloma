# Run this after prep_clinical_data.py

import os
import pandas as pd
import numpy as np

data_dir = 'data/processed'

# Load study data - to keep just the CoMMpass patients
study_data = pd.read_csv(os.path.join(data_dir, 'patient_study.csv'))
study_data.set_index('PUBLIC_ID', inplace=True)

# Load demographic data
demo_data = pd.read_csv(os.path.join(data_dir, 'patient_data.csv'))
demo_data.set_index('PUBLIC_ID', inplace=True)

# Load visit data
visit_data = pd.read_csv(os.path.join(data_dir, 'clinical_data.csv'))

# Select just the baseline/screening data
visit_data = visit_data.loc[visit_data['VISIT'] <= 0]

# Combine rows for each patient
#   Averaging the rows takes care of missing data/NaNs properly
# unique_ids = data['PUBLIC_ID'].unique()
# print(unique_ids)
visit_data = visit_data.groupby('PUBLIC_ID').mean()
visit_data.drop('VISIT', axis=1, inplace=True)

# Combine demographic and visit data
data = demo_data
data = data.join(visit_data)

# Only keep CoMMpass patients
data = data[study_data['STUDY_ID'] == 1]

# Drop cols that represent change over last visit
data.drop(['AT_INCREASEOF25F', 'AT_SERUMMCOMPONE', 'AT_URINEMCOMPONE', 'AT_ONLYINPATIENT', 'AT_ONLYINPATIENT2', 'AT_DEVELOPMENTOF'], axis=1, inplace=True)

# Keep only the cols that have >= threashold % entries present
#   Set this to determine how much we have to consider missing data
#   A smallish amount of missing data should be pretty straightforward imputation
#   More missing data is harder
#   Which cols can be used for regression of missing data? (if our method requires that)
keep_thres = 0.5
cols = list(data)
present = []
N, N_cols = data.shape
for col in cols:
    n = pd.notnull(data[col]).sum()
    present.append(float(n)/N)
present = np.array(present)
drop_cols = np.where(present < keep_thres)[0]
data.drop(data.columns[drop_cols], axis=1, inplace=True)
print('Dropped {n}/{N} cols that had less than {x} frac of values'.format(n=drop_cols.size, N=N_cols, x=keep_thres))

# Save combined baseline patient data
data.to_csv(os.path.join(data_dir, 'baseline_clinical_data.csv'))

# print(data)
