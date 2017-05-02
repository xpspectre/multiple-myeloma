# Run this after prep_clinical_data.py

import os
import pandas as pd

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

# Save combined baseline patient data
data.to_csv(os.path.join(data_dir, 'baseline_clinical_data.csv'))

print(data)
