# Run this after prep_clinical_data.py

import os
import pandas as pd

data_dir = 'data/processed'

data = pd.read_csv(os.path.join(data_dir, 'clinical_data.csv'))

# Select just the baseline/screening data
data = data.loc[data['VISIT'] <= 0]

# Combine rows for each patient
#   Averaging the rows takes care of missing data/NaNs properly
# unique_ids = data['PUBLIC_ID'].unique()
# print(unique_ids)
data = data.groupby('PUBLIC_ID').mean()
data.drop('VISIT', axis=1, inplace=True)

# Combine them
data.to_csv(os.path.join(data_dir, 'baseline_clinical_data.csv'))

print(data)
