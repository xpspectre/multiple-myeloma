# Run this after prep_clinical_data.py

import os
import pandas as pd

data_dir = 'data/processed'

data = pd.read_csv(os.path.join(data_dir, 'clinical_data.csv'))

# Select just the baseline/screening data
data = data.loc[data['VISIT'] <= 0]

# Combine them
data.to_csv(os.path.join(data_dir, 'baseline_clinical_data.csv'))

print(data)
