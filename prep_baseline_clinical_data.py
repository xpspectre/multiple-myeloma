# Run this after prep_clinical_data.py

import os
import pandas as pd
import numpy as np
# from fancyimpute import KNN, MICE

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

# Load endpoints and join/split with data
endp_data = pd.read_csv(os.path.join(data_dir, 'patient_endp.csv'))
endp_data.set_index('PUBLIC_ID', inplace=True)
endp_cols = list(endp_data)
data_ = data
data_ = data_.join(endp_data)
data = data_.drop(endp_cols, axis=1)
endp_data = data_[endp_cols]

# Save combined baseline patient data
data.to_csv(os.path.join(data_dir, 'baseline_clinical_data.csv'))
endp_data.to_csv(os.path.join(data_dir, 'baseline_clinical_endp.csv'))

# Impute missing data
#   If all the cols are allowed to be treated as numeric vals, then this is OK as is
#   Otherwise, if some cols still need to be categorical/indicator, then threshold and convert
# Not sure if these funs below are supposed to return multiple datasets?
# May want to recombine categorical cols into 1 col, then multinomial or softmax logistic regression on them in MI,
#   then resplit
do_mi = False
if do_mi:
    cols = list(data)
    inds = data.index.values

    X = data.as_matrix()

    X_filled_knn = KNN(k=3).complete(X)
    data_filled_knn = pd.DataFrame(data=X_filled_knn, columns=cols)
    data_filled_knn.insert(0, 'PUBLIC_ID', inds)
    data_filled_knn.set_index('PUBLIC_ID', inplace=True)

    X_filled_mice = MICE().complete(X)
    data_filled_mice = pd.DataFrame(data=X_filled_mice, columns=cols)
    data_filled_mice.insert(0, 'PUBLIC_ID', inds)
    data_filled_mice.set_index('PUBLIC_ID', inplace=True)

    # Save imputed data ready for standard analysis
    data_filled_knn.to_csv(os.path.join(data_dir, 'baseline_clinical_data_imputed_knn.csv'))
    data_filled_mice.to_csv(os.path.join(data_dir, 'baseline_clinical_data_imputed_mice.csv'))
