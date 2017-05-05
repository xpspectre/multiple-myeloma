from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
import pandas as pd

study_data = pd.read_csv('./data/processed/patient_study.csv', index_col=0)
#clinical_data = pd.read_csv('./data/processed/imputed_8.csv', index_col=0)
clinical_data = pd.read_csv('./data/processed/baseline_clinical_data_imputed_mice.csv', index_col=0)
treatment_data = pd.read_csv('./data/processed/patient_treat.csv', index_col=0)
patient_endp = pd.read_csv('./data/processed/patient_endp.csv', index_col=0)

# Drop column with empty values
treatment_data = treatment_data.drop('sct_bresp', 1)

# Only look at patients in all datasets, who were part of the study
all_ids = list(set(patient_endp.index).intersection(set(treatment_data.index)).intersection(set(clinical_data.index)).intersection(set((study_data.STUDY_ID == 1).index)))

y = pd.notnull(patient_endp[['D_PT_lstalive', 'D_PT_deathdy']].max(axis=1).loc[all_ids])
all_ids = list(set(all_ids).intersection(y.index))
X = clinical_data.join(treatment_data).loc[all_ids]

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=1)

logreg = LogisticRegression(random_state=1, max_iter=1000, n_jobs=3)
logreg.fit(X_train, y_train)
score_train = logreg.score(X_train, y_train)
score_test = logreg.score(X_test, y_test)

print("Training accuracy: "+str(score_train))
print("    Test accuracy: "+str(score_test))

