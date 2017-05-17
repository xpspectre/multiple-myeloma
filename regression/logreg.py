# Baseline n-yr survival
# Feed different baseline data files, each with different sets of baseline features

from sklearn.linear_model import LogisticRegression
import pandas as pd
import numpy as np
from sklearn.metrics import roc_curve, auc
import matplotlib.pyplot as plt

data_file = '../data/processed/baseline_data_all.csv'
data = pd.read_csv(data_file)

last_observed = data['LAST_OBSERVED'].values
censored = data['CENSORED'].values
train = data['TRAIN_SET'].values.astype(np.bool)
meta_cols = ['PUBLIC_ID', 'LAST_OBSERVED', 'CENSORED', 'ISS', 'TRAIN_SET']
data.drop(meta_cols, axis=1, inplace=True)
col_names = list(data)
X = data.as_matrix()

# Toss invalid endpoints - NaNs
keep = np.logical_not(np.isnan(last_observed))
X = X[keep,:]
last_observed = last_observed[keep]
censored = censored[keep]
train = train[keep]

# Make survival endpoint
n = 2
n_day = n * 365
dead = (last_observed < n_day) & (np.logical_not(censored))  # not censored means death recorded
keep = (last_observed >= n_day) | dead  # patients who don't have the full measurement timeframe are only kept if dead

y = dead.astype(np.float64)

X = X[keep,:]
y = y[keep]
train = train[keep]

X_train = X[train,:]
X_test = X[np.logical_not(train),:]
y_train = y[train]
y_test = y[np.logical_not(train)]

reg = LogisticRegression(penalty='l1', C=1, random_state=1, max_iter=10000, n_jobs=3)
reg.fit(X_train, y_train)
score_train = reg.score(X_train, y_train)
score_test = reg.score(X_test, y_test)

print("Training accuracy: "+str(score_train))
print("    Test accuracy: "+str(score_test))

# Output fit params in order of importance
coefs = reg.coef_.squeeze()
coefs_sortinds = np.argsort(np.abs(coefs))
for i in coefs_sortinds:
    print(col_names[i] + '\t' + str(coefs[i]))

yhat_train = reg.decision_function(X_train)
fpr_train, tpr_train, _ = roc_curve(y_train, yhat_train)
auc_train = auc(fpr_train, tpr_train)

yhat_test = reg.decision_function(X_test)
fpr_test, tpr_test, _ = roc_curve(y_test, yhat_test)
auc_test = auc(fpr_test, tpr_test)

plt.figure()
lw = 2
plt.plot(fpr_train, tpr_train, color='darkorange', lw=lw, label='Train (area = %0.2f)' % auc_train)
plt.plot(fpr_test, tpr_test, color='blue', lw=lw, label='Test (area = %0.2f)' % auc_test)
plt.plot([0, 1], [0, 1], color='black', lw=lw, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.0])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.legend(loc="lower right")
plt.title('Baseline 1 yr Survival, Clinical Data Only, ROC')
# plt.title('Baseline 1 yr Survival, Clinical + NS Muts Data, ROC')
# plt.title('Baseline 1 yr Survival, Clinical + NS Muts + RNAseq Expr Data, ROC')
# plt.title('Baseline 1 yr Survival, Clinical + RNAseq Expr Data, ROC')
# plt.title('Baseline 1 yr Survival, Clinical + NS Muts + RNAseq GSEA, ROC')
# plt.title('Baseline 1 yr Survival, Clinical + NS Muts + RNAseq Expr + RNAseq GSEA, ROC')
plt.draw()
plt.show()