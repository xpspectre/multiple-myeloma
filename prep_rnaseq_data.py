# Load, filter, prep RNAseq data for baseline analyses
# Note: This code is deprecated - use the Matlab version prep_rnaseq.m

from load_patient_data import load_per_visit_data
import numpy as np
import pandas as pd
import csv
import cmsgpack
import scipy.io as sio

# Flags for intermediate processing
load_data = False

# First grab per-visit data to get headings for baseline visit
per_visit_data, per_visit_dict, per_visit_fields = load_per_visit_data()
vdata = per_visit_data[['PUBLIC_ID', 'SPECTRUM_SEQ', 'VISIT']]
vdata = vdata.loc[vdata['VISIT'] <= 0]
vdata['SPECTRUM_SEQ'] = vdata['SPECTRUM_SEQ'].replace({'': np.nan})
vdata = vdata.loc[pd.notnull(vdata['SPECTRUM_SEQ'])]
vdata.drop('VISIT', axis=1, inplace=True)
vdata.set_index('PUBLIC_ID', inplace=True)
keep_studies = vdata['SPECTRUM_SEQ'].tolist()

# Load RNAseq file
# Pick the Cufflinks PFKM file as a normalized, easy-to-use starting point
# Do this in 2 stages: prep to get rows and cols, preallocate np array, then fill in data
rnaseq_file = 'data/raw/MMRF_CoMMpass_IA9_E74GTF_Cufflinks_Gene_FPKM.txt'
rnaseq_out = 'data/processed/RNAseq_Cufflinks_Gene_FPKM'

if load_data:
    n_rows = 0
    with open(rnaseq_file) as c:
        r = csv.reader(c, delimiter='\t')
        header = next(r)
        for row in r:
            n_rows += 1

    meas_id = header[2:]
    n_cols = len(meas_id)
    data = np.zeros((n_rows, n_cols), dtype=np.float64)
    gene_ids = []
    gene_pos = []
    with open(rnaseq_file) as c:
        r = csv.reader(c, delimiter='\t')
        next(r)
        for i, row in enumerate(r):
            data[i, :] = row[2:]
            gene_ids.append(row[0])
            gene_pos.append(row[1])

    # The following numpy command to do this is slow and fails for some reason on a short col?
    # data = np.loadtxt(rnaseq_file, delimiter='\t', skiprows=1, usecols=range(2, n_cols+1))

    # Save data in a more convenient format
    with open(rnaseq_out, 'wb') as f:
        f.write(cmsgpack.packb((data, meas_id, gene_ids, gene_pos)))

else:
    with open(rnaseq_out, 'rb') as f:
        data, meas_id, gene_ids, gene_pos = cmsgpack.unpackb(f.read())

# Find cols of data that match desired baseline studies
# First append suffix from the data fileheaders
keep_col_names = []
keep_col_inds = []
n_found = 0
n_miss = 0
for study in keep_studies:
    col = study + '_BM'
    try:
        keep_col_inds.append(meas_id.index(col))  # ValueError on not found
        keep_col_names.append(col)
        n_found += 1
    except ValueError:
        print('Col {col} not found in dataset headers'.format(col=col))
        n_miss += 1
print('{n} per-visit cols found in dataset'.format(n=n_found))
print('{n} per-visit cols not found in dataset'.format(n=n_miss))
# ~1/3 the patients listed SPECTRUM_SEQ not found in dataset

# The other way around
per_visit_cols = set([x + '_BM' for x in keep_studies])
rnaseq_cols = set(meas_id)
print('{n} cols in dataset not kept because not found or later times, etc'.format(n=len(list(rnaseq_cols - per_visit_cols))))
# print(len(list(per_visit_cols - rnaseq_cols)))

# Keep cols
X = data[:, keep_col_inds]

# Quick processing: remove rows that only have a single value for everyone
n_rows = X.shape[0]
all_same = np.zeros((n_rows,), dtype=bool)
for i in range(n_rows):
    if np.unique(X[i, :]).size == 1:
        all_same[i] = True
X = np.delete(X, all_same, axis=0)  # Note: the FutureWarning behavior is what we want
gene_ids = np.delete(gene_ids, all_same)
gene_pos = np.delete(gene_pos, all_same)
print(len(X))
print(len(gene_ids))
print(len(gene_pos))

# Transpose X to get patient x feature matrix
X = X.transpose()

rnaseq_baseline_out = 'data/processed/RNAseq_Cufflinks_Gene_FPKM_baseline'
with open(rnaseq_baseline_out, 'wb') as f:
    f.write(cmsgpack.packb((X, keep_col_names, gene_ids, gene_pos)))

rnaseq_baseline_out_mat = 'data/processed/RNAseq_Cufflinks_Gene_FPKM_baseline.mat'
sio.savemat(rnaseq_baseline_out_mat, {'X': X,
                                      'meas_ids': keep_col_names,
                                      'gene_ids': gene_ids,
                                      'gene_pos': gene_pos})