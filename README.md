# Multiple Myeloma Project

Note: I used Python 3 for this.

## Setup

Make a new Anaconda3 environment
```
conda create -n mm-env python=3.5 numpy scipy matplotlib pandas msgpack-python scikit-learn pillow xlrd lifelines
source activate mm-env
pip install tensorflow msgpack-numpy fancyimpute
```

Warning: `fancyimpute` has a lot of dependencies. So it's use has been suppressed for now. Use R's MICE package or similar to do imputation for now.

Make a symlink in this directory to the data store
```
ln -s data <mm data dir>
```

Inside the data store, the subdir `raw/clinical_data_tables/CoMMpass_IA9_FlatFile_Dictionaries` and `raw/clinical_data_tables/CoMMpass_IA9_FlatFiles` contains the unzipped data.

## Baseline Data Prep

The following scripts place new csv files in the data store subdir `processed/`.

1. Run `prep_patient_data.py` to clean up the per-patient cols and generate a main patient demographic data file.
2. Run `prep_clinical_data.py` to clean up the per-visit cols and generate a main clinical data file. See the comments in the file for details. A bunch of cols are dropped from the main data table and split off into their own tables, which can be further processed as needed.
3. Run `prep_baseline_clinical_data.py` to get the aggregated baseline/screening data.
4. Decide which endpoint you're interested in and grab it from the result of `prep_patient_data.py` or `prep_clinical_data.py`.
5. Run `split_train_test.m` to generate single list of patients and train-test indices for all analyses.

## Baseline Genomics Data Prep

1. Run `prep_mut.m` to get incidence of common mutations. Generates mutations-per-patient for GO enrichment analysis (which may be useful).
2. Run `prep_rnaseq.m` to process RNAseq expression data and generate input files and script for GSEA. To run GSEA, download the gene lists and gene chip from the Broad site.
3. Run `analyze_rnaseq.m` to collect the RNAseq results run on the cluster into a table of patient x gene_set enriched.
4. Run `analyze_rnaseq_individual.m` to get a subset of genes whose expression individually correlate to survival and turn those into features directly. 

## Timeseries Data Prep

This requires all the baseline data.

1. Run `analyze_treatments.py` to collect treatments over time and turn into a tensor of patient x treatment x time interval.
2. Run `analyze_timeseries.m` to get a patients x features + outcomes matrix suitable for simple prediction tasks (like disease progression at the next interval given the current interval and treatment regimen).

## Survival Analysis

Some basic survival analysis from baseline data is in `analyze_endpoints.py`.