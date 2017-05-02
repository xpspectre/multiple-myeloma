# Multiple Myeloma Project

Note: I used Python 3 for this.

## Setup

Make a new Anaconda3 environment
```
conda create -n mm-env python=3.5 numpy scipy matplotlib pandas msgpack-python scikit-learn pillow xlrd lifelines
source activate mm-env
pip install tensorflow msgpack-numpy
```

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

## Survival Analysis

Some basic survival analysis from baseline data is in `analyze_endpoints.py`.