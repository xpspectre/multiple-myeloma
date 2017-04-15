# Load and preprocess cytogenetic data pulled from IA9 explorer
# http://nikgrozev.com/2015/07/01/reshaping-in-pandas-pivot-pivot-table-stack-and-unstack-explained-with-pictures/

import os
import pandas as pd
import numpy as np


def pivot_dataset(df, col, exclude=[], missing=[]):
    """Pivot dataframe df of Patient with column col of categoricals into n_categories columns of 1/0. Specify
    list of cols to exclude (such as redundant ones). Specify list of missing col names that indicate the patient
    data is missing or noninformative."""
    x = df.pivot(index='Patient', columns=col, values=col)

    # Drop redundant columns
    x.drop(exclude, axis=1, inplace=True)

    # Replace strs with numbers
    x.replace({None: 0}, inplace=True)
    x[x != 0] = 1
    x = x.apply(pd.to_numeric, errors='ignore', downcast='float')

    # Set patients with Missing = 1 to all NaNs and drop the Missing row
    for mis in missing:
        x.loc[x[mis] == 1, :] = np.nan
        x.drop(mis, axis=1, inplace=True)

    return x


def merge_datasets_merge_cols(df1, df2):
    """Merge (outer join) dataframe df2 onto df1, combining columns present in both dataframs using logical OR
    (if at least one of the dataframes has it being 1, set to 1 in the merged dataframe)"""
    x = df1.join(df2, how='outer', rsuffix='_meta')

    abn_fish_cols = list(df1)
    abn_meta_cols = list(df2)
    for col in abn_meta_cols:
        if col in abn_fish_cols:
            meta_col = col + '_meta'
            both_nans = pd.isnull(x[col]) & pd.isnull(x[meta_col])
            x[col] = (pd.notnull(x[col]) & x[col] == 1) | (pd.notnull(x[meta_col]) & x[meta_col] == 1)  # not NaN and one of the cols is true
            x.loc[both_nans, col] = np.nan  # restore NaNs if both dataframes were unknown
            x.drop(meta_col, axis=1, inplace=True)

    for col in list(x):
        x[col] = x[col].astype(float)

    return x


if __name__ == '__main__':
    # Files
    data_dir = 'data/raw/cytogenetic_data'
    results_dir = 'data/processed'

    abn_fish_base = 'mmrf_explore_ia9_cytogenetics_abnormalities_fish_'
    abn_fish_times = [
        'baseline',
        'month3',
        'month6',
        'month9',
        'month12'
    ]

    abn_meta_base = 'mmrf_explore_ia9_cytogenetics_abnormalities_meta_'
    abn_meta_times = [
        'baseline'
    ]

    tris_base = 'mmrf_explore_ia9_cytogenetics_trisomies_'
    tris_times = [
        'baseline'
    ]

    aneu_base = 'mmrf_explore_ia9_cytogenetics_aneuploidy_'
    aneu_times = [
        'baseline'
    ]

    # Load baseline FISH data
    abn_fish_raw = pd.read_csv(os.path.join(data_dir, abn_fish_base + abn_fish_times[0] + '.csv'))
    abn_fish = pivot_dataset(abn_fish_raw, 'Abnormalities (FISH) - Baseline', exclude=['Multiple', 'None'], missing=['Missing'])

    # Load baseline META data
    abn_meta_raw = pd.read_csv(os.path.join(data_dir, abn_meta_base + abn_meta_times[0] + '.csv'))
    abn_meta = pivot_dataset(abn_meta_raw, 'Abnormalities (Meta) - Baseline', exclude=['Other'])

    # May need special handling of "del17 or p53" and "del 17" cols

    # Combine FISH and META. An abnormality is present if it's a 1 in either of the datasets
    abn_comb = merge_datasets_merge_cols(abn_fish, abn_meta)

    # Grab the Other abnormalities in the META dataset
    abn_meta_other = abn_meta_raw.loc[pd.notnull(abn_meta_raw['Other']), ['Patient', 'Other']]
    abn_meta_other.set_index('Patient', inplace=True)

    # Load trisomies
    tris_raw = pd.read_csv(os.path.join(data_dir, tris_base + tris_times[0] + '.csv'))
    tris_raw.drop_duplicates(subset=['Patient', 'Trisomies - Baseline'], take_last=True, inplace=True)  # Delete duplicate info like in row 114 (of the csv)
    tris = pivot_dataset(tris_raw, 'Trisomies - Baseline', exclude=[], missing=['Not Done'])

    # Rename trisomies columns
    tris.rename(columns=dict([(i, 'Trisomy_' + i) for i in list(tris)]), inplace=True)

    # Combine FISH+META and trisomies
    abn_comb = abn_comb.join(tris, how='outer')

    # Load aneuploidies
    aneu_raw = pd.read_csv(os.path.join(data_dir, aneu_base + aneu_times[0] + '.csv'))
    aneu = pivot_dataset(aneu_raw, 'Aneuploidy - Baseline', exclude=['None'], missing=['Not Done/Unknown'])

    # Combine
    abn_comb = abn_comb.join(aneu, how='outer')

    # Add back in the meta Other text col
    abn_comb = abn_comb.join(abn_meta_other, how='outer')

    # Save results
    print(abn_comb)
    abn_comb.to_csv(os.path.join(results_dir, 'cytogenetic_data_details_processed' + '.csv'))

    # Process summary file separately - coarser features
    summ_base = 'mmrf_explore_ia9_cytogenetics_cytogenetic_summary_'
    summ_times = [
        'baseline'
    ]
    summ_raw = pd.read_csv(os.path.join(data_dir, summ_base + summ_times[0] + '.csv'))
    summ = pivot_dataset(summ_raw, 'Cytogenetic Summary - Baseline', exclude=['Normal'], missing=['Not Done/Unknown', 'Indeterminant'])
    print(summ)
    summ.to_csv(os.path.join(results_dir, 'cytogenetic_data_summary_processed' + '.csv'))
