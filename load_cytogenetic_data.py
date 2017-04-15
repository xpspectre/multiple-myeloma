# Load and preprocess cytogenetic data pulled from IA9 explorer
# http://nikgrozev.com/2015/07/01/reshaping-in-pandas-pivot-pivot-table-stack-and-unstack-explained-with-pictures/

import os
import pandas as pd
import numpy as np


def pivot_dataset(df, col, exclude=[], missing=None):
    """Pivot dataframe df of Patient with column col of categoricals into n_categories columns of 1/0. Specify
    list of cols to exclude (such as redundant ones). missing is the value that indicates no data for that
    patient."""
    x = df.pivot(index='Patient', columns=col, values=col)

    # Drop redundant columns
    x.drop(exclude, axis=1, inplace=True)

    # Replace strs with numbers
    x.replace({None: 0}, inplace=True)
    x[x != 0] = 1
    x = x.apply(pd.to_numeric, errors='ignore', downcast='float')

    # Set patients with Missing = 1 to all NaNs and drop the Missing row
    if missing is not None:
        x.loc[x[missing] == 1, :] = np.nan
        x.drop(missing, axis=1, inplace=True)

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

    trisomies_base = 'mmrf_explore_ia9_cytogenetics_trisomies_'
    trisomies_times = [
        'baseline'
    ]

    # Load baseline FISH data
    abn_fish_raw = pd.read_csv(os.path.join(data_dir, abn_fish_base + abn_fish_times[0] + '.csv'))
    abn_fish = pivot_dataset(abn_fish_raw, 'Abnormalities (FISH) - Baseline', exclude=['Multiple', 'None'], missing='Missing')

    # Load baseline META data
    abn_meta_raw = pd.read_csv(os.path.join(data_dir, abn_meta_base + abn_meta_times[0] + '.csv'))
    abn_meta = pivot_dataset(abn_meta_raw, 'Abnormalities (Meta) - Baseline', exclude=['Other'])

    # May need special handling of "del17 or p53" and "del 17" cols

    # Combine FISH and META. An abnormality is present if it's a 1 in either of the datasets
    abn_comb = merge_datasets_merge_cols(abn_fish, abn_meta)

    # TODO: Grab the Other abnormalities in the META dataset

    # Load trisomies
    trisomies_raw = pd.read_csv(os.path.join(data_dir, trisomies_base + trisomies_times[0] + '.csv'))
    trisomies_raw.drop_duplicates(subset=['Patient', 'Trisomies - Baseline'], take_last=True, inplace=True)  # Delete duplicate info like in row 114 (of the csv)
    trisomies = pivot_dataset(trisomies_raw, 'Trisomies - Baseline', exclude=[], missing='Not Done')

    # Rename trisomies columns
    trisomies.rename(columns=dict([(i, 'Trisomy_' + i) for i in list(trisomies)]), inplace=True)

    # Combine FISH+META and trisomies.
    abn_comb = abn_comb.join(trisomies, how='outer')

    # Save results
    results_file = os.path.join(data_dir, 'data_test' + '.csv')
    abn_comb.to_csv(results_file)

    print(abn_comb)

    1
