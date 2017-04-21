# Load various patient data, combine, preprocess
# Note/Warning: Are these files latin1 encoded?

import os
import numpy as np
import pandas as pd

# Global constants
data_dir = 'data/raw/clinical_data_tables'
flat_files_dir = 'CoMMpass_IA9_FlatFiles'
flat_files_dicts_dir = 'CoMMpass_IA9_FlatFile_Dictionaries'


def load_per_patient_data():
    """Load main patient info file"""
    dict_file = os.path.join(data_dir, flat_files_dicts_dir, 'PER_PATIENT.xlsx')
    dict_data = pd.read_excel(dict_file)

    fields = {}
    for i, row in dict_data.iterrows():
        type_str = row['vartype']
        if type_str == 'Char':
            type = str
        elif type_str == 'Num':
            type = float
        else:
            raise ValueError('Unknown type {type_str} detected. Fix this.'.format(type_str=type_str))

        fields[row['name']] = type

    data_file = os.path.join(data_dir, flat_files_dir, 'PER_PATIENT.csv')
    data = pd.read_csv(data_file, engine='c', dtype=fields)
    data.set_index('PUBLIC_ID')

    # Fill in NaN cells in str columns with blank strings
    for label, type in fields.items():
        if type == str:
            data[label].fillna('', inplace=True)

    return data, dict_data, fields


def load_per_visit_data():
    """Load individual visit data"""
    dict_file = os.path.join(data_dir, flat_files_dicts_dir, 'PER_PATIENT_VISIT.xlsx')
    dict_data = pd.read_excel(dict_file)

    fields = {}
    for i, row in dict_data.iterrows():
        label = row['name']
        type_str = row['vartype']
        if type_str == 'Char':
            type = str
        elif type_str == 'Num':
            type = np.float64
        else:
            raise ValueError('Unknown type {type_str} detected. Fix this.'.format(type_str=type_str))

        # Handle stupid corner case for the D_CM_ABNORMALPERCE col, which is numeric but contains 2 rows with '**'
        if label == 'D_CM_ABNORMALPERCE':
            type = str

        fields[label] = type

    data_file = os.path.join(data_dir, flat_files_dir, 'PER_PATIENT_VISIT.csv')
    data = pd.read_csv(data_file, engine='c', encoding='latin1', dtype=fields,
                                 warn_bad_lines=True)
    data.set_index('PUBLIC_ID')

    # Fix D_CM_ABNORMALPERCE
    data['D_CM_ABNORMALPERCE'] = pd.to_numeric(data['D_CM_ABNORMALPERCE'], errors='coerce')
    fields['D_CM_ABNORMALPERCE'] = np.float64

    # Fill in NaN cells in str columns with blank strings
    for label, type in fields.items():
        if type == str:
            data[label].fillna('', inplace=True)

    return data, dict_data, fields


def load_treatment_regimen():
    dict_file = os.path.join(data_dir, flat_files_dicts_dir, 'STAND_ALONE_TREATMENT_REGIMEN.xlsx')
    dict_data = pd.read_excel(dict_file)

    fields = {}
    for i, row in dict_data.iterrows():
        type_str = row['vartype']
        if type_str == 'Char':
            type = str
        elif type_str == 'Num':
            type = float
        else:
            raise ValueError('Unknown type {type_str} detected. Fix this.'.format(type_str=type_str))

        fields[row['name']] = type

    data_file = os.path.join(data_dir, flat_files_dir, 'STAND_ALONE_TREATMENT_REGIMEN.csv')
    data = pd.read_csv(data_file, engine='c', dtype=fields)

    # Rename public_id field, since it's lowercase for some reason in this file
    data.rename(columns={'public_id': 'PUBLIC_ID'}, inplace=True)
    fields['PUBLIC_ID'] = fields.pop('public_id')

    data.set_index('PUBLIC_ID')

    # Fill in NaN cells in str columns with blank strings
    for label, type in fields.items():
        if type == str:
            data[label].fillna('', inplace=True)

    return data, dict_data, fields


def load_treatment_resp():
    dict_file = os.path.join(data_dir, flat_files_dicts_dir, 'STAND_ALONE_TRTRESP.xlsx')
    dict_data = pd.read_excel(dict_file)

    fields = {}
    for i, row in dict_data.iterrows():
        type_str = row['vartype']
        if type_str == 'Char':
            type = str
        elif type_str == 'Num':
            type = float
        else:
            raise ValueError('Unknown type {type_str} detected. Fix this.'.format(type_str=type_str))

        fields[row['name']] = type

    data_file = os.path.join(data_dir, flat_files_dir, 'STAND_ALONE_TRTRESP.csv')
    data = pd.read_csv(data_file, engine='c', dtype=fields)

    # Rename public_id field, since it's lowercase for some reason in this file
    data.rename(columns={'public_id': 'PUBLIC_ID'}, inplace=True)
    fields['PUBLIC_ID'] = fields.pop('public_id')

    data.set_index('PUBLIC_ID')

    # Fill in NaN cells in str columns with blank strings
    for label, type in fields.items():
        if type == str:
            data[label].fillna('', inplace=True)

    return data, dict_data, fields


def load_survival():
    dict_file = os.path.join(data_dir, flat_files_dicts_dir, 'STAND_ALONE_SURVIVAL.xlsx')
    dict_data = pd.read_excel(dict_file)

    fields = {}
    for i, row in dict_data.iterrows():
        type_str = row['vartype']
        if type_str == 'Char':
            type = str
        elif type_str == 'Num':
            type = float
        else:
            raise ValueError('Unknown type {type_str} detected. Fix this.'.format(type_str=type_str))

        fields[row['name']] = type

    data_file = os.path.join(data_dir, flat_files_dir, 'STAND_ALONE_SURVIVAL.csv')
    data = pd.read_csv(data_file, engine='c', dtype=fields)

    # Rename public_id field, since it's lowercase for some reason in this file
    data.rename(columns={'public_id': 'PUBLIC_ID'}, inplace=True)
    fields['PUBLIC_ID'] = fields.pop('public_id')

    data.set_index('PUBLIC_ID')

    # Fill in NaN cells in str columns with blank strings
    for label, type in fields.items():
        if type == str:
            data[label].fillna('', inplace=True)

    return data, dict_data, fields


if __name__ == '__main__':
    per_patient_data, per_patient_dict, per_patient_fields = load_per_patient_data()
    per_visit_data, per_visit_dict, per_visit_fields = load_per_visit_data()
    treatment_data, treatment_dict, treatment_fields = load_treatment_regimen()
    treatment_resp_data, treatment_resp_dict, treatment_resp_fields = load_treatment_resp()
    survival_data, survival_dict, survival_fields = load_survival()

    # Grab all info on 1 patient
    patient = 'MMRF_1014'
    demographics = per_patient_data.loc[per_patient_data['PUBLIC_ID'] == patient]
    data = per_visit_data.loc[per_visit_data['PUBLIC_ID'] == patient]
    treatment = treatment_data.loc[treatment_data['PUBLIC_ID'] == patient]
    treatment_resp = treatment_resp_data.loc[treatment_resp_data['PUBLIC_ID'] == patient]
    survival = survival_data.loc[survival_data['PUBLIC_ID'] == patient]

    print(demographics)
    print(data)
    print(treatment)
    print(treatment_resp)
    print(survival)

    1
