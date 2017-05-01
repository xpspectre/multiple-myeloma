# Prepare per-patient clinical data
# TODO: Possibly incorporate Palumbo data

import os
from load_patient_data import load_per_patient_data
import numpy as np
import pandas as pd

pd.options.mode.chained_assignment = None

# Load input data
data, per_patient_dict, per_patient_fields = load_per_patient_data()

# Dictionary of categorical mappings
cats = {}

# Build a cleaned up dataset data
#   Lots of the cols are redundant

# Build a separate dataset of endpoints endp
endp = data[['PUBLIC_ID']]

# Study ID - CoMMpass vs Palumbo
#   For now, we'll just take the CoMMpass patients
study = data[['PUBLIC_ID']]

# Therapy info - this is a simplified version of the stuff in the other table, but is enough to get started with coarse
#   classes
treat = data[['PUBLIC_ID']]

# Keep PUBLIC_ID as the main index

# Drop informed consent date D_PT_ic_day
# Drop reason for dropping study D_PT_disc_com
data.drop(['D_PT_ic_day', 'D_PT_disc_com'], axis=1, inplace=True)

# Study end to endp
endp = endp.join(data['D_PT_lastdy'])
data.drop('D_PT_lastdy', axis=1, inplace=True)

# Nobody has completed the study, they're either still in it or dropped for some reason (including death?)
complete_map = {
    '': np.nan,
    'No': 0
}
cats['D_PT_DIDPATIENTCOM'] = complete_map
endp = endp.join(data['D_PT_DIDPATIENTCOM'].replace(complete_map))
data.drop('D_PT_DIDPATIENTCOM', axis=1, inplace=True)

# Primary reason for patient drop - death vs other
drop_reason_map = {
    '': np.nan,
    'Death': 1,
    'Other': 0,
    'Patient no longer consents to participate in the study': 0,
    'Patient is lost to follow-up': 0,
    'Inter-current illness that interferes with study assessments': 0,
    'Noncompliance with study procedures': 0
}
cats['D_PT_PRIMARYREASON'] = drop_reason_map
endp = endp.join(data['D_PT_PRIMARYREASON'].replace(drop_reason_map))
data.drop('D_PT_PRIMARYREASON', axis=1, inplace=True)

# Cause of death due to MM or other
#   For some analyses, other is a form of right-censoring
death_reason_map = {
    '': np.nan,
    'Disease Progression': 1,
    'Other': 0
}
cats['D_PT_CAUSEOFDEATH'] = death_reason_map
endp = endp.join(data['D_PT_CAUSEOFDEATH'].replace(death_reason_map))
data.drop('D_PT_CAUSEOFDEATH', axis=1, inplace=True)

# Date of death
endp = endp.join(data['D_PT_deathdy'])
data.drop('D_PT_deathdy', axis=1, inplace=True)

# Drop some redundant cols that are hard to interpret
#   A bunch of these are coded versions of more descriptive cols
#   D_PT_trtstdy is just 1 for everyone
data.drop(['D_PT_complete', 'D_PT_discont', 'D_PT_DISCREAS', 'D_PT_dthreas', 'D_PT_raceoth', 'D_PT_race',
           'D_PT_ethnic', 'D_PT_gender', 'D_PT_DIDTHEPATIENT', 'D_PT_screen', 'D_PT_trtstdy', 'D_PT_sdeathdy',
           'D_PT_enr', 'D_PT_lvisit', 'D_PT_lvisitdy', 'D_PT_lvisitc'], axis=1, inplace=True)

# Last day seen alive is important for right-censoring
endp = endp.join(data['D_PT_lstalive'])
data.drop('D_PT_lstalive', axis=1, inplace=True)

# Keep age D_PT_age

# Drop Palumbo stuff for now
data.drop(['CLINICAL', 'RANDOM', 'gender_char', 'race_char', 'informed_consent_version', 'Date_of_diagnosis',
           'ENROLLED'], axis=1, inplace=True)

# Keep apparant Palumbo col but actually applies to everyone: demog_height (in) and demog_weight (lb)

# Drop units for height and weight
data.drop(['DEMOG_HEIGHTUNITOFM', 'DEMOG_WEIGHTUNITOFM'], axis=1, inplace=True)

# ISS disease stage "endpoint": D_PT_iss
endp = endp.join(data['D_PT_iss'])
data.drop('D_PT_iss', axis=1, inplace=True)

# CoMMpass vs Palumbo patients
study_map = {
    'CoMMpass': 1,
    'Palumbo': 0
}
cats['STUDY_ID'] = study_map
study = study.join(data['STUDY_ID'].replace(study_map))
data.drop('STUDY_ID', axis=1, inplace=True)

# Drop redundant stage info
data.drop(['D_PT_issstage_char', 'D_PT_issstage'], axis=1, inplace=True)

# Preprocess basic treatment info
#    3 individual treatments: Bortezomib, Carfilzomib, IMIDs
#    3 dual treatments: bortezomib/carfilzomib, bortezomib/IMIDs, IMIDs/carfilzomib
#    1 triple treatment: bortezomib/IMIDs/carfilzomib
# Make indicators of presence of each treatment from D_PT_therclass col
#   Could use a regex...
has_bor = (data['D_PT_therclass'] == 'Bortezomib-based') | \
          (data['D_PT_therclass'] == 'combined bortezomib/carfilzomib-based') | \
          (data['D_PT_therclass'] == 'combined bortezomib/IMIDs-based') | \
          (data['D_PT_therclass'] == 'combined bortezomib/IMIDs/carfilzomib-based')
has_car = (data['D_PT_therclass'] == 'Carfilzomib-based') | \
          (data['D_PT_therclass'] == 'combined bortezomib/carfilzomib-based') | \
          (data['D_PT_therclass'] == 'combined IMIDs/carfilzomib-based') | \
          (data['D_PT_therclass'] == 'combined bortezomib/IMIDs/carfilzomib-based')
has_imi = (data['D_PT_therclass'] == 'IMIDs-based') | \
          (data['D_PT_therclass'] == 'combined bortezomib/IMIDs-based') | \
          (data['D_PT_therclass'] == 'combined IMIDs/carfilzomib-based') | \
          (data['D_PT_therclass'] == 'combined bortezomib/IMIDs/carfilzomib-based')
treat['TREAT_BOR'] = has_bor.astype(int)  # True/False -> 1/0 map
treat['TREAT_CAR'] = has_car.astype(int)
treat['TREAT_IMI'] = has_imi.astype(int)

# Drop the rest of the treatment cols
data.drop(['D_PT_therclass', 'D_PT_therfstn', 'D_PT_therclassn', 'D_PT_maxline', 'ftrttrpl'], axis=1, inplace=True)

# Copy over the SCT (stem cell transplant) codes, but I don't know what they are
treat = treat.join(data[['sct_bresp', 'line1sct']])
data.drop(['sct_bresp', 'line1sct'], axis=1, inplace=True)

# What is PD? It sounds like a response, so move to endp table
endp = endp.join(data[['D_PT_pddy', 'D_PT_pdflag', 'D_PT_ttfpdw', 'D_PT_respdur', 'D_PT_mmstatus', 'D_PT_mmstatus1', 'D_PT_mmstatus2', 'D_PT_mmstatus3', 'D_PT_rapd', 'D_PT_dresp']])
data.drop(['D_PT_pddy', 'D_PT_pdflag', 'D_PT_ttfpdw', 'D_PT_respdur', 'D_PT_mmstatus', 'D_PT_mmstatus1', 'D_PT_mmstatus2', 'D_PT_mmstatus3', 'D_PT_rapd', 'D_PT_dresp'], axis=1, inplace=True)

# Drop redundant cols
data.drop(['demog_vj_interval', 'demog_visitdy'], axis=1, inplace=True)

# Keep race cols, lose other+unknown cols to get a linearly independent categorical set
# Keep DEMOG_AMERICANINDIA, DEMOG_BLACKORAFRICA, DEMOG_NATIVEHAWAIIA, DEMOG_WHITE, DEMOG_ASIAN and convert checked
checked_map = {
    '': 0,
    'Checked': 1
}
data['DEMOG_AMERICANINDIA'] = data['DEMOG_AMERICANINDIA'].replace(checked_map)
data['DEMOG_BLACKORAFRICA'] = data['DEMOG_BLACKORAFRICA'].replace(checked_map)
data['DEMOG_NATIVEHAWAIIA'] = data['DEMOG_NATIVEHAWAIIA'].replace(checked_map)
data['DEMOG_WHITE'] = data['DEMOG_WHITE'].replace(checked_map)
data['DEMOG_ASIAN'] = data['DEMOG_ASIAN'].replace(checked_map)

data.drop(['DEMOG_OTHER', 'DEMOG_SPECIFY'], axis=1, inplace=True)

# Gender - use this col since we know/control the coding
gender_map = {
    'Male': 1,
    'Female': 0,
    '': np.nan
}
cats['DEMOG_GENDER'] = gender_map
data['DEMOG_GENDER'] = data['DEMOG_GENDER'].replace(gender_map)

# Ethnicity: Hispanic/Latino or not
eth_map = {
    'Hispanic or Latino': 1,
    'Not Hispanic or Latino': 0,
    'Other': 0,
    '': 0
}
cats['DEMOG_ETHNICITY'] = eth_map
data['DEMOG_ETHNICITY'] = data['DEMOG_ETHNICITY'].replace(eth_map)
data.drop(['DEMOG_SPECIFY2'], axis=1, inplace=True)

# Drop redundant visit and age cols
data.drop(['DEMOG_DAYOFVISIT', 'DEMOG_DAYOFBIRTH', 'DEMOG_PATIENTAGE', 'demog_visit', 'enr'], axis=1, inplace=True)

# print(data['DEMOG_ETHNICITY'].unique())
# print(data)
# print(treat)

# Save processed tables
output_dir = 'data/processed'

data.set_index('PUBLIC_ID', inplace=True)
data.to_csv(os.path.join(output_dir, 'patient_data.csv'))

endp.set_index('PUBLIC_ID', inplace=True)
endp.to_csv(os.path.join(output_dir, 'patient_endp.csv'))

treat.set_index('PUBLIC_ID', inplace=True)
treat.to_csv(os.path.join(output_dir, 'patient_treat.csv'))

