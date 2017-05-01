# Main function that goes thru PER_PATIENT_VISIT clinical data table and turns it into a usable form
# Goes thru columns manually and assess quality
# A unique index is made from PUBLIC_ID + VISIT
# TODO: Apply cols in pres to turn various rows of data into nan (especially those converted from 'Checked'/blank)
#   This is done for data table only. If we want to look at other tables, apply this to them as well
# Note: This dataset has 2 cytogenetics blocks: conventional and fluorescence
#   Combine them somehow

from load_patient_data import load_per_visit_data
import numpy as np
import pandas as pd
import os
import cmsgpack


yn_map = {
    'Yes': 1,
    'No': 0,
    '': np.nan
}

# Check marks indicate yes, blanks indicate no, but note that other metadata cols may indicate these should be NaN
checked_map = {
    'Checked': 1,
    '': 0
}


# Helper functions
def replace_map(x, col, map):
    """General function to replace col vals with map"""
    d = x[col].replace(map)
    x.drop(col, axis=1, inplace=True)
    x[col] = d
    return x


def replace_yn(x, cols):
    """Replace a yes/no/blank with 1/0/nan col. Returns x."""
    for col in cols:
        d = x[col].replace(yn_map)
        x.drop(col, axis=1, inplace=True)
        x[col] = d
    return x


def replace_checked(x, cols):
    """Replace checked/blank with 1/0. Returns x"""
    for col in cols:
        d = x[col].replace(checked_map)
        x.drop(col, axis=1, inplace=True)
        x[col] = d
    return x


def add_checked_cols(x, cols):
    """Convert list of cols that contains 'Checked' or '' (blank) into 1 or 0 and add to x. Returns x.
    Note: this is a closure that depends on data"""
    for col in cols:
        x = x.join(data[col].replace(checked_map))
    return x


# Load input data
data, per_visit_dict, per_visit_fields = load_per_visit_data()

# Keep track of categorical cols in data
#   key = col name, val = dict of mapping
cats = {}

# Modify data df in place as main data store
# Build pres df to hold cols solely used to indicate whether something was done or not
# Note: Lots of the cols specify yes/no using checked/blank. Assume checked = yes, blank = no. These should be in data, not in pres
#   0 = No (explicit no or blank, indicating unknown/not applicable which is equiv to no here), 1 = Yes
pres = data[['PUBLIC_ID', 'VISIT']]

# Dates - assessment dates that may be different from visit dates - may use them or just coarsen to use visit date
date = data[['PUBLIC_ID', 'VISIT']]

# Misc data like assessment dates that may be used, categories that provide additional context, etc.
# The idea is to put stuff here that we don't know how to put into the main feature tables yet
misc = data[['PUBLIC_ID', 'VISIT']]

# Text cols of various types that should be followed up on
#   Ex: Radiology image interpretations. These are short and we can get basic sentiments like normal, insignificant,
#   significant, abnormal, significant but unrelated to MM etc. from simple text parsing
text = data[['PUBLIC_ID', 'VISIT']]

# Treatment data in the per patient visit dataset
#   Includes stem cell transplant status and whether patient is eligible for various treatments
treat = data[['PUBLIC_ID', 'VISIT']]

# Endpoints like death
endp = data[['PUBLIC_ID', 'VISIT']]

# Survey questions: bare q's and aggregated results
survey = data[['PUBLIC_ID', 'VISIT']]
survey_agg = data[['PUBLIC_ID', 'VISIT']]

########################################################################################################################
# Keep PUBLIC_ID as main index

# Drop spectrum id
data.drop('SPECTRUM_SEQ', axis=1, inplace=True)

# Keep visit day but drop coarse-grained visit interval
data.drop('VJ_INTERVAL', axis=1, inplace=True)

# Drop informed consent day
data.drop('VISITDY', axis=1, inplace=True)

########################################################################################################################

# Bone lesions check: Either 'checked' -> 1 or blank -> NaN
pres = pres.join(data['AT_DEFINITEDEVEL'].replace(yn_map))
data.drop('AT_DEFINITEDEVEL', axis=1, inplace=True)

# Treatment response checked at last visit
pres = pres.join(data['AT_WASANASSESSME'].replace(yn_map))
data.drop('AT_WASANASSESSME', axis=1, inplace=True)

# Response assessment date
date = date.join(data['AT_RESPONSEASSES'])
data.drop('AT_RESPONSEASSES', axis=1, inplace=True)

# Convert treatment response to numeric index and put into end results table
response_map = {
    'Stringent Complete Response (sCR)': 6,
    'Complete Response': 5,
    'Very Good Partial Response (VGPR)': 4,
    'Partial Response': 3,
    'Stable Disease': 2,
    'Progressive Disease': 1,
    '': np.nan
}
endp = endp.join(data['AT_TREATMENTRESP'].replace(response_map))
data.drop('AT_TREATMENTRESP', axis=1, inplace=True)
cats['AT_TREATMENTRESP'] = response_map

# Basic checks on worsening conditions
data = replace_checked(data, ['AT_INCREASEOF25F', 'AT_SERUMMCOMPONE', 'AT_URINEMCOMPONE', 'AT_ONLYINPATIENT', 'AT_ONLYINPATIENT2', 'AT_DEVELOPMENTOF'])

# Go back and ensure:
#   NaN to unassessed patients
assessed = pres['AT_WASANASSESSME'] == 1
replace_cols = ['AT_DEFINITEDEVEL', 'AT_TREATMENTRESP', 'AT_INCREASEOF25F', 'AT_SERUMMCOMPONE', 'AT_URINEMCOMPONE', 'AT_ONLYINPATIENT', 'AT_ONLYINPATIENT2', 'AT_DEVELOPMENTOF']
for col in replace_cols:
    data.loc[assessed == 0, col] = np.nan

########################################################################################################################
# Bone marrow aspriate useful?
pres = pres.join(data['BMA_WASTHESCREENI'].replace(yn_map))
data.drop('BMA_WASTHESCREENI', axis=1, inplace=True)

########################################################################################################################
# Bone marrow assessment done?
pres = pres.join(data['BA_WASABONEASSES'].replace(yn_map))
data.drop('BA_WASABONEASSES', axis=1, inplace=True)

# Bone assessment fields
# Keep bone assessment date
date = date.join(data['BA_DAYOFASSESSM'])
data.drop('BA_DAYOFASSESSM', axis=1, inplace=True)

# Type of bone assessment - don't know if it will be useful right now
misc = misc.join(data[['BA_TYPEOFBONEASS', 'BA_SPECIFY']])
data.drop(['BA_TYPEOFBONEASS', 'BA_SPECIFY'], axis=1, inplace=True)

# Assign bone interpretation categories
bone_interp_map = {
    'Clinically significant abnormality related to MM': 2,
    'Clinically significant abnormality unrelated to MM': 1,
    'Normal or clinically insignificant abnormality': 0,
    '': np.nan,
}
data = replace_map(data, 'BA_IMAGEINTERPRE', bone_interp_map)
cats['BA_IMAGEINTERPRE'] = bone_interp_map

# Lytic lesions
data = replace_checked(data, ['BA_LYTICLESIONS'])
lesions_map = {
    '': np.nan,
    '1': 1,
    '2': 2,
    '>=3': 3
}
data = replace_map(data, 'BA_OFLYTICLESION', lesions_map)
cats['BA_OFLYTICLESION'] = lesions_map

# Checks on a bunch of bones
data = replace_checked(data, ['BA_OSTEOPENIAOST', 'BA_PATHOLOGICFRA', 'BA_MULTIPLEDISSE', 'BA_SKULL', 'BA_SPINE', 'BA_CERVICAL', 'BA_THORACIC', 'BA_LUMBAR', 'BA_SACRAL', 'BA_PELVIS', 'BA_RIGHT', 'BA_LEFT', 'BA_FEMUR', 'BA_RIGHT2', 'BA_LEFT2', 'BA_HUMERUS', 'BA_RIGHT3', 'BA_LEFT3', 'BA_RIBS', 'BA_RIGHT4', 'BA_LEFT4', 'BA_OTHER'])
# Keep number of ribs on each side affected: BA_NUMBERAFFECTE and BA_NUMBERAFFECTE2

# Text data on other myeloma-related bone damage
text = text.join(data['BA_SPECIFY2'])
data.drop('BA_SPECIFY2', axis=1, inplace=True)

# Go back and ensure:
#   0 to assessed patients with no abnormality
#   NaN to unassessed patients
assessed = pres['BA_WASABONEASSES'] == 1
no_abnormal = data['BA_IMAGEINTERPRE'] == 0
replace_cols = ['BA_LYTICLESIONS', 'BA_OFLYTICLESION', 'BA_OSTEOPENIAOST',  'BA_PATHOLOGICFRA', 'BA_MULTIPLEDISSE', 'BA_SKULL', 'BA_SPINE', 'BA_CERVICAL', 'BA_THORACIC', 'BA_LUMBAR', 'BA_SACRAL', 'BA_PELVIS', 'BA_RIGHT', 'BA_LEFT', 'BA_FEMUR', 'BA_RIGHT2', 'BA_LEFT2', 'BA_HUMERUS', 'BA_RIGHT3', 'BA_LEFT3', 'BA_RIBS', 'BA_RIGHT4', 'BA_LEFT4', 'BA_OTHER', 'BA_NUMBERAFFECTE', 'BA_NUMBERAFFECTE2']
for col in replace_cols:
    data.loc[assessed == 0, col] = np.nan  # does the work
    data.loc[no_abnormal == 1, col] = 0  # shouldn't be necessary

########################################################################################################################
# Other reason for procedure. AML seems significant and weird but others seem like normal
misc = misc.join(data[['BB_SPECIFY', 'BB_REASONFORPROC']])
data.drop(['BB_SPECIFY', 'BB_REASONFORPROC'], axis=1, inplace=True)

# Handle cellularity cols
# http://stackoverflow.com/questions/16689514/how-to-get-the-average-of-dataframe-column-values
#   If both are present, average them
#   If 1 is present, take it
#   If neither are present, leave as nan
data['BB_CELLULARITY'] = data[['BB_CELLULARITY', 'BB_CELLULARITY2']].mean(axis=1)
data.drop('BB_CELLULARITY2', axis=1, inplace=True)

# Fix negative pct plasma cells with nans
data.loc[(data['BB_PERCENTOFPLAS'] < 0) , 'BB_PERCENTOFPLAS'] = np.nan

# Bone marrow biopsy collected?
pres = pres.join(data['BB_WASABONEMARRO'].replace(yn_map))
data.drop('BB_WASABONEMARRO', axis=1, inplace=True)
pres = pres.join(data['BB_WASBIOPSYSPEC'].replace(yn_map))
data.drop('BB_WASBIOPSYSPEC', axis=1, inplace=True)

# Bone marrow biopsy date of collection
date = date.join(data['BB_DAYOFCOLLECT'])
data.drop('BB_DAYOFCOLLECT', axis=1, inplace=True)

# Reason for invalid bone biopsy/aspirate
misc = misc.join(data[['BB_IFNOPLEASESPE']])
data.drop(['BB_IFNOPLEASESPE'], axis=1, inplace=True)

# Go back and ensure:
#   NaN to unassessed and invalid results patients
assessed = (pres['BB_WASABONEMARRO'] == 1) & (pres['BB_WASBIOPSYSPEC'] == 1)
replace_cols = ['BB_CELLULARITY', 'BB_PERCENTOFPLAS']
for col in replace_cols:
    data.loc[assessed == 0, col] = np.nan

########################################################################################################################
misc = misc.join(data[['BONE_SPECIFY', 'BONE_IFNOPLEASESPE', 'BONE_IFNOPLEASESPE2']])
data.drop(['BONE_SPECIFY', 'BONE_IFNOPLEASESPE', 'BONE_IFNOPLEASESPE2'], axis=1, inplace=True)

pres = pres.join(data['BONE_WASTHEASPIRAT'].replace(yn_map))
data.drop('BONE_WASTHEASPIRAT', axis=1, inplace=True)
pres = pres.join(data['BONE_WASABONEMARRO'].replace(yn_map))
data.drop('BONE_WASABONEMARRO', axis=1, inplace=True)
date = date.join(data['BONE_DAYOFCOLLECT'])
data.drop('BONE_DAYOFCOLLECT', axis=1, inplace=True)
pres = pres.join(data['BONE_WASBONEMARROW'].replace(yn_map))  # Warning, this col looks like the one 2 above it
data.drop('BONE_WASBONEMARROW', axis=1, inplace=True)

# Bone biopsy metadata
misc = misc.join(data[['BONE_VOLUMEOFBMASP', 'BONE_PLASMAUNK', 'BONE_BMAVOLUMESENT', 'BONE_BMAUNK', 'BONE_WASPERIPHERAL']])
data.drop(['BONE_VOLUMEOFBMASP', 'BONE_PLASMAUNK', 'BONE_BMAVOLUMESENT', 'BONE_BMAUNK', 'BONE_WASPERIPHERAL'], axis=1, inplace=True)
date = date.join(data['BONE_DAYOFCOLLECT2'])
data.drop('BONE_DAYOFCOLLECT2', axis=1, inplace=True)

# Note: For the next few sections, the "reported as %PC" cols go into pres - they're meta-features
# %PC positive for CD138 by IHC checked
pres = add_checked_cols(pres, ['BONE_PC6'])
data.drop('BONE_PC6', axis=1, inplace=True)
# Keep measurement of BONE_PERCENTOFPLAS

# Reason for bone procedure
misc = misc.join(data[['BONE_REASONFORPROC', 'BONE_SPECIFY2']])
data.drop(['BONE_REASONFORPROC', 'BONE_SPECIFY2'], axis=1, inplace=True)

# Plasma cells in bone marrow aspirate
data = replace_checked(data, ['BONE_PERCENTOFPLAS3'])

# Lab for bone marrow aspirate
misc = misc.join(data['BONE_PERCENTOFPLAS2'])
data.drop('BONE_PERCENTOFPLAS2', axis=1, inplace=True)

# BMA: Reported as %PC by cytomorphology
pres = add_checked_cols(pres, ['BONE_REPORTEDASPER'])
data.drop('BONE_REPORTEDASPER', axis=1, inplace=True)

# BMA: %PC positive for CD138 by IHC
pres = add_checked_cols(pres, ['BONE_RANGEOFPLASMA2', 'BONE_RANGEUNK'])
data.drop(['BONE_RANGEOFPLASMA2', 'BONE_RANGEUNK'], axis=1, inplace=True)
# Keep lo and hi measurements BONE_RANGELOW and BONE_RANGEHIGH

# BMA: %PC by cytomorphology (smear)
pres = add_checked_cols(pres, ['BONE_BCYTOMORPHOLO', 'BONE_PC7', 'BONE_RANGEOFPLASMA3', 'BONE_FCAUNK'])
data.drop(['BONE_BCYTOMORPHOLO', 'BONE_PC7', 'BONE_RANGEOFPLASMA3', 'BONE_FCAUNK'], axis=1, inplace=True)
# Keep actual val BONE_PC2
# Keep lo and hi measurements BONE_RANGELOW2 and BONE_RANGEHIGH2

# BMA: Reported as %cells positive by flow cytometry(FLOW)
pres = add_checked_cols(pres, ['BONE_AFLOWCYTOMETR', 'BONE_PC8', 'BONE_RANGEOFPLASMA4', 'BONE_UNKNOWN'])
data.drop(['BONE_AFLOWCYTOMETR', 'BONE_PC8', 'BONE_RANGEOFPLASMA4', 'BONE_UNKNOWN'], axis=1, inplace=True)
# Keep actual val BONE_REPORTEDBYCYT
# Keep lo and hi measurements BONE_RANGELOW3 and BONE_RANGEHIGH3

# Plasma cells reported in a clot section or bone marrow biopsy: (%PC positive for CD138 by IHC)
data = replace_checked(data, ['BONE_PC'])
pres = add_checked_cols(pres, ['BONE_REPORTEDASPER2', 'BONE_PC9', 'BONE_RANGEOFPLASMA5', 'BONE_LOWVALUE'])
data.drop(['BONE_REPORTEDASPER2', 'BONE_PC9', 'BONE_RANGEOFPLASMA5', 'BONE_LOWVALUE'], axis=1, inplace=True)
# Keep actual val BONE_PC3
# Keep lo and hi measurements BONE_RANGELOW4 and BONE_RANGEHIGH4

# Lab for clot section
misc = misc.join(data['BONE_LOCALLABNAME'])
data.drop('BONE_LOCALLABNAME', axis=1, inplace=True)

# PC in clot section or BM biopsy: Reported as %PC by cytomorphology (smear)
pres = add_checked_cols(pres, ['BONE_REPORTEDASPER3', 'BONE_PC10', 'BONE_RANGEOFPLASMA6', 'BONE_HIGHVALUE'])
data.drop(['BONE_REPORTEDASPER3', 'BONE_PC10', 'BONE_RANGEOFPLASMA6', 'BONE_HIGHVALUE'], axis=1, inplace=True)
# Keep actual val BONE_PC4
# Keep lo and hi measurements BONE_RANGELOW5 and BONE_RANGEHIGH5

# PC in clot section or BM biopsy: Reported as %cells positive by flow cytometry(FLOW)
pres = add_checked_cols(pres, ['BONE_RANGEOFPLASMA', 'BONE_PC11', 'BONE_RANGEOFPLASMA7', 'BONE_NOTDONEUNKNOW'])
data.drop(['BONE_RANGEOFPLASMA', 'BONE_PC11', 'BONE_RANGEOFPLASMA7', 'BONE_NOTDONEUNKNOW'], axis=1, inplace=True)
# Keep actual val BONE_PC5
# Keep lo and hi measurements BONE_RANGELOW6 and BONE_RANGEHIGH6

misc = misc.join(data['BONE_DOESTHEPCRESU'])
data.drop('BONE_DOESTHEPCRESU', axis=1, inplace=True)

# Bone specimen adequate for assessment? May have other implications...
pres = pres.join(data['BONE_WASTHESPECIME'].replace(yn_map))
data.drop('BONE_WASTHESPECIME', axis=1, inplace=True)

# Go back and ensure:
#   NaN to unassessed patients for significant cols handled by replace_checked
assessed = pres['BONE_WASABONEMARRO'] == 1
replace_cols = ['BONE_PERCENTOFPLAS3', 'BONE_PC']
for col in replace_cols:
    data.loc[assessed == 0, col] = np.nan

########################################################################################################################
# Bone marrow transplant
treat = treat.join(data['BMT_WASABONEMARRO'].replace(yn_map))
data.drop('BMT_WASABONEMARRO', axis=1, inplace=True)
treat = treat.join(data['BMT_DAYOFTRANSPL'])
data.drop('BMT_DAYOFTRANSPL', axis=1, inplace=True)
treat = add_checked_cols(treat, ['BMT_AUTOLOGOUS', 'BMT_ALLOGENIC', 'BMT_STEMCELLS', 'BMT_BONEMARROW'])
data.drop(['BMT_AUTOLOGOUS', 'BMT_ALLOGENIC', 'BMT_STEMCELLS', 'BMT_BONEMARROW'], axis=1, inplace=True)
treat = treat.join(data['BMT_PLEASESPECIFY'])
data.drop('BMT_PLEASESPECIFY', axis=1, inplace=True)
data.drop('BMT_SPECIFY', axis=1, inplace=True)  # All empty
treat = treat.join(data['BMT_NUMBER'])
data.drop('BMT_NUMBER', axis=1, inplace=True)

########################################################################################################################
# ECOG: patient function assessments
ecog_map = {
    '': np.nan,
    '0 = Fully Active': 0,
    '1 = Restricted in physically strenuous activity': 1,
    '2 = Ambulatory and capable of all selfcare': 2,
    '3 = Capable of only limited selfcare': 3,
    '4 = Completely disabled': 4,
    '5 = Dead': 5
}
data = replace_map(data, 'ECOG_PERFORMANCEST', ecog_map)
cats['ECOG_PERFORMANCEST'] = ecog_map

pres = pres.join(data['ECOG_WASANECOGASSE'].replace(yn_map))
data.drop('ECOG_WASANECOGASSE', axis=1, inplace=True)

date = date.join(data['ECOG_ASSESSMENTDAT'])
data.drop('ECOG_ASSESSMENTDAT', axis=1, inplace=True)

# ECOG_PERFORMANCEST col contains number or NaN as appropriate if the test was taken

########################################################################################################################
# Drop informed consent cols
data.drop(['IC_DAYOFINFORME', 'IC_DAYOFVISIT', 'IC_DAYUNIVERSAL', 'IC_HASTHEPATIENT', 'IC_IFYESPLEASEPR', 'IC_IPERMITMYSAMP', 'IC_IPERMITMYSAMP2', 'IC_IPERMITRESEAR', 'IC_PATIENTHASREA', 'IC_WASTHISPATIEN2'], axis=1, inplace=True)
# Drop a couple more cols that are probably not important...
data.drop(['IC_DIDTHEPATIENT', 'IC_OTHREASON'], axis=1, inplace=True)
# All the patients are over 18 yo
data.drop(['IC_PATIENTISATLE'], axis=1, inplace=True)

# Keep the patient support system enrollment question
data = replace_yn(data, ['IC_PARTICIPATEPSS'])

# Keep the important cols mixed in with them...
#   Last 2 of these are elevated/not?
data = replace_yn(data, ['IC_INVOLVEDFREEL', 'IC_PATIENTHADANO', 'IC_SERUMMPROTEIN', 'IC_URINEMPROTEIN'])

treat = treat.join(data['IC_NOMORETHAN30D'].replace(yn_map))
data.drop('IC_NOMORETHAN30D', axis=1, inplace=True)

treat = treat.join(data['IC_PATIENTISALRE'].replace(yn_map))
data.drop('IC_PATIENTISALRE', axis=1, inplace=True)

# All of these are no or unknown
treat = treat.join(data['IC_PATIENTISENRO'].replace(yn_map))
data.drop('IC_PATIENTISENRO', axis=1, inplace=True)

# Drop - all empty
data.drop('IC_SPECREASON', axis=1, inplace=True)

treat = treat.join(data['IC_THEPATIENTISA'].replace(yn_map))
data.drop('IC_THEPATIENTISA', axis=1, inplace=True)

########################################################################################################################
# Other surgery
treat = treat.join(data['MMSURG_SPECIFY'])
data.drop('MMSURG_SPECIFY', axis=1, inplace=True)
treat = treat.join(data[['MMSURG_NONE', 'MMSURG_NONE', 'MMSURG_VERTEBROPLAST', 'MMSURG_KYPHOPLASTY', 'MMSURG_OTHER']].replace(yn_map))
data.drop(['MMSURG_NONE', 'MMSURG_NONE', 'MMSURG_VERTEBROPLAST', 'MMSURG_KYPHOPLASTY', 'MMSURG_OTHER'], axis=1, inplace=True)

########################################################################################################################
# Peripheral neuropathy
pres = pres.join(data['PN_WASANEUROLOGI'].replace(yn_map))
data.drop('PN_WASANEUROLOGI', axis=1, inplace=True)
date = date.join(data['PN_ASSESSMENTDAT'])
data.drop('PN_ASSESSMENTDAT', axis=1, inplace=True)

pn_present_map = {
    '': np.nan,
    'Not Done': np.nan,
    'No': 0,
    'Yes': 1
}
data = replace_map(data, 'PN_PERIPHERALSEN', pn_present_map)
cats['PN_PERIPHERALSEN'] = pn_present_map

pn_grade_map = {
    '': np.nan,
    'Grade 1: Asymptomatic; loss of deep tendon reflexes or parathesia.': 1,
    'Grade 2: Moderate symptoms; limiting instrumental ADL.': 2,
    'Grade 3: Severe symptoms; limiting self care ADL.': 3
}
data = replace_map(data, 'PN_IFYESINDICATE', pn_grade_map)
cats['PN_IFYESINDICATE'] = pn_grade_map

data = replace_map(data, 'PN_PERIPHERALMOT', pn_present_map)  # same map as base peripheral neuropathy present
cats['PN_PERIPHERALMOT'] = pn_present_map

pmn_grade_map = {
    '': np.nan,
    'Grade 1: Asymptomatic; clinical or diagnostic observations only; intervention not indicated': 1,
    'Grade 2: Moderate symptoms; limiting instrumental ADL': 2,
    'Grade 3: Severe symptoms; limiting self care ADL; assistive device indicated': 3
}
data = replace_map(data, 'PN_IFYESINDICATE2', pmn_grade_map)
cats['PN_IFYESINDICATE2'] = pmn_grade_map

# These cols contain number or NaN as appropriate if the test was taken

########################################################################################################################
# Emergency room visits since last study visit - keep these as features
data = replace_yn(data, ['RU_DIDTHEPATIENT2', 'RU_DIDTHEPATIENT'])

# Death and endpoint
death_map = {
    '': np.nan,
    'Other': 0,
    'Disease Progression': 1
}
endp = endp.join(data['SE_CAUSEOFDEATH'].replace(death_map))
cats['SE_CAUSEOFDEATH'] = death_map
endp = endp.join(data[['SE_DAYOFDEATH', 'SE_DAYPATIENTCO', 'SE_LASTKNOWNDAY']])
endp = endp.join(data['SE_DIDPATIENTCOM'].replace(yn_map))
data['SE_PRIMARYREASON'] = data['SE_PRIMARYREASON'] + ':' + data['SE_SPECIFY'] + ':' + data['SE_SPECIFY2']  # Join all-cause leaving trial
endp = endp.join(data['SE_PRIMARYREASON'])
data.drop(['SE_CAUSEOFDEATH', 'SE_DAYOFDEATH', 'SE_DAYPATIENTCO', 'SE_DIDPATIENTCOM', 'SE_LASTKNOWNDAY', 'SE_PRIMARYREASON', 'SE_SPECIFY', 'SE_SPECIFY2'], axis=1, inplace=True)

########################################################################################################################
# Supplemental treatments
supp_treatments = ['SUPP_ANTICOAGULATI', 'SUPP_NONE', 'SUPP_TRANSFUSION', 'SUPP_PACKEDREDBLOO', 'SUPP_PLATELETS', 'SUPP_OTHER', 'SUPP_RADIOTHERAPY', 'SUPP_DIALYSIS', 'SUPP_BISPHOSPHONAT', 'SUPP_WBCGROWTHFACT', 'SUPP_ERYTHROPOIESI', 'SUPP_ANTIEMETIC', 'SUPP_ANALGESICSFOR', 'SUPP_OPIOID', 'SUPP_OTHER2', 'SUPP_ANTIVIRAL', 'SUPP_MEDICATIONSFO', 'SUPP_MEDICATIONSFO2']
treat = treat.join(data[supp_treatments].replace(yn_map))
data.drop(supp_treatments, axis=1, inplace=True)
treat = treat.join(data[['SUPP_SPECIFY', 'SUPP_SPECIFYSITE']])
data.drop(['SUPP_SPECIFY', 'SUPP_SPECIFYSITE'], axis=1, inplace=True)

########################################################################################################################
# Myeloma-specific related symptoms
#   SS_ISTHEPATIENTR reports follow-up visits. This block includes follow-up visits and baseline.
#   A good col to indicate whether results should be present is SS_DAYOFVISIT not NaN
date = date.join(data['SS_DAYOFVISIT'])
data.drop('SS_DAYOFVISIT', axis=1, inplace=True)
pres = pres.join(data['SS_ISTHEPATIENTR'].replace(yn_map))
data.drop('SS_ISTHEPATIENTR', axis=1, inplace=True)

data = replace_yn(data, ['SS_DOESTHEPATIEN'])
data = replace_checked(data, ['SS_SPINALCORDCOM', 'SS_BONEPAIN', 'SS_FATIGUE', 'SS_HYPERCALCEMIA', 'SS_RENALINSUFFIC', 'SS_ANEMIAHEMOGLO', 'SS_BONELESIONSLY', 'SS_BONELESIONSOS', 'SS_SOFTTISSUEPLA', 'SS_OFTTISSUEPLAS', 'SS_RECURRENTBACT', 'SS_SYMPTOMATICHY', 'SS_AMYLOIDOSIS', 'SS_OTHER'])
text = text.join(data['SS_SPECIFY'])
data.drop('SS_SPECIFY', axis=1, inplace=True)

# Go back and ensure:
#   NaN to unassessed patients
assessed = pd.notnull(date['SS_DAYOFVISIT'])
replace_cols = ['SS_SPINALCORDCOM', 'SS_BONEPAIN', 'SS_FATIGUE', 'SS_HYPERCALCEMIA', 'SS_RENALINSUFFIC', 'SS_ANEMIAHEMOGLO', 'SS_BONELESIONSLY', 'SS_BONELESIONSOS', 'SS_SOFTTISSUEPLA', 'SS_OFTTISSUEPLAS', 'SS_RECURRENTBACT', 'SS_SYMPTOMATICHY', 'SS_AMYLOIDOSIS', 'SS_OTHER']
for col in replace_cols:
    data.loc[assessed == 0, col] = np.nan

########################################################################################################################
# Plasmacytoma/soft tissue
# Combine these 2
data['ST_PLASMACYTOMAN'] = data[['ST_PLASMACYTOMAN', 'ST_PLASMACYTOMAN2']].mean(axis=1)
data.drop('ST_PLASMACYTOMAN2', axis=1, inplace=True)

st_map = {
    '': np.nan,
    'No Change': 0,
    'Yes': 1
}
data = replace_map(data, 'ST_INCREASEINNUM', st_map)
cats['ST_INCREASEINNUM'] = st_map

text = text.join(data['ST_RESULTOFEXAMI2'])
data.drop('ST_RESULTOFEXAMI2', axis=1, inplace=True)

pres = pres.join(data['ST_WASANASSESSME'].replace(yn_map))
data.drop('ST_WASANASSESSME', axis=1, inplace=True)
date = date.join(data['ST_ASSESSMENTDAT'])
data.drop('ST_ASSESSMENTDAT', axis=1, inplace=True)

data = replace_yn(data, ['ST_WEREANYSOFTTI'])
st_num_map = {
    '': np.nan,
    'Single': 1,
    'Multiple': 2
}
data = replace_map(data, 'ST_NUMBEROFPLASM', st_num_map)
cats['ST_NUMBEROFPLASM'] = st_num_map
data = replace_map(data, 'ST_WASTHEREACHAN', st_map)  # same map as above
cats['ST_WASTHEREACHAN'] = st_map
data = replace_map(data, 'ST_INCREASEINSIZE', st_map)  # same map as above
cats['ST_INCREASEINSIZE'] = st_map
st_result_map = {
    '': np.nan,
    'Reduction in size': -1,
    'No change': 0,
    'Increase in size': 1
}
data = replace_map(data, 'ST_RESULTOFEXAMI', st_result_map)
cats['ST_RESULTOFEXAMI'] = st_result_map

# These cols contain number or NaN as appropriate if the test was taken

########################################################################################################################
# Conventional cytogenetics/chromosomal damage tests
# Drop unknown col
data.drop('D_CM_enr', axis=1, inplace=True)

aneu_cat_map = {  # this is sort of ordinal, but unknown...
    '': np.nan,
    'Unknown': np.nan,
    'Not Done': np.nan,
    'None': 0,
    'Hypodiploid(44-45 Chromosomes)': -1,
    'Near diploid (44/45 to 46/47 Chromosomes)': 1,  # not "normal"
    'Hyperdiploid(>46/47 Chromosomes)': 2,
    'Hypotetraploid (near tetraploid) (>75 Chromosomes)': 3
}
data = replace_map(data, 'D_CM_ANEUPLOIDYCAT', aneu_cat_map)
cats['D_CM_ANEUPLOIDYCAT'] = aneu_cat_map

pres = pres.join(data['D_CM_WASCONVENTION'].replace(yn_map))
data.drop('D_CM_WASCONVENTION', axis=1, inplace=True)

date = date.join(data['D_CM_DAYPERFORMED'])
data.drop('D_CM_DAYPERFORMED', axis=1, inplace=True)

misc = misc.join(data['D_CM_NUMBEROFCELLS'])
data.drop('D_CM_NUMBEROFCELLS', axis=1, inplace=True)

text = text.join(data[['D_CM_KARYOTYPE', 'D_CM_SPECIFIY']])
data.drop(['D_CM_KARYOTYPE', 'D_CM_SPECIFIY'], axis=1, inplace=True)

cyto_map = {
    '': np.nan,
    'Test not perf': np.nan,
    'Indeterminant': np.nan,
    'Normal': 0,
    'Abnormal': 1
}
data = replace_map(data, 'D_CM_CYTOGENICSRES', cyto_map)
cats['D_CM_CYTOGENICSRES'] = cyto_map

# Keep D_CM_ABNORMALPERCE, D_CM_TOTALNUMBEROF, D_CM_abnres

data = replace_checked(data, ['D_CM_DEL13', 'D_CM_DEL17', 'D_CM_T614', 'D_CM_1PDELETION', 'D_CM_T814', 'D_CM_1QAMPLIFICATI', 'D_CM_T1114', 'D_CM_OTHER'])

pres = pres.join(data[['D_CM_cm', 'D_CM_abnmiss']])
data.drop(['D_CM_cm', 'D_CM_abnmiss'], axis=1, inplace=True)

# Drop redundant cols
data.drop(['D_CM_cres', 'D_CM_aneu', 'D_CM_abnres'], axis=1, inplace=True)

# Go back and ensure:
#   NaN to unassessed patients
assessed = pres['D_CM_WASCONVENTION'] == 1
replace_cols = ['D_CM_DEL13', 'D_CM_DEL17', 'D_CM_T614', 'D_CM_1PDELETION', 'D_CM_T814', 'D_CM_1QAMPLIFICATI', 'D_CM_T1114', 'D_CM_OTHER']
for col in replace_cols:
    data.loc[assessed == 0, col] = np.nan

########################################################################################################################
# D_IM (immunofluorescence?) stuff
# Drop redundant cols
data.drop(['D_IM_DCL_PATIENT_ID'], axis=1, inplace=True)

date = date.join(data[['D_IM_COLLECTION_DAY', 'D_IM_PROCESSING_DAY']])
data.drop(['D_IM_COLLECTION_DAY', 'D_IM_PROCESSING_DAY'], axis=1, inplace=True)

misc = misc.join(data[['D_IM_PATIENTS__MMRF_ELIGIBILITY']])
data.drop(['D_IM_PATIENTS__MMRF_ELIGIBILITY'], axis=1, inplace=True)

# Keep the other CD's as is
# Keep FGFR3

# Keep D_IM_FLOWCYT_PCT_ANEUPLOID_POPUL
# Keep D_IM_DNA_INDEX, D_IM_BRAF_STATUS

# Drop IM cols that are otherwise coded
#   It's possible we may use the readings, but unlikely
data.drop(['D_IM_IGH_SITE', 'D_IM_IGL_SITE', 'D_IM_CD38_DESCRIPTION', 'D_IM_CD138_DESCRIPTION', 'D_IM_CD45_DESCRIPTION', 'D_IM_CD56_DESCRIPTION', 'D_IM_LIGHT_CHAIN_BY_FLOW'], axis=1, inplace=True)

data = replace_map(data, 'D_IM_BRAF_CDNA_VARIANT', {'': np.nan})

braf_map = {
    '': np.nan,
    'V600E': 1
}
data = replace_map(data, 'D_IM_BRAF_PROTEIN_VARIANT', braf_map)
cats['D_IM_BRAF_PROTEIN_VARIANT'] = braf_map

# Keep coded descs

date = date.join(data[['D_IM_mrgday']])
data.drop(['D_IM_mrgday'], axis=1, inplace=True)

misc = misc.join(data[['D_IM_BONE_REASONFORPROC', 'D_IM_BONE_SPECIFY2']])
data.drop(['D_IM_BONE_REASONFORPROC', 'D_IM_BONE_SPECIFY2'], axis=1, inplace=True)

# Keep correlated disease progression?

data.drop(['D_IM_enr'], axis=1, inplace=True)

# Keep derived reason for sample?

date = date.join(data[['D_IM_crdiff', 'D_IM_pddiff']])
data.drop(['D_IM_crdiff', 'D_IM_pddiff'], axis=1, inplace=True)

# Flags for acceptable samples
pres = pres.join(data[['D_IM_det_plasma_cells', 'D_IM_detectable', 'D_IM_kaplam', 'D_IM_crhit', 'D_IM_pdhit', 'D_IM_hit']])
data.drop(['D_IM_det_plasma_cells', 'D_IM_detectable', 'D_IM_kaplam', 'D_IM_crhit', 'D_IM_pdhit', 'D_IM_hit'], axis=1, inplace=True)

########################################################################################################################
# Lab results - keep all numeric lab results

########################################################################################################################
# Survey questions
col_names = list(data)

q1_ind = col_names.index('D_QOL_Q1')
q50_ind = col_names.index('D_QOL_Q50')
survey = survey.join(data.ix[:, q1_ind:q50_ind+1])

q_sum_start_ind = col_names.index('D_QOL_mPF2')
q_sum_end_ind = col_names.index('D_QOL_MYSE')
survey_agg = survey_agg.join(data.ix[:, q_sum_start_ind:q_sum_end_ind+1])

pres = pres.join(data[['D_QOL_QLQ_MY20', 'D_QOL_QLQ_C30']].replace(yn_map))
date = date.join(data[['D_QOL_QOL_DY']])

q_start_ind = col_names.index('D_QOL_Q1')
q_end_ind = col_names.index('D_QOL_enr')
data.drop(data.columns[q_start_ind:q_end_ind+1], axis=1, inplace=True)

########################################################################################################################
# Keep cols with percent abonrmal cells
#   D_TRI_CF_T1420ABNORMAL, D_TRI_CF_1PDELETIONABN, D_TRI_CF_1PAMPLIFICATI2, D_TRI_CF_P53LOCUSCHROM2,
#   D_TRI_CF_OTHERSPECIFYA, D_TRI_CF_OTHERSPECIFYA2

# Keep cols with number abnormal cells
#   D_TRI_CF_1PDELETIONOFA, D_TRI_CF_1QAMPLIFICATI, D_TRI_CF_P53LOCUSCHROM3, D_TRI_CF_OTHEROFABNORM

# Drop cols that just say the result is unknown - these don't seem to correspond to anything
data.drop(['D_TRI_CF_T1420', 'D_TRI_CF_CCOUNT1PUN', 'D_TRI_CF_ICELLS1PUN', 'D_TRI_CF_1PDELETIONUN',
           'D_TRI_CF_CCOUNT1QUN', 'D_TRI_CF_ICELLS1QUN', 'D_TRI_CF_1PAMPUN', 'D_TRI_CF_CCOUNTP53UN',
           'D_TRI_CF_ICELLSP53UN', 'D_TRI_CF_P53LOCUN', 'D_TRI_CF_CCOUNTOTHERUN', 'D_TRI_CF_ICELLSOTHERUN',
           'D_TRI_CF_OTHERSPECUN', 'D_TRI_CF_CCOUNTOTHER2UN', 'D_TRI_CF_ICELLSOTHER2UN', 'D_TRI_CF_OTHERSPEC2UN',
           ], axis=1, inplace=True)

# Drop cols of counts of cells - not useful features
data.drop(['D_TRI_CF_1PDELETIONOFC', 'D_TRI_CF_1PAMPLIFICATI', 'D_TRI_CF_P53LOCUSCHROM', 'D_TRI_CF_OTHERSPECIFYO',
           'D_TRI_CF_OTHEROFABNORM2', 'D_TRI_CF_OTHERSPECIFYO2'], axis=1, inplace=True)

# Presence/absence of abnormalities
p1del_map = {
    '': np.nan,
    '.': np.nan,  # not sure about this
    'Not Done': np.nan,
    'No': 0,
    'Yes': 1,
}
data = replace_map(data, 'D_TRI_CF_ABNORMALITYPR12', p1del_map)
cats['D_TRI_CF_ABNORMALITYPR12'] = p1del_map

data = replace_map(data, 'D_TRI_CF_ABNORMALITYPR13', p1del_map)
cats['D_TRI_CF_ABNORMALITYPR13'] = p1del_map

data = replace_map(data, 'D_TRI_CF_ABNORMAILITYP', p1del_map)
cats['D_TRI_CF_ABNORMAILITYP'] = p1del_map

data = replace_map(data, 'D_TRI_CF_ABNORMALITYPR14', p1del_map)
cats['D_TRI_CF_ABNORMALITYPR14'] = p1del_map

data = replace_map(data, 'D_TRI_CF_ABNORMALITYPR15', p1del_map)
cats['D_TRI_CF_ABNORMALITYPR15'] = p1del_map

# Freeform text
text = text.join(data[['D_TRI_CF_SPECIFY2', 'D_TRI_CF_SPECIFY3']])
data.drop(['D_TRI_CF_SPECIFY2', 'D_TRI_CF_SPECIFY3'], axis=1, inplace=True)

# Trisomies
data = replace_checked(data, ['D_TRI_CF_TRISOMIES3', 'D_TRI_CF_TRISOMIES5', 'D_TRI_CF_TRISOMIES7', 'D_TRI_CF_TRISOMIES9', 'D_TRI_CF_TRISOMIES11', 'D_TRI_CF_TRISOMIES15', 'D_TRI_CF_TRISOMIES19', 'D_TRI_CF_TRISOMIES21', 'D_TRI_CF_TRISOMIESOTH'])
text = text.join(data[['D_TRI_CF_SPECIFY']])
data.drop(['D_TRI_CF_SPECIFY'], axis=1, inplace=True)

# Drop redundant no or unknown trisomies - they will be nans in the above cols
data.drop(['D_TRI_CF_TRISOMIESNONE', 'D_TRI_CF_TRISOMIESUNK', 'D_TRI_CF_TRISOMIESNTRPTD'], axis=1, inplace=True)

# Drop redundant trisomy col that reports presence
data.drop(['D_TRI_CF_TRISOMIES'], axis=1, inplace=True)

# Metadata
pres = pres.join(data[['D_TRI_CF_WASCYTOGENICS', 'D_TRI_CF_WASCLGFISHORP']].replace(yn_map))
data.drop(['D_TRI_CF_WASCYTOGENICS', 'D_TRI_CF_WASCLGFISHORP'], axis=1, inplace=True)

date = date.join(data['D_TRI_CF_DAYPERFORMED'])
data.drop('D_TRI_CF_DAYPERFORMED', axis=1, inplace=True)

# Another block of classes of abnormalities
# Keep cols with percent abonrmal cells
#   D_TRI_CF_DEL13ABNORMAL, D_TRI_CF_13QABNORMALCE, D_TRI_CF_DEL17ABNORMAL, D_TRI_CF_17PABNORMALCE,
#   D_TRI_CF_T414ABNORMALC, D_TRI_CF_T614ABNORMALC, D_TRI_CF_T814ABNORMALC, D_TRI_CF_T1114OFABNORM,
#   D_TRI_CF_T1114ABNORMAL, D_TRI_CF_T1214ABNORMAL, D_TRI_CF_T1416OFABNORM, D_TRI_CF_T1416ABNORMAL,
#

# Keep cols with number abnormal cells
#   D_TRI_CF_DEL13QOFABNOR, D_TRI_CF_DEL17OFABNORM, D_TRI_CF_DEL17POFABNOR, D_TRI_CF_T414OFABNORMA,
#   D_TRI_CF_T814OFABNORMA, D_TRI_CF_T1214OFABNORM, D_TRI_CF_T1420OFABNORM,

# Drop cols that just say the result is unknown - these don't seem to correspond to anything
data.drop(['D_TRI_CF_CCOUNTDEL13UN', 'D_TRI_CF_ICELLDEL13UN', 'D_TRI_CF_DEL13UN', 'D_TRI_CF_CCOUNTDEL13QUN',
           'D_TRI_CF_ICELLSDEL13UN', 'D_TRI_CF_DEL12QUN', 'D_TRI_CF_CCOUNTDEL17UN', 'D_TRI_CF_ICELLSDEL17UN',
           'D_TRI_CF_DEL17UN', 'D_TRI_CF_CCOUNTDEL17PUN', 'D_TRI_CF_ICELLSDEL17PUN', 'D_TRI_CF_DEL17P',
           'D_TRI_CF_CCOUNTT414UN', 'D_TRI_CF_ICELLST414UN', 'D_TRI_CF_T414UN', 'D_TRI_CF_CCOUNTT614UN',
           'D_TRI_CF_ICELLST614UN', 'D_TRI_CF_T614UN', 'D_TRI_CF_CCOUNTT814UN', 'D_TRI_CF_ICELLST814UN',
           'D_TRI_CF_T814UN', 'D_TRI_CF_CCOUNTT1114UN', 'D_TRI_CF_ICELLST1114UN', 'D_TRI_CF_T1114UN',
           'D_TRI_CF_CCOUNTT1214UN', 'D_TRI_CF_ICELLST1214UN', 'D_TRI_CF_T1214UN', 'D_TRI_CF_CCOUNTT1416UN',
           'D_TRI_CF_ICELLST1416', 'D_TRI_CF_T1416', 'D_TRI_CF_CCOUNTT1420UN', 'D_TRI_CF_ICELLST1420'
           ], axis=1, inplace=True)

# Drop cols of counts of cells - not useful features
data.drop(['D_TRI_CF_DEL13OFABNORM', 'D_TRI_CF_DEL13OFCELLSE', 'D_TRI_CF_13QOFCELLSEXA', 'D_TRI_CF_DEL17OFCELLSE',
           'D_TRI_CF_17POFCELLSEXA', 'D_TRI_CF_T414OFCELLSEX', 'D_TRI_CF_T614OFABNORMA', 'D_TRI_CF_T614OFCELLSEX',
           'D_TRI_CF_T814OFCELLSEX', 'D_TRI_CF_T1114OFCELLSE', 'D_TRI_CF_T1214OFCELLSE', 'D_TRI_CF_T1416OFCELLSE',
           'D_TRI_CF_T1420OFCELLSE',
           ], axis=1, inplace=True)

# Presence/absence of abnormalities
data = replace_map(data, 'D_TRI_CF_ABNORMALITYPR', p1del_map)
cats['D_TRI_CF_ABNORMALITYPR'] = p1del_map

data = replace_map(data, 'D_TRI_CF_ABNORMALITYPR10', p1del_map)
cats['D_TRI_CF_ABNORMALITYPR10'] = p1del_map

data = replace_map(data, 'D_TRI_CF_ABNORMALITYPR2', p1del_map)
cats['D_TRI_CF_ABNORMALITYPR2'] = p1del_map

data = replace_map(data, 'D_TRI_CF_ABNORMALITYPR11', p1del_map)
cats['D_TRI_CF_ABNORMALITYPR11'] = p1del_map

data = replace_map(data, 'D_TRI_CF_ABNORMALITYPR3', p1del_map)
cats['D_TRI_CF_ABNORMALITYPR3'] = p1del_map

data = replace_map(data, 'D_TRI_CF_ABNORMALITYPR4', p1del_map)
cats['D_TRI_CF_ABNORMALITYPR4'] = p1del_map

data = replace_map(data, 'D_TRI_CF_ABNORMALITYPR5', p1del_map)
cats['D_TRI_CF_ABNORMALITYPR5'] = p1del_map

data = replace_map(data, 'D_TRI_CF_ABNORMALITYPR6', p1del_map)
cats['D_TRI_CF_ABNORMALITYPR6'] = p1del_map

data = replace_map(data, 'D_TRI_CF_ABNORMALITYPR7', p1del_map)
cats['D_TRI_CF_ABNORMALITYPR7'] = p1del_map

data = replace_map(data, 'D_TRI_CF_ABNORMALITYPR8', p1del_map)
cats['D_TRI_CF_ABNORMALITYPR8'] = p1del_map

data = replace_map(data, 'D_TRI_CF_ABNORMALITYPR9', p1del_map)
cats['D_TRI_CF_ABNORMALITYPR9'] = p1del_map

# More metadata - number versions of above
pres = pres.join(data[['D_TRI_cf', 'D_TRI_clg']])
data.drop(['D_TRI_cf', 'D_TRI_clg'], axis=1, inplace=True)

# Keep cols for present/absent
# Note/TODO: These semi-redundant with above - keep the ones at the level of granularity we want
#   D_TRI_abn1, D_TRI_abn3, D_TRI_abn4, D_TRI_abn5, D_TRI_abn6, D_TRI_abn7, D_TRI_abn8, D_TRI_abn9, D_TRI_abn10,
#   D_TRI_abn12, D_TRI_abn13, D_TRI_abnoth
# Drop redundant cols
data.drop(['D_TRI_missabn', 'D_TRI_noneabn'], axis=1, inplace=True)
# Drop more redundant trisomies data
data.drop(['D_TRI_trisomies', 'D_TRI_trisnd', 'D_TRI_trisnone'], axis=1, inplace=True)

# Go back and ensure:
#   NaN to unassessed patients
assessed = pres['D_TRI_CF_WASCYTOGENICS'] == 1
replace_cols = ['D_TRI_CF_TRISOMIES3', 'D_TRI_CF_TRISOMIES5', 'D_TRI_CF_TRISOMIES7', 'D_TRI_CF_TRISOMIES9', 'D_TRI_CF_TRISOMIES11', 'D_TRI_CF_TRISOMIES15', 'D_TRI_CF_TRISOMIES19', 'D_TRI_CF_TRISOMIES21', 'D_TRI_CF_TRISOMIESOTH']
for col in replace_cols:
    data.loc[assessed == 0, col] = np.nan

########################################################################################################################
# Another hyperdiploid measure
data = replace_yn(data, ['Hyperdiploid'])

# Drop missing CM versions of the translocation tests
data.drop(['CM_ABNORMALITYPR9', 'CM_ABNORMALITYPR3', 'CM_ABNORMALITYPR8'], axis=1, inplace=True)

# Drop study name data
data.drop(['STUDY_ID'], axis=1, inplace=True)

# Grab Palumbo cytogenetic tests
data = replace_map(data, 'Deletion_17q13', p1del_map)
cats['Deletion_17q13'] = p1del_map

# Drop a col with only 2 entries
data.drop(['Other_cytogenetic_analysis'], axis=1, inplace=True)

# AT data?
data = replace_map(data, 'AT_ABSENCEOFCLON', p1del_map)
cats['AT_ABSENCEOFCLON'] = p1del_map

# Keep 2nd col of disease progression
#   TODO: compare to other col, merge?
endp = endp.join(data['AT_CDRREVIEWPERI'].replace(response_map))
data.drop('AT_CDRREVIEWPERI', axis=1, inplace=True)
cats['AT_CDRREVIEWPERI'] = response_map

# Semi-structured text
text = text.join(data[['AT_CDRCOMMENTPER']])
data.drop(['AT_CDRCOMMENTPER'], axis=1, inplace=True)

# Palumbo bone data
# Drop for now - there's only 24 measurements, their code is indecipherable, don't know how to merge
data.drop(['Bone_lytic_evaluation', 'MYELOMA_INVOLVEMENT_IN_THE_BONE_', 'diagn_Plasmocitoma_number', 'diagn_Plasmocitoma_size', 'diagn_Plasmocitoma_Site', 'diagn_test_plasmocitoma'], axis=1, inplace=True)

# Treat CMMC effect as endpoint
cmmc_map = {
    'Relapse/Progression': 10,
    'Remission/Response': 9,
    'complete response': 8,
    'Restaging': 7,
    'Post-transplant': 6,
    'Pre-transplant': 5,
    'follow up': 4,
    'Confirmation': 3,
    'Baseline': 2,
    'Screening': 1,
    'Other': 0,
    '': np.nan
}
endp = endp.join(data['CMMC_VISIT_NAME'].replace(cmmc_map))
data.drop('CMMC_VISIT_NAME', axis=1, inplace=True)
cats['CMMC_VISIT_NAME'] = cmmc_map

# Text indicates "extrapolated result" and similar
text = text.join(data['CMMC_COMMENTS'])
data.drop('CMMC_COMMENTS', axis=1, inplace=True)

# Keep CMMC causes
cmmc2_map = {
    'Confirm Progression': 10,
    'Confirm Response/sC': 8,
    'Restaging/disease s': 7,
    'Post transplant': 6,
    'Pre transplant': 5,
    'Baseline': 2,
    'Other': 0,
    '': np.nan
}
data = replace_map(data, 'CMMC_REASONFORPROC', cmmc2_map)
cats['CMMC_REASONFORPROC'] = cmmc2_map

cml_map = {
    'Confirm myeloid leukemia': 5,
    'S/P BMT': 4,
    'per physician': 3,
    'response to date PR...confirm or review for pd': 2,
    'Follow Up': 1,
    '': np.nan
}
data = replace_map(data, 'CMMC_BONE_SPECIFY2', cml_map)
cats['CMMC_BONE_SPECIFY2'] = cml_map

# Keep CMMC_REASONCODE, CMMC_CRHIT, CMMC_PDHIT, CMMC_RECDY

# Drop enr
data.drop('enr', axis=1, inplace=True)

# Save processed tables
output_dir = 'data/processed'

data.set_index(['PUBLIC_ID', 'VISIT'], inplace=True)
data.to_csv(os.path.join(output_dir, 'clinical_data.csv'))

endp.set_index(['PUBLIC_ID', 'VISIT'], inplace=True)
endp.to_csv(os.path.join(output_dir, 'clinical_endp.csv'))

# pres.to_csv(os.path.join(output_dir, 'clinical_pres.csv'))
# text.to_csv(os.path.join(output_dir, 'clinical_text.csv'))
# misc.to_csv(os.path.join(output_dir, 'clinical_misc.csv'))
# date.to_csv(os.path.join(output_dir, 'clinical_date.csv'))
# treat.to_csv(os.path.join(output_dir, 'clinical_treat.csv'))
# survey.to_csv(os.path.join(output_dir, 'clinical_survey.csv'))
# survey_agg.to_csv(os.path.join(output_dir, 'clinical_survey_agg.csv'))

