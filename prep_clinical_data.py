# Main function that goes thru PER_PATIENT_VISIT clinical data table and turns it into a usable form
# Goes thru columns manually and assess quality
# A unique index is made from PUBLIC_ID + VISIT

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

checked_map = {
    'Checked': 1,
    '': 0  # Note: This may be changed to nan, if that makes more sense
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

########################################################################################################################
# Keep PUBLIC_ID as main index

# Drop spectrum id
data.drop('SPECTRUM_SEQ', axis=1, inplace=True)

# Keep visit day but drop coarse-grained visit interval
data.drop('VJ_INTERVAL', axis=1, inplace=True)

# Drop informed consent day
data.drop('VISITDY', axis=1, inplace=True)

# Bone lesions check: Either 'checked' -> 1 or blank -> NaN
pres = pres.join(data['AT_DEFINITEDEVEL'].replace(yn_map))
data.drop('AT_DEFINITEDEVEL', axis=1, inplace=True)

# Treatment response checked at last visit
pres = pres.join(data['AT_WASANASSESSME'].replace(yn_map))
data.drop('AT_WASANASSESSME', axis=1, inplace=True)

# Response assessment date
date = date.join(data['AT_RESPONSEASSES'])
data.drop('AT_RESPONSEASSES', axis=1, inplace=True)

# Convert treatment response to numeric index
response_map = {
    'Stringent Complete Response (sCR)': 6,
    'Complete Response': 5,
    'Very Good Partial Response (VGPR)': 4,
    'Partial Response': 3,
    'Stable Disease': 2,
    'Progressive Disease': 1,
    '': np.nan
}
data = replace_map(data, 'AT_TREATMENTRESP', response_map)
cats['AT_TREATMENTRESP'] = response_map

# Basic checks on worsening conditions
data = replace_checked(data, ['AT_INCREASEOF25F', 'AT_SERUMMCOMPONE', 'AT_URINEMCOMPONE', 'AT_ONLYINPATIENT', 'AT_ONLYINPATIENT2', 'AT_DEVELOPMENTOF'])

# Bone marrow aspriate useful?
pres = add_checked_cols(pres, ['BMA_WASTHESCREENI'])
data.drop('BMA_WASTHESCREENI', axis=1, inplace=True)

# Bone marrow assessment done?
pres = add_checked_cols(pres, ['BA_WASABONEASSES'])
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
    '': np.nan,
    'Clinically significant abnormality related to MM': 1,
    'Clinically significant abnormality unrelated to MM': 2,
    'Normal or clinically insignificant abnormality': 3
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
misc = misc.join(data[['BB_IFNOPLEASESPE', 'BONE_SPECIFY', 'BONE_IFNOPLEASESPE', 'BONE_IFNOPLEASESPE2']])
data.drop(['BB_IFNOPLEASESPE', 'BONE_SPECIFY', 'BONE_IFNOPLEASESPE', 'BONE_IFNOPLEASESPE2'], axis=1, inplace=True)
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

# Drop informed consent cols
data.drop(['IC_DAYOFINFORME', 'IC_DAYOFVISIT', 'IC_DAYUNIVERSAL', 'IC_HASTHEPATIENT', 'IC_IFYESPLEASEPR', 'IC_IPERMITMYSAMP', 'IC_IPERMITMYSAMP2', 'IC_IPERMITRESEAR', 'IC_PATIENTHASREA', 'IC_WASTHISPATIEN2'], axis=1, inplace=True)
# Drop a couple more cols that are probably not important...
data.drop(['IC_DIDTHEPATIENT', 'IC_OTHREASON'], axis=1, inplace=True)
# All the patients are over 18 yo
data.drop(['IC_PATIENTISATLE'], axis=1, inplace=True)

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

# Other surgery
treat = treat.join(data['MMSURG_SPECIFY'])
data.drop('MMSURG_SPECIFY', axis=1, inplace=True)
treat = treat.join(data[['MMSURG_NONE', 'MMSURG_NONE', 'MMSURG_VERTEBROPLAST', 'MMSURG_KYPHOPLASTY', 'MMSURG_OTHER']].replace(yn_map))
data.drop(['MMSURG_NONE', 'MMSURG_NONE', 'MMSURG_VERTEBROPLAST', 'MMSURG_KYPHOPLASTY', 'MMSURG_OTHER'], axis=1, inplace=True)

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

# Supplemental treatments
supp_treatments = ['SUPP_ANTICOAGULATI', 'SUPP_NONE', 'SUPP_TRANSFUSION', 'SUPP_PACKEDREDBLOO', 'SUPP_PLATELETS', 'SUPP_OTHER', 'SUPP_RADIOTHERAPY', 'SUPP_DIALYSIS', 'SUPP_BISPHOSPHONAT', 'SUPP_WBCGROWTHFACT', 'SUPP_ERYTHROPOIESI', 'SUPP_ANTIEMETIC', 'SUPP_ANALGESICSFOR', 'SUPP_OPIOID', 'SUPP_OTHER2', 'SUPP_ANTIVIRAL', 'SUPP_MEDICATIONSFO', 'SUPP_MEDICATIONSFO2']
treat = treat.join(data[supp_treatments].replace(yn_map))
data.drop(supp_treatments, axis=1, inplace=True)
treat = treat.join(data[['SUPP_SPECIFY', 'SUPP_SPECIFYSITE']])
data.drop(['SUPP_SPECIFY', 'SUPP_SPECIFYSITE'], axis=1, inplace=True)

# Myeloma-specific related symptoms
date = date.join(data['SS_DAYOFVISIT'])
data.drop('SS_DAYOFVISIT', axis=1, inplace=True)
data = replace_yn(data, ['SS_ISTHEPATIENTR', 'SS_DOESTHEPATIEN'])
data = replace_checked(data, ['SS_SPINALCORDCOM', 'SS_BONEPAIN', 'SS_FATIGUE', 'SS_HYPERCALCEMIA', 'SS_RENALINSUFFIC', 'SS_ANEMIAHEMOGLO', 'SS_BONELESIONSLY', 'SS_BONELESIONSOS', 'SS_SOFTTISSUEPLA', 'SS_OFTTISSUEPLAS', 'SS_RECURRENTBACT', 'SS_SYMPTOMATICHY', 'SS_AMYLOIDOSIS', 'SS_OTHER'])
text = text.join(data['SS_SPECIFY'])
data.drop('SS_SPECIFY', axis=1, inplace=True)

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

# Cytogenetics/chromosomal damage - this is the same as in the separate cytogenetic_data folder

print(data['D_CM_enr'].unique())
print(data['D_CM_ANEUPLOIDYCAT'].unique())
print(data['D_CM_WASCONVENTION'].unique())

# TODO: Left off on col 227

# print(pres)
# print(misc)
# print(text)
# print(data)
# data.to_csv('test.csv')
