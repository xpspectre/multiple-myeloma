# Look at staging, survival, progression
# TODO: Right now, overall survival is assessed. Do MM-caused death survival as well.

import os
import pandas as pd
from lifelines import KaplanMeierFitter
import matplotlib.pyplot as plt

data_dir = 'data/processed'

# Load study data - to keep just the CoMMpass patients
study_data = pd.read_csv(os.path.join(data_dir, 'patient_study.csv'))
study_data.set_index('PUBLIC_ID', inplace=True)

endp_data = pd.read_csv(os.path.join(data_dir, 'patient_endp.csv'))
endp_data.set_index('PUBLIC_ID', inplace=True)

endp_data = endp_data[study_data['STUDY_ID'] == 1]

# Survival analysis
died = pd.notnull(endp_data['D_PT_deathdy']).astype(int)
last_observed = endp_data[['D_PT_lstalive', 'D_PT_deathdy']].max(axis=1)  # death date should be after last alive date, probably

# Overall survival
kmf = KaplanMeierFitter()
kmf.fit(last_observed, event_observed=died, label='Overall')

kmf.plot()
plt.xlabel('Time (days)')
plt.ylabel('Survival')
plt.title('CoMMpass Overall Survival')
plt.draw()

# Survival by stage
iss_stage = endp_data['D_PT_iss']
died_1 = died[iss_stage == 1]
died_2 = died[iss_stage == 2]
died_3 = died[iss_stage == 3]
last_observed_1 = last_observed[iss_stage == 1]
last_observed_2 = last_observed[iss_stage == 2]
last_observed_3 = last_observed[iss_stage == 3]

kmf = KaplanMeierFitter()
kmf.fit(last_observed_1, event_observed=died_1, label='1')
ax = kmf.plot()
kmf.fit(last_observed_2, event_observed=died_2, label='2')
kmf.plot(ax=ax)
kmf.fit(last_observed_3, event_observed=died_3, label='3')
kmf.plot(ax=ax)
plt.xlabel('Time (days)')
plt.ylabel('Survival')
plt.title('CoMMpass Survival by ISS Stage')
plt.draw()

# Survival by treatment
#   Regardless of stage (confounding: stage may determine treatment aggressiveness - massive chance of Simpson's paradox)
#   Regardless of combo treatments
#   Regardless of everything else...
# Since the treatments overlap, just the 4 main treatments isn't fair...
# Load treatment data
treat_data = pd.read_csv(os.path.join(data_dir, 'patient_treat.csv'))
treat_data.set_index('PUBLIC_ID', inplace=True)

# Go thru the formalism of consistently joining the tables and separating out the components as needed
data = endp_data
data = data.join(treat_data)

died = pd.notnull(data['D_PT_deathdy']).astype(int)
last_observed = data[['D_PT_lstalive', 'D_PT_deathdy']].max(axis=1)
cols = ['TREAT_BOR', 'TREAT_CAR', 'TREAT_IMI', 'line1sct']
treatments = ['Bortezomib', 'Carfilzomib', 'IMIDs', 'SCT']

kmf = KaplanMeierFitter()
for i, col in enumerate(cols):
    died_i = died[data[col] == 1]
    last_observed_i = last_observed[data[col] == 1]
    kmf.fit(last_observed_i, event_observed=died_i, label=treatments[i])
    if i == 0:
        ax = kmf.plot()
    else:
        kmf.plot(ax=ax)
    plt.xlabel('Time (days)')
    plt.ylabel('Survival')
    plt.title('CoMMpass Survival by Treatment\n(Treatments may overlap, may be dependent on severity)')
    plt.draw()

# 1st line stem cell transfer vs no
died_sct = died[data['line1sct'] == 1]
died_no_sct = died[data['line1sct'] == 0]
last_observed_sct = last_observed[data['line1sct'] == 1]
last_observed_no_sct = last_observed[data['line1sct'] == 0]

kmf = KaplanMeierFitter()
kmf.fit(last_observed_sct, event_observed=died_sct, label='SCT')
ax = kmf.plot()
kmf.fit(last_observed_no_sct, event_observed=died_no_sct, label='No SCT')
kmf.plot(ax=ax)
plt.xlabel('Time (days)')
plt.ylabel('Survival')
plt.title('CoMMpass Survival by 1st Line Stem Cell Transfer Treatment')
plt.draw()

plt.show()
