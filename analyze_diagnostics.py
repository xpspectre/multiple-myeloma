# Look at measures of patient health

from load_patient_data import load_per_visit_data
import numpy as np
import os
import cmsgpack

if __name__ == '__main__':
    per_visit_data, per_visit_dict, per_visit_fields = load_per_visit_data()

    data_dir = 'data/processed'
    diags_field_counts_file = 'patient_diagnostics_field_counts'

    # Get useful diagnostics
    # TODO: Filter to get list of only the diagnostics

    # Get diagnostic results that are:
    #   1) numeric
    #   2) >= x% of patients have measurements at some point
    frac_patient_present = 0.8

    patients = np.unique(per_visit_data['PUBLIC_ID'])
    n_patients = len(patients)

    field_counts = {}
    for field, type in per_visit_fields.items():
        if np.issubdtype(type, np.number):
            field_counts[field] = 0

    print('Found {n_num_fields}/{n_fields} numeric fields'.format(n_num_fields=len(field_counts), n_fields=len(per_visit_fields)))

    for patient in patients:
        patient_data = per_visit_data[per_visit_data['PUBLIC_ID'] == patient]
        for field in field_counts:
            col = np.asarray(patient_data[field])
            if not np.all(np.isnan(col)):
                field_counts[field] += 1

    with open(os.path.join(data_dir, diags_field_counts_file), 'wb') as f:
        f.write(cmsgpack.packb(field_counts))

    # with open(os.path.join(data_dir, diags_field_counts_file), 'rb') as f:
    #     field_counts = cmsgpack.unpackb(f.read())

    # Filter for fields that satisfy presence cutoff
    keep_fields = {}
    for field, count in field_counts.items():
        frac = count / n_patients
        if frac > frac_patient_present:
            keep_fields[field] = frac
    print('{n_keep}/{n} fields kept with > {frac:.2f} frac of patients having some value'.format(n_keep=len(keep_fields), n=len(field_counts), frac=frac_patient_present))

    # Grab datapoints for kept fields from each patient
    per_visit_data_filtered = {}
    for patient in patients:
        patient_data_filtered = {}
        patient_data = per_visit_data[per_visit_data['PUBLIC_ID'] == patient]
        dates = np.asarray(patient_data['VISITDY'])
        for field in keep_fields:
            col = np.asarray(patient_data[field])
            valid = np.isnan(col) == False
            x = col[valid]
            date = dates[valid]
            patient_data_filtered[field] = (date, x)
        per_visit_data_filtered[patient] = patient_data_filtered

    # TODO: Separate out diagnostics that were taken at negative days
    # TODO: Turn into feature vectors
    diags_filtered_file = 'patient_diagnostics_filtered'
    with open(os.path.join(data_dir, diags_filtered_file), 'wb') as f:
        f.write(cmsgpack.packb(per_visit_data_filtered))

    # with open(os.path.join(data_dir, diags_filtered_file), 'rb') as f:
    #     per_visit_data_filtered = cmsgpack.unpackb(f.read())
