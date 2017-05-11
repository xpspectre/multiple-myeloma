# Analyze treatment regimens by patient
# Note: The treatment_resp dataset seems to be more complete, and the treatment_regimen dataset seems to have info
#   on dosage (which we're tentatively ignoring for now)

import os
from load_patient_data import load_treatment_regimen, load_treatment_resp
import numpy as np
import scipy.io as sio
import cmsgpack


def merge_intervals(intervals):
    """Merge list of possibly overlapping (a, b) intervals. Requires a < b already for each pair.
    https://codereview.stackexchange.com/questions/69242/merging-overlapping-intervals
    """
    sorted_by_lower_bound = sorted(intervals, key=lambda tup: tup[0])
    merged = []
    for higher in sorted_by_lower_bound:
        if not merged:
            merged.append(higher)
        else:
            lower = merged[-1]
            # test for intersection between lower and higher:
            # we know via sorting that lower[0] <= higher[0]
            if higher[0] <= lower[1]:
                upper_bound = max(lower[1], higher[1])
                merged[-1] = (lower[0], upper_bound)  # replace by merged interval
            else:
                merged.append(higher)
    return merged


def build_treatment_regimens():
    """Use STAND_ALONE_TREATMENT_REGIMEN data"""
    treatment_data, treatment_dict, treatment_fields = load_treatment_regimen()

    # First build list of all the therapies
    therapies = np.unique(treatment_data['MMTX_THERAPY'])
    print('Therapies: {therapies}'.format(therapies=therapies))

    patients = np.unique(treatment_data['PUBLIC_ID'])

    # Create a condensed listing of patient therapies:
    #   dict of key = patient ID, val = dict of key = therapy name, val = list of tuples of (start,end) times
    # What do negative values for start and stop of therapy mean?
    patient_therapies = {}
    for patient in patients:
        patient_data = treatment_data[treatment_data['PUBLIC_ID'] == patient]

        patient_therapy = {}
        for i, row in patient_data.iterrows():
            therapy = row['MMTX_THERAPY']
            start = row['startday']
            stop = row['stopday']  # this can be NaN (blank in the original dataset) indicating not known

            if start > stop:  # just in case?
                start, stop = stop, start

            # Skip invalid (blank therapy)
            if len(therapy) == 0:
                continue

            if therapy not in patient_therapy:
                patient_therapy[therapy] = [(start, stop)]
            else:
                patient_therapy[therapy].append((start,stop))

        patient_therapies[patient] = patient_therapy

    return patient_therapies, therapies


def build_treatment_resps():
    treatment_resp_data, treatment_resp_dict, treatment_resp_fields = load_treatment_resp()

    # Collect all the treatments, separating out multiple treatments into the individual therapies (drugs)
    mixed_therapies = treatment_resp_data['trtname']
    therapies = set()
    for i, entry in mixed_therapies.iteritems():
        entries = entry.split('/')
        for therapy in entries:
            therapies.add(therapy)
    print('{n} different therapies found'.format(n=len(therapies)))

    # Get each patient's treatment intervals with each therapy
    patients = np.unique(treatment_resp_data['PUBLIC_ID'])
    patient_therapies = {}
    for patient in patients:
        patient_data = treatment_resp_data[treatment_resp_data['PUBLIC_ID'] == patient]

        patient_therapy = {}
        for i, row in patient_data.iterrows():
            therapy = row['trtname']
            start = row['trtstdy']
            stop = row['trtendy']

            entries = therapy.split('/')
            for entry in entries:
                # Skip invalid (blank therapy)
                if len(therapy) == 0:
                    continue

                if entry not in patient_therapy:
                    patient_therapy[entry] = [(start, stop)]
                else:
                    patient_therapy[entry].append((start, stop))

        # Clean up overlapping therapies
        merged_patient_therapies = {}
        for therapy, times in patient_therapy.items():
            merged_patient_therapies[therapy] = merge_intervals(times)

        patient_therapies[patient] = merged_patient_therapies

    return patient_therapies, therapies


def cleanup_treatments(patient_timeseries, therapies):
    """Convert somewhat messy list of therapies to better-defined therapies and replace in patient profiles.
    Some therapies have different names (pick a definitive one)
    Combine uncommon therapies into an 'Other' category"""

    # First get the number of patients on each therapy and look at the list to decide below
    t = dict.fromkeys(therapies, 0)
    for patient, tis in patient_timeseries.items():
        for ti in tis:
            t[ti] += 1

    # Look at the main ones and decide by eye
    #   All other therapies are put into a catch-all called 'Other'
    main_therapies = ['Bortezomib', 'Carfilzomib', 'Cyclophosphamide', 'Dexamethasone', 'Lenalidomide' ,'Melphalan']

    # Map misspellings of main therapies to above
    map_misspell = {
        'Melfalan': 'Melphalan'
    }

    # Map all thalidomide derivatives to Lenalidamide - there's a lot less of them but they're probably significant for the same reason
    #   And just OR the treatment statuses
    map_lidomides = {
        'Pomalidomide': 'Lenalidomide',
        'Thalidomide': 'Lenalidomide'
    }

    # Remap all patients
    n_times = patient_timeseries['MMRF_1011']['Bortezomib'].size  # sample that should be present
    patient_therapies_ = {}
    for patient, tis in patient_timeseries.items():
        tis_ = {}
        val_lidomides = np.zeros((n_times,), dtype=np.int8)
        val_others = np.zeros((n_times,), dtype=np.int8)

        for ti, vals in tis.items():
            if ti in map_misspell:
                tis_[map_misspell[ti]] = vals
            elif ti in map_lidomides or ti == 'Lenalidomide':
                val_lidomides = np.logical_or(val_lidomides, vals)
            elif ti in main_therapies:
                tis_[ti] = vals
            else:
                val_others = np.logical_or(val_others, vals)

        tis_['Lenalidomide'] = val_lidomides.astype(np.int8)
        tis_['Other'] = val_others.astype(np.int8)

        patient_therapies_[patient] = tis_

    main_therapies.append('Other')
    return patient_therapies_, main_therapies


def build_treatment_features(patient_therapies, time_interval=90, time_final=3650):
    """Convert treatment data in the form of intervals to timeseries features. Inputs are in days."""
    # Create timeseries for each patient
    #   If the patient was on a therapy at all in the interval [0, n), mark as 1, else 0. Interval granularity in days
    # This is inefficient, but ultimately only takes a few sec and only has to be run once
    intervals = np.arange(0, time_final, time_interval, dtype=np.float32)
    n_intervals = intervals.size - 1
    patient_timeseries = {}

    for patient, therapies in patient_therapies.items():
        patient_timeserie = {}
        for therapy, times in therapies.items():
            x = np.zeros((n_intervals,), dtype=np.int8)
            for i in range(1, n_intervals):
                i_start = intervals[i-1]
                i_stop = intervals[i]
                for time in times:
                    start = time[0]
                    stop = time[1]
                    if (start <= i_start < stop) or (start <= i_stop < stop) or (i_start <= start and stop < i_stop):
                        x[i-1] = 1
                        break
            patient_timeserie[therapy] = x
        patient_timeseries[patient] = patient_timeserie

    # Toss last timepoint, which represents an interval starting at the end time
    intervals = intervals[:-1]

    return patient_timeseries, intervals

if __name__ == '__main__':
    patient_therapies, therapies = build_treatment_resps()
    patient_timeseries, intervals = build_treatment_features(patient_therapies)
    patient_timeseries, therapies = cleanup_treatments(patient_timeseries, therapies)

    # Build patients x therapy x time matrix
    n_patients = len(patient_timeseries)
    n_therapies = len(therapies)
    n_times = len(intervals)
    patients = sorted(list(patient_timeseries.keys()))
    data = np.zeros((n_patients, n_therapies, n_times), dtype=np.int8)
    for i in range(n_patients):
        patient = patients[i]
        for j in range(n_therapies):
            therapy = therapies[j]
            if therapy in patient_timeseries[patient]:
                data[i, j, :] = patient_timeseries[patient][therapy]

    data_dir = 'data/processed'
    treatment_features_file = 'patient_treatment_timeseries'
    with open(os.path.join(data_dir, treatment_features_file), 'wb') as f:
        f.write(cmsgpack.packb((patient_timeseries, intervals)))

    treatment_matrix_file = os.path.join(data_dir, 'patient_treatment_matrix.mat')
    sio.savemat(treatment_matrix_file, {'data': data, 'patients': patients, 'therapies': therapies, 'intervals': intervals})
