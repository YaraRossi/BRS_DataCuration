import matplotlib.pyplot as plt
import pandas as pd
import os
import math

def angular_diff_deg(est, true):
    """
    Return signed minimal difference est - true in degrees, range [-180, 180).
    est and true can be scalars, numpy arrays, or pandas Series.
    """
    diff = est - true
    # wrap to [-180, 180)
    diff_wrapped = ((diff + 180) % 360) - 180
    return diff_wrapped

root = "/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data/correlation_baz_csv"

data_all_events = {}

for csv in os.listdir(root):
    if csv.endswith(".csv"):
        if '0.2' in csv: # only choose one frequency type
            df = pd.read_csv(os.path.join(root, csv))

            # Extract key (before first underscore)
            folder = csv.split('_')[0]
            freqmin = csv.split('_')[1]
            freqmax = csv.split('_')[2]

            # Add to dictionary
            data_all_events[folder] = df

# plot data

fig, axs = plt.subplots(1,1,figsize=(5,5))
fig.suptitle(f'Cross-correlation value for al events | {freqmin} - {freqmax} Hz')
for event in data_all_events:
    event_data = data_all_events[event]

    rel_time = event_data['time']/max(event_data['time'])
    axs.plot(rel_time, event_data['corr_value_midpoint'], label=event)

axs.set_xlabel('Relative time')
axs.set_ylabel('CC')
axs.legend(loc='lower right')

plt.show()

