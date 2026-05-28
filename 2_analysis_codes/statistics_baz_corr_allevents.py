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
root_savefig = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/2_Figures/Backazimuth_estimations/statistics'
data_all_events = {}

#location = 'LIV_internal' # LIV_external, LIV_internal, HAN_external
for location in ['LIV_external', 'LIV_internal', 'HAN_external']:
    for method in ['mid','max']:
        data_all_events = {}
        for csv in os.listdir(root):
            if csv.endswith(".csv"):
                if '0.01_0.1' in csv:  # only choose one frequency type
                    if location in csv:
                        df = pd.read_csv(os.path.join(root, csv))

                        # Extract key (before first underscore)
                        folder = csv.split('_')[0]
                        freqmin = csv.split('_')[-6]
                        freqmax = csv.split('_')[-5]

                        # Add to dictionary
                        data_all_events[folder] = df

        # plot data
        number_events = len(data_all_events)+1 # +1 to have one with all events
        n = math.ceil(math.sqrt(number_events))   # rows
        m = math.ceil(number_events / n)
        nbins=20# columns
        nbins_hc=10
        mincorr = 0.8
        mincorr2 = 0.9

        all_baz_midpoint = []   # collect all midpoint errors over all events
        all_baz_midpoint2 = []   # collect all midpoint errors over all events
        all_baz_max = []        # if you want also max errors
        all_baz_max2 = []        # if you want also max errors

        fig, axs = plt.subplots(n,m,figsize=(11,11), sharex=True)
        fig.suptitle(f'{location} BAz error - {method}points for CC>{mincorr} | {freqmin} - {freqmax} Hz')
        for i, event in zip(range(number_events-1), data_all_events):
            ax = axs.flat[i]
            event_data = data_all_events[event]
            #baz_midpoint = angular_diff_deg(event_data['baz_value_midpoint'], event_data['true_baz'])
            #baz_max = angular_diff_deg(event_data['baz_value_max'], event_data['true_baz'])
            #ax.hist(baz_midpoint,nbins, fill=False, ec='lightblue')
            #ax.hist(baz_max,nbins, fill=False, ec='pink')

            # only where correlation is high:
            event_data_hc = event_data[event_data['corr_value_max']>mincorr]
            baz_midpoint = angular_diff_deg(event_data_hc['baz_value_midpoint'], event_data_hc['true_baz'])
            event_data_hc = event_data[event_data['corr_value_midpoint']>mincorr]
            baz_max = angular_diff_deg(event_data_hc['baz_value_max'], event_data_hc['true_baz'])
            if method == 'mid':
                ax.hist(baz_midpoint,nbins_hc, fill=True, ec='lightblue', alpha=0.5, label=f'CC>{mincorr}')
            elif method == 'max':
                ax.hist(baz_max,nbins_hc, fill=True, ec='lightblue', alpha=0.5, label=f'CC>{mincorr}')

            #ax.hist(baz_max,nbins_hc, fill=True, ec='pink',alpha=0.5)

            # only where correlation is high:
            event_data_hc = event_data[event_data['corr_value_max']>mincorr2]
            baz_midpoint2 = angular_diff_deg(event_data_hc['baz_value_midpoint'], event_data_hc['true_baz'])
            event_data_hc = event_data[event_data['corr_value_midpoint']>mincorr2]
            baz_max2 = angular_diff_deg(event_data_hc['baz_value_max'], event_data_hc['true_baz'])
            if method == 'mid':
                ax.hist(baz_midpoint2,nbins_hc, fill=True, ec='darkblue', alpha=0.5, label=f'CC>{mincorr2}')
            elif method == 'max':
                ax.hist(baz_max2,nbins_hc, fill=True, ec='darkblue', alpha=0.5, label=f'CC>{mincorr2}')


            ax.text(0.05,0.9,event, transform=ax.transAxes)
            ax.legend(loc='center left')

            all_baz_midpoint.extend(baz_midpoint.values)
            all_baz_midpoint2.extend(baz_midpoint2.values)
            all_baz_max.extend(baz_max.values)
            all_baz_max2.extend(baz_max2.values)

        if method == 'mid':
            axs.flat[-1].hist(all_baz_midpoint, nbins, fill=True, ec='k', label=f'CC>{mincorr}') # all event in one histogram
            axs.flat[-1].hist(all_baz_midpoint2, nbins, fill=True, ec='k', label=f'CC>{mincorr2}') # all event in one histogram
        elif method == 'max':
            axs.flat[-1].hist(all_baz_max, nbins, fill=True, ec='k', label=f'CC>{mincorr}') # all event in one histogram
            axs.flat[-1].hist(all_baz_max2, nbins, fill=True, ec='k', label=f'CC>{mincorr2}') # all event in one histogram
        axs.flat[-1].text(0.02,0.9,f'all events', transform=axs.flat[-1].transAxes)
        axs.flat[-1].legend(loc='center left')

        for i in range(m):
            axs[-1,i].set_xlabel('BAz error [°]')
            axs[i,0].set_ylabel(f'# of windows')
        axs[0,0].set_xlim(left=-180,right=180)

        plt.savefig(f'{root_savefig}/statistics_BAz_Corr_{location}_{method}_allevents.png', dpi=300, bbox_inches='tight')
plt.show()

