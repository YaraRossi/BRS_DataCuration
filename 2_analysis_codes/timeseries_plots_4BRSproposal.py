import os
import pandas as pd
import matplotlib.pyplot as plt
from obspy import read, UTCDateTime
from correlation_function import correlate2

root_data = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data/miniseed_LIGO'
root_savefig = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/2_Figures/LIGO_analysis'
root_event = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data'




for folder in os.listdir(root_data):
    if folder in ['.DS_Store', '2023-08-20', '2025-02-15']:
        continue
    '''if '2024-04-02' not in folder:
       continue'''

    # make a folder to save all the figures
    newpath = r'%s/%s' % (root_savefig, folder)
    if not os.path.exists(newpath):
        os.makedirs(newpath)



    LIV, HAN = True, True

    if folder in ['2024-07-31','2024-09-23','2024-05-28','2024-07-02','2025-01-20','2022-01-15','2023-08-20']:
        LIV = False
    if folder in ['2024-05-13']:
        HAN = False

    if LIV:
        event_times = pd.read_csv('%s/event_times_LIV.txt' % root_event, parse_dates=True)

        selected_rows = event_times.loc[event_times['location'] == folder]
        starttime = UTCDateTime(str(selected_rows['eventstart'].iloc[0]))
        endtime = UTCDateTime(str(selected_rows['eventend'].iloc[0]))
        label = selected_rows['location'].iloc[0]

        # BRS
        LIV_EY_DJX = read('%s/%s/LI_LIV_EY_DJX.ms' %(root_data,folder))
        LIV_EX_DJY = read('%s/%s/LI_LIV_EX_DJY.ms' %(root_data,folder))
        LIV_IY_DJX = read('%s/%s/LI_LIV_IY_DJX.ms' %(root_data,folder))
        LIV_IX_DJY = read('%s/%s/LI_LIV_IX_DJY.ms' %(root_data,folder))

        # seis
        LIV_EY_DHX = read('%s/%s/LI_LIV_EY_DHX.ms' %(root_data,folder))
        LIV_EY_DHY = read('%s/%s/LI_LIV_EY_DHY.ms' %(root_data,folder))
        LIV_EY_DHZ = read('%s/%s/LI_LIV_EY_DHZ.ms' %(root_data,folder))
        LIV_EX_DHX = read('%s/%s/LI_LIV_EX_DHX.ms' %(root_data,folder))
        LIV_EX_DHY = read('%s/%s/LI_LIV_EX_DHY.ms' %(root_data,folder))
        LIV_EX_DHZ = read('%s/%s/LI_LIV_EX_DHZ.ms' %(root_data,folder))

        LIV_IY_DHX = read('%s/%s/LI_LIV_IY_DHX.ms' %(root_data,folder))
        LIV_IY_DHY = read('%s/%s/LI_LIV_IY_DHY.ms' %(root_data,folder))
        LIV_IY_DHZ = read('%s/%s/LI_LIV_IY_DHZ.ms' %(root_data,folder))
        LIV_IX_DHX = read('%s/%s/LI_LIV_IX_DHX.ms' %(root_data,folder))
        LIV_IX_DHY = read('%s/%s/LI_LIV_IX_DHY.ms' %(root_data,folder))
        LIV_IX_DHZ = read('%s/%s/LI_LIV_IX_DHZ.ms' %(root_data,folder))

        fig, axs = plt.subplots(8,1, figsize=(11,8), sharex=True)
        plt.subplots_adjust(hspace=0, top=0.98, bottom=0.02)
        labels = ['LIV_EY_DJX', 'LIV_EY_DHX', 'LIV_EY_DHY', 'LIV_EY_DHZ',
                  'LIV_EX_DJY', 'LIV_EX_DHX', 'LIV_EX_DHY', 'LIV_EX_DHZ',
                  'LIV_IY_DJX', 'LIV_IY_DHX', 'LIV_IY_DHY', 'LIV_IY_DHZ',
                  'LIV_IX_DJY', 'LIV_IX_DHX', 'LIV_IX_DHY', 'LIV_IX_DHZ']
        axs_new = [axs[0], axs[1], axs[1], axs[1],
                   axs[2], axs[3], axs[3], axs[3],
                   axs[4], axs[5], axs[5], axs[5],
                   axs[6], axs[7], axs[7], axs[7]]
        for data, ax, label in zip([LIV_EY_DJX, LIV_EY_DHX, LIV_EY_DHY, LIV_EY_DHZ,
                                    LIV_EX_DJY, LIV_EX_DHX, LIV_EX_DHY, LIV_EX_DHZ,
                                    LIV_IY_DJX, LIV_IY_DHX, LIV_IY_DHY, LIV_IY_DHZ,
                                    LIV_IX_DJY, LIV_IX_DHX, LIV_IX_DHY, LIV_IX_DHZ], axs_new, labels):
            if 'DJ' in label:
                color='red'
                linestyle='-'
                ax.set_ylabel('[nrad]')
            elif 'DHX' in label:
                color='black'
                linestyle='-'
                ax.set_ylabel('[nm/s]')
            elif 'DHY' in label:
                color='grey'
                linestyle='--'
            elif 'DHZ' in label:
                color='lightgrey'
                linestyle='dotted'

            data_proc = data.copy().detrend('simple').filter('bandpass', freqmin=0.008, freqmax=20, zerophase=True).slice(starttime, endtime)
            ax.plot(data_proc[0].times('matplotlib'), data_proc[0].data, color=color, label=label, linestyle=linestyle)
            ax.legend(loc='upper right')
            ax.margins(0)
            ax.xaxis_date()
        start = data_proc[0].stats.starttime
        axs[0].set_title('LIV %s, freq=0.008-20 Hz' %folder)
        axs[-1].set_xlabel('Time since %s UTC [s]' %start)
        plt.savefig('%s/%s/LI_LIV_%s_all_stations_4BRSproposal.png' %(root_savefig,folder,folder), dpi=300, bbox_inches='tight')
        #plt.show()
        #exit()

    if HAN:
        event_times = pd.read_csv('%s/event_times_HAN.txt' % root_event, parse_dates=True)

        selected_rows = event_times.loc[event_times['location'] == folder]
        starttime = UTCDateTime(str(selected_rows['eventstart'].iloc[0]))
        endtime = UTCDateTime(str(selected_rows['eventend'].iloc[0]))
        label = selected_rows['location'].iloc[0]

        # BRS
        HAN_EY_DJX = read('%s/%s/LI_HAN_EY_DJX.ms' % (root_data, folder))
        HAN_EX_DJY = read('%s/%s/LI_HAN_EX_DJY.ms' % (root_data, folder))

        # seis
        HAN_EY_DHX = read('%s/%s/LI_HAN_EY_DHX.ms' % (root_data, folder))
        HAN_EY_DHY = read('%s/%s/LI_HAN_EY_DHY.ms' % (root_data, folder))
        HAN_EY_DHZ = read('%s/%s/LI_HAN_EY_DHZ.ms' % (root_data, folder))
        HAN_EX_DHX = read('%s/%s/LI_HAN_EX_DHX.ms' % (root_data, folder))
        HAN_EX_DHY = read('%s/%s/LI_HAN_EX_DHY.ms' % (root_data, folder))
        HAN_EX_DHZ = read('%s/%s/LI_HAN_EX_DHZ.ms' % (root_data, folder))

        HAN_IY_DHX = read('%s/%s/LI_HAN_IY_DHX.ms' % (root_data, folder))
        HAN_IY_DHY = read('%s/%s/LI_HAN_IY_DHY.ms' % (root_data, folder))
        HAN_IY_DHZ = read('%s/%s/LI_HAN_IY_DHZ.ms' % (root_data, folder))
        HAN_VA_DHX = read('%s/%s/LI_HAN_VA_DHX.ms' % (root_data, folder))
        HAN_VA_DHY = read('%s/%s/LI_HAN_VA_DHY.ms' % (root_data, folder))
        #HAN_IX_DHZ = read('%s/%s/LI_HAN_IX_DHZ.ms' % (root_data, folder))

        channels = [HAN_EY_DJX, HAN_EY_DHX, HAN_EY_DHY, HAN_EY_DHZ,
                                    HAN_EX_DJY, HAN_EX_DHX, HAN_EX_DHY, HAN_EX_DHZ,
                                     HAN_IY_DHX, HAN_IY_DHY, HAN_IY_DHZ,
                                     HAN_VA_DHX, HAN_VA_DHY #, HAN_IX_DHZ
                    ]

        fig, axs = plt.subplots(6, 1, figsize=(11, 6), sharex=True)
        plt.subplots_adjust(hspace=0, top=0.98, bottom=0.02)
        labels = ['HAN_EY_DJX', 'HAN_EY_DHX', 'HAN_EY_DHY', 'HAN_EY_DHZ',
                  'HAN_EX_DJY', 'HAN_EX_DHX', 'HAN_EX_DHY', 'HAN_EX_DHZ',
                  'HAN_IY_DHX', 'HAN_IY_DHY', 'HAN_IY_DHZ',
                  'HAN_VA_DHX', 'HAN_VA_DHY' #, 'HAN_IX_DHZ'
                  ]
        axs_new = [axs[0], axs[1], axs[1], axs[1], axs[2], axs[3], axs[3], axs[3], axs[4], axs[4], axs[4], axs[5], axs[5], axs[5]]
        for data, ax, label in zip(channels, axs_new, labels):
            if 'DJ' in label:
                color = 'red'
                linestyle = '-'
                ax.set_ylabel('[nrad]')
            elif 'DHX' in label:
                color = 'black'
                linestyle = '-'
                ax.set_ylabel('[nm/s]')
            elif 'DHY' in label:
                color = 'grey'
                linestyle = '--'
            elif 'DHZ' in label:
                color = 'lightgrey'
                linestyle = 'dotted'

            data_proc = data.copy().detrend('simple').filter('bandpass', freqmin=0.008, freqmax=20, zerophase=True).slice(starttime, endtime)
            ax.plot(data_proc[0].times('matplotlib'), data_proc[0].data, color=color, label=label)
            ax.legend(loc='upper right')
            ax.margins(0)
            ax.xaxis_date()
        start = data_proc[0].stats.starttime
        axs[0].set_title('HAN %s, freq=0.008-20 Hz' %folder)
        axs[-1].set_xlabel('Time since %s UTC [s]' % start)
        plt.savefig('%s/%s/LI_HAN_%s_all_stations_4BRSproposal.png' % (root_savefig, folder, folder), dpi=300, bbox_inches='tight')
plt.show()

