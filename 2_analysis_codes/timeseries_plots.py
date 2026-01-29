import os

import matplotlib.pyplot as plt
from obspy import read, UTCDateTime
from correlation_function import correlate2

root_data = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data/miniseed_LIGO'
root_savefig = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/2_Figures/LIGO_analysis'

for folder in os.listdir(root_data):
    if folder in ['.DS_Store', '2023-08-20', '2025-02-15']:
        continue
    '''if '2024-04-02' not in folder:
       continue'''

    # make a folder to save all the figures
    newpath = r'%s/%s' % (root_savefig, folder)
    if not os.path.exists(newpath):
        os.makedirs(newpath)

    LIV, HAN = True, False

    if folder in ['2024-07-31','2024-09-23','2024-05-28','2024-07-02','2025-01-20','2022-01-15','2023-08-20']:
        LIV = False
    if folder in ['2024-05-13']:
        HAN = False

    if LIV:
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

        fig, axs = plt.subplots(16,1, figsize=(11,16), sharex=True)
        plt.subplots_adjust(hspace=0, top=0.98, bottom=0.02)
        labels = ['LIV_EY_DJX', 'LIV_EY_DHX', 'LIV_EY_DHY', 'LIV_EY_DHZ',
                  'LIV_EX_DJY', 'LIV_EX_DHX', 'LIV_EX_DHY', 'LIV_EX_DHZ',
                  'LIV_IY_DJX', 'LIV_IY_DHX', 'LIV_IY_DHY', 'LIV_IY_DHZ',
                  'LIV_IX_DJY', 'LIV_IX_DHX', 'LIV_IX_DHY', 'LIV_IX_DHZ']

        for data, ax, label in zip([LIV_EY_DJX, LIV_EY_DHX, LIV_EY_DHY, LIV_EY_DHZ,
                                    LIV_EX_DJY, LIV_EX_DHX, LIV_EX_DHY, LIV_EX_DHZ,
                                    LIV_IY_DJX, LIV_IY_DHX, LIV_IY_DHY, LIV_IY_DHZ,
                                    LIV_IX_DJY, LIV_IX_DHX, LIV_IX_DHY, LIV_IX_DHZ], axs, labels):
            if 'DJ' in label:
                color='red'
                ax.set_ylabel('[nrad]')
            elif 'DH' in label:
                color='black'
                ax.set_ylabel('[nm/s]')

            data_proc = data.copy().detrend('simple').filter('highpass', freq=0.008, zerophase=True)#.filter('bandpass',freqmin=0.008,freqmax=100, zerophase=True)
            ax.plot(data_proc[0].times('matplotlib'), data_proc[0].data, color=color, label=label)
            ax.legend(loc='upper right')
            ax.margins(0)
            ax.xaxis_date()
        start = data_proc[0].stats.starttime
        axs[0].set_title('LIV_%s_all_stations: highpass, freq=0.008' %folder)
        axs[-1].set_xlabel('Time since %s UTC [s]' %start)
        plt.savefig('%s/%s/LI_LIV_%s_all_stations.png' %(root_savefig,folder,folder), dpi=300, bbox_inches='tight')
        plt.show()
        #exit()

    if HAN:
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

        fig, axs = plt.subplots(len(channels), 1, figsize=(11, len(channels)), sharex=True)
        plt.subplots_adjust(hspace=0, top=0.98, bottom=0.02)
        labels = ['HAN_EY_DJX', 'HAN_EY_DHX', 'HAN_EY_DHY', 'HAN_EY_DHZ',
                  'HAN_EX_DJY', 'HAN_EX_DHX', 'HAN_EX_DHY', 'HAN_EX_DHZ',
                  'HAN_IY_DHX', 'HAN_IY_DHY', 'HAN_IY_DHZ',
                  'HAN_VA_DHX', 'HAN_VA_DHY' #, 'HAN_IX_DHZ'
                  ]

        for data, ax, label in zip(channels, axs, labels):
            if 'DJ' in label:
                color = 'red'
                ax.set_ylabel('[nrad]')
            elif 'DH' in label:
                color = 'black'
                ax.set_ylabel('[nm/s]')

            data_proc = data.copy().detrend('simple').filter('highpass', freq=0.008, zerophase=True)
            ax.plot(data_proc[0].times('matplotlib'), data_proc[0].data, color=color, label=label)
            ax.legend(loc='upper right')
            ax.margins(0)
            ax.xaxis_date()
        start = data_proc[0].stats.starttime
        axs[0].set_title('HAN_%s_all_stations, highpass, freq=0.008' % folder)
        axs[-1].set_xlabel('Time since %s UTC [s]' % start)
        plt.savefig('%s/%s/LI_HAN_%s_all_stations.png' % (root_savefig, folder, folder), dpi=300, bbox_inches='tight')
        plt.show()

