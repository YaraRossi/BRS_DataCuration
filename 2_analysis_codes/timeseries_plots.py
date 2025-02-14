import os

import matplotlib.pyplot as plt
from obspy import read, UTCDateTime
from correlation_function import correlate2

root_data = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data/miniseed_LIGO'
root_savefig = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/2_Figures/LIGO_analysis'

for folder in os.listdir(root_data):
    if folder == '.DS_Store':
        continue

    # make a folder to save all the figures
    newpath = r'%s/%s' % (root_savefig, folder)
    if not os.path.exists(newpath):
        os.makedirs(newpath)

    LIV, HAN = False, True

    if folder in ['July31.2024','Sept23.2024','May28.2024','July02.2024','Jan20.2025']:
        LIV = False
    if folder in ['May13.2024']:
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
            elif 'DH' in label:
                color='black'

            data_proc = data.copy().detrend('simple').filter('highpass', freq=0.008, zerophase=True)
            ax.plot(data_proc[0].times('matplotlib'), data_proc[0].data, color=color, label=label)
            ax.legend(loc='upper right')
            ax.margins(0)
            ax.xaxis_date()
        start = data_proc[0].stats.starttime
        axs[0].set_title('LIV_%s_all_stations, .detrend(simple).filter(highpass, freq=0.008, zerophase=True)' %folder)
        axs[-1].set_xlabel('Time since %s UTC [s]' %start)
        plt.savefig('%s/%s/LI_LIV_%s_all_stations.png' %(root_savefig,folder,folder), dpi=300, bbox_inches='tight')
        #plt.show()

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
            elif 'DH' in label:
                color = 'black'

            data_proc = data.copy().detrend('simple').filter('highpass', freq=0.008, zerophase=True)
            ax.plot(data_proc[0].times('matplotlib'), data_proc[0].data, color=color, label=label)
            ax.legend(loc='upper right')
            ax.margins(0)
            ax.xaxis_date()
        start = data_proc[0].stats.starttime
        axs[0].set_title('HAN_%s_all_stations, .detrend(simple).filter(highpass, freq=0.008, zerophase=True)' % folder)
        axs[-1].set_xlabel('Time since %s UTC [s]' % start)
        plt.savefig('%s/%s/LI_HAN_%s_all_stations.png' % (root_savefig, folder, folder), dpi=300, bbox_inches='tight')
        plt.show()

    '''
    
    fig, axs = plt.subplots(8,1, figsize=(11,16))
    plt.subplots_adjust(hspace=0)
    labels = ['LIV_EY_DJX - LIV_IY_DJX', 'LIV_EX_DJY - LIV_IX_DJY',
              'LIV_EY_DHX - LIV_IY_DHX', 'LIV_EY_DHY - LIV_IY_DHY', 'LIV_EY_DHZ - LIV_IY_DHZ',
              'LIV_EX_DHX - LIV_IX_DHX', 'LIV_EX_DHY - LIV_IX_DHY', 'LIV_EX_DHZ - LIV_IX_DHZ']
    
    for data1, data2, ax, label in zip([LIV_EY_DJX, LIV_EX_DJY,
                                        LIV_EY_DHX, LIV_EY_DHY, LIV_EY_DHZ,
                                        LIV_EX_DHX, LIV_EX_DHY, LIV_EX_DHZ],
                                       [LIV_IY_DJX, LIV_IX_DJY,
                                        LIV_IY_DHX, LIV_IY_DHY, LIV_IY_DHZ,
                                        LIV_IX_DHX, LIV_IX_DHY, LIV_IX_DHZ], axs, labels):
        if 'DJ' in label:
            color='red'
        elif 'DH' in label:
            color='black'
    
        data1_proc = data1.copy().detrend('simple').filter('highpass', freq=0.008, zerophase=True).slice(starttime=start, endtime=end)
        data2_proc = data2.copy().detrend('simple').filter('highpass', freq=0.008, zerophase=True).slice(starttime=start, endtime=end)

    
        ax.axvline(0, color='grey')
    
        lag, corr, lag_part, corr_part = correlate2(data1_proc[0].data, data2_proc[0].data, method='same', sps=data2_proc[0].stats.sampling_rate, onlypos=False)
        ax.plot(lag, corr, color=color, label=label)
        ax.legend(loc='upper right')
        ax.margins(0)
        print('done with: %s' % label)
    
    axs[0].set_title('LIV_April02_crosscorrelation, .detrend(simple).filter(highpass, freq=0.008, zerophase=True)')
    axs[-1].set_xlabel('Lag [s]')
    plt.savefig('%s/LI_LIV_April02_correlation_Bodywave_all_LIGO_stations.png' % root_savefig, dpi=300, bbox_inches='tight')
    
    plt.show()'''