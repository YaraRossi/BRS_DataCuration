import os
from scipy import signal
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

    LIV, HAN = True, False

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

        data1 = [LIV_EY_DJX, LIV_EY_DHY]
        data2 = [LIV_IY_DJX, LIV_IY_DHY]
        labels12 = ['EY_DJX - IY_DJX', 'EY_DHY - IY_DHY']

        data3 = [LIV_EX_DJY, LIV_EX_DHY]
        data4 = [LIV_IX_DJY, LIV_IX_DHY]
        labels34 = ['EX_DJY - IX_DJY', 'EX_DHX - IX_DHX']

        data5 = [LIV_EY_DJX, LIV_EX_DJY]
        data6 = [LIV_EY_DHY, LIV_EX_DHX]
        labels56 = ['EY_DJX - EY_DHY', 'EX_DJY - EX_DHX']

        data7 = [LIV_IY_DJX, LIV_IX_DJY]
        data8 = [LIV_IY_DHY, LIV_IX_DHX]
        labels78 = ['IY_DJX - IY_DHY', 'IX_DJY - IX_DHX']

        fig, axs = plt.subplots(2,2, figsize=(11,6), sharex=True)
        plt.subplots_adjust(hspace=0, top=0.92, bottom=0.1)

        ax = axs[0,0]
        for dataa, datab, label in zip(data1, data2, labels12):

            if 'DJ' in label:
                color='red'
            elif 'DH' in label:
                color='black'

            #a = dataa.copy().detrend('simple').filter('highpass', freq=0.008, zerophase=True)[0].data
            #b = datab.copy().detrend('simple').filter('highpass', freq=0.008, zerophase=True)[0].data

            a = dataa.copy().detrend('simple')[0].data
            b = datab.copy().detrend('simple')[0].data
            N = len(dataa[0].data)
            fs = dataa[0].stats.sampling_rate
            f, Cxy = signal.coherence(a, b, fs=fs, window='hann', nperseg=int(N / 10), nfft=None, detrend=False,
                                      axis=-1)
            ax.semilogx(f, Cxy, color=color, label=label)
            ax.legend(loc='upper right')
            #ax.margins(0)
            #ax.xaxis_date()
            ax.set_ylabel(r'$C_{IE}$')

        ax = axs[1,0]
        for dataa, datab, label in zip(data3, data4, labels34):

            if 'DJ' in label:
                color = 'red'
            elif 'DH' in label:
                color = 'black'

            #a = dataa.copy().detrend('simple').filter('highpass', freq=0.008, zerophase=True)[0].data
            #b = datab.copy().detrend('simple').filter('highpass', freq=0.008, zerophase=True)[0].data

            a = dataa.copy().detrend('simple')[0].data
            b = datab.copy().detrend('simple')[0].data

            N = len(dataa[0].data)
            fs = dataa[0].stats.sampling_rate
            f, Cxy = signal.coherence(a, b, fs=fs, window='hann', nperseg=int(N / 10), nfft=None,
                                      detrend=False,
                                      axis=-1)
            ax.semilogx(f, Cxy, color=color, label=label)
            ax.legend(loc='upper right')
            # ax.margins(0)
            # ax.xaxis_date()
            ax.set_ylabel(r'$C_{IE}$')

        ax = axs[0, 1]
        for dataa, datab, label, linestyle, color in zip(data5, data6, labels56, ['-','--'], ['black', 'grey']):

            a = dataa.copy().detrend('simple')[0].data
            b = datab.copy().detrend('simple')[0].data
            N = len(dataa[0].data)
            fs = dataa[0].stats.sampling_rate
            f, Cxy = signal.coherence(a, b, fs=fs, window='hann', nperseg=int(N / 10), nfft=None, detrend=False,
                                      axis=-1)
            ax.semilogx(f, Cxy, color=color, label=label)
            ax.legend(loc='upper right')
            # ax.margins(0)
            # ax.xaxis_date()
            ax.set_ylabel(r'$C_{RT}$')

        ax = axs[1, 1]
        for dataa, datab, label, linestyle, color in zip(data7, data8, labels78, ['-','--'], ['black', 'grey']):

            a = dataa.copy().detrend('simple')[0].data
            b = datab.copy().detrend('simple')[0].data

            N = len(dataa[0].data)
            fs = dataa[0].stats.sampling_rate
            f, Cxy = signal.coherence(a, b, fs=fs, window='hann', nperseg=int(N / 10), nfft=None, detrend=False,
                                      axis=-1)
            ax.semilogx(f, Cxy, color=color, label=label)
            ax.legend(loc='upper right')
            # ax.margins(0)
            # ax.xaxis_date()
            ax.set_ylabel(r'$C_{RT}$')


        start = dataa[0].stats.starttime
        plt.suptitle('LIV_%s_all_stations, filter(highpass, freq=0.008)' %folder)
        axs[0,0].set_title('interval vs external')
        axs[0,1].set_title('rotation vs translation')
        axs[1,0].set_xlabel('Frequency [Hz]')
        axs[1,1].set_xlabel('Frequency [Hz]')

        plt.savefig('%s/%s/LI_LIV_%s_coherence.png' %(root_savefig,folder,folder), dpi=300, bbox_inches='tight')


        #plt.show()

