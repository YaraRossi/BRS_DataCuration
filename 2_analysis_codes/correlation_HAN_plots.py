import os
import numpy
from scipy import signal
import matplotlib.pyplot as plt
from functions import get_ligo_data
from obspy import read, UTCDateTime
from correlation_function import correlate2

root_data = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data/miniseed_LIGO'
root_savefig = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/2_Figures/LIGO_analysis'
freqmin, freqmax = 0.01,3
for folder in os.listdir(root_data):
    if folder == '.DS_Store':
        continue

    #if 'Apr' not in folder:
    #    continue

    # make a folder to save all the figures
    newpath = r'%s/%s' % (root_savefig, folder)
    if not os.path.exists(newpath):
        os.makedirs(newpath)

    LIV, HAN = False, True

    if folder in ['2024-07-31','2024-09-23','2024-05-28','2024-07-02','2025-01-20','2023-08-20','2022-01-15']:
        LIV = False
    if folder in ['2024-05-13', '2025-02-15']:
        HAN = False

    if HAN:
        location = 'HAN'
        HAN_EY_DJX, HAN_EX_DJY, HAN_EY_DH, HAN_EX_DH, HAN_VA_DH = get_ligo_data(folder, location)
        HAN_EY_DHX, HAN_EY_DHY, HAN_EY_DHZ = HAN_EY_DH.select(channel='DHX'), HAN_EY_DH.select(channel='DHY'),HAN_EY_DH.select(channel='DHZ')
        HAN_EX_DHX, HAN_EX_DHY, HAN_EX_DHZ = HAN_EX_DH.select(channel='DHX'), HAN_EX_DH.select(channel='DHY'),HAN_EX_DH.select(channel='DHZ')
        '''# BRS
        HAN_EY_DJX = read('%s/%s/LI_HAN_EY_DJX.ms' % (root_data, folder))
        HAN_EX_DJY = read('%s/%s/LI_HAN_EX_DJY.ms' % (root_data, folder))

        # seis
        HAN_EY_DHX = read('%s/%s/LI_HAN_EY_DHX.ms' % (root_data, folder))
        HAN_EY_DHY = read('%s/%s/LI_HAN_EY_DHY.ms' % (root_data, folder))
        HAN_EY_DHZ = read('%s/%s/LI_HAN_EY_DHZ.ms' % (root_data, folder))
        HAN_EX_DHX = read('%s/%s/LI_HAN_EX_DHX.ms' % (root_data, folder))
        HAN_EX_DHY = read('%s/%s/LI_HAN_EX_DHY.ms' % (root_data, folder))
        HAN_EX_DHZ = read('%s/%s/LI_HAN_EX_DHZ.ms' % (root_data, folder))

        #HAN_IY_DHX = read('%s/%s/LI_HAN_IY_DHX.ms' % (root_data, folder))
        #HAN_IY_DHY = read('%s/%s/LI_HAN_IY_DHY.ms' % (root_data, folder))
        #HAN_IY_DHZ = read('%s/%s/LI_HAN_IY_DHZ.ms' % (root_data, folder))
        HAN_VA_DHX = read('%s/%s/LI_HAN_VA_DHX.ms' % (root_data, folder))
        HAN_VA_DHY = read('%s/%s/LI_HAN_VA_DHY.ms' % (root_data, folder))'''

        for trace in [HAN_EY_DJX, HAN_EY_DHX, HAN_EY_DHY, HAN_EY_DHZ,
                      HAN_EX_DJY, HAN_EX_DHX, HAN_EX_DHY, HAN_EX_DHZ]:
            trace = trace.detrend('demean').filter('lowpass', freq=10, zerophase=True).resample(20).filter('bandpass', freqmin=freqmin, freqmax=freqmax, zerophase=True)

        data5 = [HAN_EY_DJX, HAN_EX_DJY]
        data6 = [HAN_EY_DHZ, HAN_EX_DHZ]
        labels56 = ['EY_DJX - EY_DHZ', 'EX_DJY - EX_DHZ']

        data7 = [HAN_EY_DJX, HAN_EX_DJY]
        data8 = [HAN_EY_DHY, HAN_EX_DHX]
        labels78 = ['EY_DJX - EY_DHY', 'EX_DJY - EX_DHX']

        fig, axs = plt.subplots(2,1, figsize=(6,6), sharex=True)
        plt.subplots_adjust(hspace=0, top=0.92, bottom=0.1)
        for ax in axs.flatten():
            ax.set_xlim(left=-100, right=100)
            ax.set_ylim(top=1, bottom=-1)
            ax.axvline(0, color='lightblue')

        ax = axs[0]
        for dataa, datab, label, linestyle, color in zip(data5, data6, labels56, ['-','--'], ['black', 'grey']):

            a = dataa.copy()[0].data
            b = datab.copy()[0].data

            N = len(dataa[0].data)
            fs = dataa[0].stats.sampling_rate
            #Cxy = signal.correlate(a, b, mode='same', method='auto')
            norm = numpy.sqrt(numpy.sum((a - numpy.mean(a)) ** 2) * numpy.sum((b - numpy.mean(b)) ** 2))
            Cxy = numpy.correlate(a, b, 'same')
            lag = numpy.linspace(-(len(a) - 1), (len(b) - 1), max(len(a), len(b))) / fs
            ax.plot(lag, Cxy / norm, color=color, label=label)
            ax.legend(loc='upper right')
            # ax.margins(0)
            # ax.xaxis_date()
            ax.set_ylabel(r'$C_{RT}$')

        ax = axs[1]
        for dataa, datab, label, linestyle, color in zip(data7, data8, labels78, ['-','--'], ['black', 'grey']):

            a = dataa.copy()[0].data
            b = datab.copy()[0].data

            N = len(dataa[0].data)
            fs = dataa[0].stats.sampling_rate
            # Cxy = signal.correlate(a, b, mode='same', method='auto')
            norm = numpy.sqrt(numpy.sum((a - numpy.mean(a)) ** 2) * numpy.sum((b - numpy.mean(b)) ** 2))
            Cxy = numpy.correlate(a, b, 'same')
            lag = numpy.linspace(-(len(a) - 1), (len(b) - 1), max(len(a), len(b))) / fs
            ax.plot(lag, Cxy / norm, color=color, label=label)
            ax.legend(loc='upper right')
            # ax.margins(0)
            # ax.xaxis_date()
            ax.set_ylabel(r'$C_{RT}$')


        start = dataa[0].stats.starttime
        plt.suptitle(f'{location}_%s_all_stations, f: {freqmin} - {freqmax}Hz' %folder)
        axs[0].set_title('rotation vs translation')
        axs[1].set_xlabel('Lag [s]')

        plt.savefig(f'%s/Correlation/LI_{location}_%s_correlation.png' %(root_savefig,folder), dpi=300, bbox_inches='tight')


        #plt.show()
        #exit()

