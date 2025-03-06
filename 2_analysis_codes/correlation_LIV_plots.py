import os
import numpy
from scipy import signal
import matplotlib.pyplot as plt
from obspy import read, UTCDateTime
from correlation_function import correlate2

root_data = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data/miniseed_LIGO'
root_savefig = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/2_Figures/LIGO_analysis'

for folder in os.listdir(root_data):
    if folder == '.DS_Store':
        continue

    #if 'Apr' not in folder:
    #    continue

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

        for trace in [LIV_EX_DJY, LIV_EX_DHX, LIV_EX_DHY, LIV_EX_DHZ,
                      LIV_IX_DJY, LIV_IX_DHX, LIV_IX_DHY, LIV_IX_DHZ,
                      LIV_EY_DJX, LIV_EY_DHX, LIV_EY_DHY, LIV_EY_DHZ,
                      LIV_IY_DJX, LIV_IY_DHX, LIV_IY_DHY, LIV_IY_DHZ]:
            trace = trace.filter('lowpass', freq=10, zerophase=True).resample(20).detrend('simple').filter('bandpass', freqmin=0.01, freqmax = 1, zerophase=True)
            #if trace[0].stats.channel == 'DH*':
            #    #trace = trace.integrate().filter('bandpass', freqmin=0.008, freqmax = 120, zerophase=True)
            #    trace = trace.filter('bandpass', freqmin=0.008, freqmax=120, zerophase=True)

        data1 = [LIV_EY_DJX, LIV_EY_DHZ]
        data2 = [LIV_IY_DJX, LIV_IY_DHZ]
        labels12 = ['EY_DJX - IY_DJX', 'EY_DHZ - IY_DHZ']

        data3 = [LIV_EX_DJY, LIV_EX_DHZ]
        data4 = [LIV_IX_DJY, LIV_IX_DHZ]
        labels34 = ['EX_DJY - IX_DJY', 'EX_DHZ - IX_DHZ']

        data5 = [LIV_EY_DJX, LIV_EX_DJY]
        data6 = [LIV_EY_DHZ, LIV_EX_DHZ]
        labels56 = ['EY_DJX - EY_DHZ', 'EX_DJY - EX_DHZ']

        data7 = [LIV_IY_DJX, LIV_IX_DJY]
        data8 = [LIV_IY_DHZ, LIV_IX_DHZ]
        labels78 = ['IY_DJX - IY_DHZ', 'IX_DJY - IX_DHZ']

        fig, axs = plt.subplots(2,2, figsize=(11,6), sharex=True)
        plt.subplots_adjust(hspace=0, top=0.92, bottom=0.1)
        for ax in axs.flatten():
            ax.set_xlim(left=-100, right=100)
            ax.set_ylim(top=1, bottom=-1)
            ax.axvline(0, color='lightblue')

        ax = axs[0,0]
        for dataa, datab, label in zip(data1, data2, labels12):

            if 'DJ' in label:
                color='red'
            elif 'DH' in label:
                color='black'

            a = dataa.copy()[0].data
            b = datab.copy()[0].data

            N = len(dataa[0].data)
            fs = dataa[0].stats.sampling_rate
            #Cxy = signal.correlate(a, b, mode='same', method='auto')
            norm = numpy.sqrt(numpy.sum((a - numpy.mean(a)) ** 2) * numpy.sum((b - numpy.mean(b)) ** 2))
            Cxy = numpy.correlate(a, b, 'same')
            lag = numpy.linspace(-(len(a) - 1), (len(b) - 1), max(len(a), len(b))) / fs
            ax.plot(lag, Cxy/norm, color=color, label=label)
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
            ax.set_ylabel(r'$C_{IE}$')

        ax = axs[0, 1]
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

        ax = axs[1, 1]
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
        plt.suptitle('LIV_%s_all_stations, f: 0.01 - 1Hz' %folder)
        axs[0,0].set_title('interval vs external')
        axs[0,1].set_title('rotation vs translation')
        axs[1,0].set_xlabel('Lag [s]')
        axs[1,1].set_xlabel('Lag [s]')

        plt.savefig('%s/%s/LI_LIV_%s_correlation.png' %(root_savefig,folder,folder), dpi=300, bbox_inches='tight')


        #plt.show()
        #exit()

