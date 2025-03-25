import os
import numpy
import pandas as pd
from scipy import signal
import matplotlib.pyplot as plt
from obspy import read, UTCDateTime

root_data = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data/miniseed_LIGO'
root_savefig = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/2_Figures/LIGO_analysis'
root_event = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data'


event_times = pd.read_csv('%s/event_times_LIV.txt' %root_event, parse_dates=True)
#Eventtime,tstart,tend,eventstart,eventend,location

signal_corr = True
numpy_corr = False

for folder in os.listdir(root_data):
    if folder == '.DS_Store':
        continue

    '''if 'Apr' not in folder:
        continue'''

    # make a folder to save all the figures
    newpath = r'%s/%s' % (root_savefig, folder)
    if not os.path.exists(newpath):
        os.makedirs(newpath)

    if folder in ['July31.2024','Sept23.2024','May28.2024','July02.2024','Jan20.2025','Aug20.2023','Jan15.2022']:
        continue

    selected_rows = event_times.loc[event_times['location'] == folder]
    do_slice = True
    if selected_rows.empty:
        do_slice = False
        print(f"No matching rows found for folder: {folder}")

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


    fig, axs = plt.subplots(2,2, figsize=(11,6), sharex=True)
    plt.subplots_adjust(hspace=0, top=0.92, bottom=0.1)
    for ax in axs.flatten():
        ax.set_xlim(left=-5, right=15)
        ax.set_ylim(top=1, bottom=-1)
        ax.axvline(0, color='lightblue')


    for freqmin, freqmax, color in zip([0.01,0.1,1,10],[0.1,1,10, 49], ['black', 'grey', 'lightgreen', 'green']):
        LIV_EX_DJY_proc, LIV_EX_DHX_proc, LIV_EX_DHY_proc, LIV_EX_DHZ_proc = LIV_EX_DJY.copy(), LIV_EX_DHX.copy(), LIV_EX_DHY.copy(), LIV_EX_DHZ.copy()
        LIV_IX_DJY_proc, LIV_IX_DHX_proc, LIV_IX_DHY_proc, LIV_IX_DHZ_proc = LIV_IX_DJY.copy(), LIV_IX_DHX.copy(), LIV_IX_DHY.copy(), LIV_IX_DHZ.copy()
        LIV_EY_DJX_proc, LIV_EY_DHX_proc, LIV_EY_DHY_proc, LIV_EY_DHZ_proc = LIV_EY_DJX.copy(), LIV_EY_DHX.copy(), LIV_EY_DHY.copy(), LIV_EY_DHZ.copy()
        LIV_IY_DJX_proc, LIV_IY_DHX_proc, LIV_IY_DHY_proc, LIV_IY_DHZ_proc = LIV_IY_DJX.copy(), LIV_IY_DHX.copy(), LIV_IY_DHY.copy(), LIV_IY_DHZ.copy()

        trace_list = [LIV_EX_DJY_proc, LIV_EX_DHX_proc, LIV_EX_DHY_proc, LIV_EX_DHZ_proc,
                      LIV_IX_DJY_proc, LIV_IX_DHX_proc, LIV_IX_DHY_proc, LIV_IX_DHZ_proc,
                      LIV_EY_DJX_proc, LIV_EY_DHX_proc, LIV_EY_DHY_proc, LIV_EY_DHZ_proc,
                      LIV_IY_DJX_proc, LIV_IY_DHX_proc, LIV_IY_DHY_proc, LIV_IY_DHZ_proc]
        processed_traces = []
        for trace in trace_list:
            # Select first matching row or adjust based on your logic
            if do_slice:
                starttime = UTCDateTime(str(selected_rows['eventstart'].iloc[0]))
                endtime = UTCDateTime(str(selected_rows['eventend'].iloc[0]))
                trace = trace.slice(starttime, endtime)

            trace.filter('lowpass',freq=49,zerophase=True).resample(100).detrend('simple').filter('bandpass',freqmin=freqmin,freqmax=freqmax, zerophase=True)
            processed_traces.append(trace)

        (LIV_EX_DJY_proc, LIV_EX_DHX_proc, LIV_EX_DHY_proc, LIV_EX_DHZ_proc,
        LIV_IX_DJY_proc, LIV_IX_DHX_proc, LIV_IX_DHY_proc, LIV_IX_DHZ_proc,
        LIV_EY_DJX_proc, LIV_EY_DHX_proc, LIV_EY_DHY_proc, LIV_EY_DHZ_proc,
        LIV_IY_DJX_proc, LIV_IY_DHX_proc, LIV_IY_DHY_proc, LIV_IY_DHZ_proc) = processed_traces

        data1 = [LIV_EY_DJX_proc, LIV_EY_DHZ_proc]
        data2 = [LIV_IY_DJX_proc, LIV_IY_DHZ_proc]
        labels12 = ['EY_DJX - IY_DJX, %s Hz' %[freqmin,freqmax], 'EY_DHZ - IY_DHZ, %s Hz' %[freqmin,freqmax]]

        data3 = [LIV_EX_DJY_proc, LIV_EX_DHZ_proc]
        data4 = [LIV_IX_DJY_proc, LIV_IX_DHZ_proc]
        labels34 = ['EX_DJY - IX_DJY, %s Hz' %[freqmin,freqmax], 'EX_DHZ - IX_DHZ, %s Hz' %[freqmin,freqmax]]

        data5 = [LIV_EY_DJX_proc, LIV_EX_DJY_proc]
        data6 = [LIV_EY_DHZ_proc, LIV_EX_DHZ_proc]
        labels56 = ['EY_DJX - EY_DHZ, %s Hz' %[freqmin,freqmax], 'EX_DJY - EX_DHZ, %s Hz' %[freqmin,freqmax]]

        data7 = [LIV_IY_DJX_proc, LIV_IX_DJY_proc]
        data8 = [LIV_IY_DHZ_proc, LIV_IX_DHZ_proc]
        labels78 = ['IY_DJX - IY_DHZ, %s Hz' %[freqmin,freqmax], 'IX_DJY - IX_DHZ, %s Hz' %[freqmin,freqmax]]


        ax = axs[0,0]
        for dataa, datab, label, linestyle in zip(data1, data2, labels12, ['-','--']):

            a = dataa.copy()[0].data
            b = datab.copy()[0].data

            N = len(dataa[0].data)
            fs = dataa[0].stats.sampling_rate
            if signal_corr:
                Cxy = signal.correlate(a, b, mode='same', method='auto')
                lag = numpy.linspace(-N // 2, N // 2, len(Cxy)) / fs

            norm = numpy.sqrt(numpy.sum((a - numpy.mean(a)) ** 2) * numpy.sum((b - numpy.mean(b)) ** 2))

            if numpy_corr:
                Cxy = numpy.correlate(a, b, 'same')
                lag = numpy.linspace(-(len(a) - 1), (len(b) - 1), max(len(a), len(b))) / fs
            ax.plot(lag, Cxy/norm, color=color, label=label, linestyle=linestyle)
            ax.legend(loc='upper right', fontsize='x-small')
            #ax.margins(0)
            #ax.xaxis_date()
            ax.set_ylabel(r'$C_{IE}$')

        ax = axs[1,0]
        for dataa, datab, label, linestyle in zip(data3, data4, labels34, ['-','--']):

            '''if 'DJ' in label:
                color = 'red'
            elif 'DH' in label:
                color = 'black' '''

            a = dataa.copy()[0].data
            b = datab.copy()[0].data

            N = len(dataa[0].data)
            fs = dataa[0].stats.sampling_rate
            if signal_corr:
                Cxy = signal.correlate(a, b, mode='same', method='auto')
                lag = numpy.linspace(-N // 2, N // 2, len(Cxy)) / fs

            norm = numpy.sqrt(numpy.sum((a - numpy.mean(a)) ** 2) * numpy.sum((b - numpy.mean(b)) ** 2))

            if numpy_corr:
                Cxy = numpy.correlate(a, b, 'same')
                lag = numpy.linspace(-(len(a) - 1), (len(b) - 1), max(len(a), len(b))) / fs
            ax.plot(lag, Cxy / norm, color=color, label=label, linestyle=linestyle)
            ax.legend(loc='upper right', fontsize='x-small')
            # ax.margins(0)
            # ax.xaxis_date()
            ax.set_ylabel(r'$C_{IE}$')

        ax = axs[0, 1]
        for dataa, datab, label, linestyle, color1 in zip(data5, data6, labels56, ['-','--'], ['black', 'grey']):

            a = dataa.copy()[0].data
            b = datab.copy()[0].data

            N = len(dataa[0].data)
            fs = dataa[0].stats.sampling_rate
            if signal_corr:
                Cxy = signal.correlate(a, b, mode='same', method='auto')
                lag = numpy.linspace(-N // 2, N // 2, len(Cxy)) / fs

            norm = numpy.sqrt(numpy.sum((a - numpy.mean(a)) ** 2) * numpy.sum((b - numpy.mean(b)) ** 2))

            if numpy_corr:
                Cxy = numpy.correlate(a, b, 'same')
                lag = numpy.linspace(-(len(a) - 1), (len(b) - 1), max(len(a), len(b))) / fs
            ax.plot(lag, Cxy / norm, color=color, label=label, linestyle=linestyle)
            ax.legend(loc='upper right', fontsize='x-small')
            # ax.margins(0)
            # ax.xaxis_date()
            ax.set_ylabel(r'$C_{RT}$')

        ax = axs[1, 1]
        for dataa, datab, label, linestyle, color1 in zip(data7, data8, labels78, ['-','--'], ['black', 'grey']):

            a = dataa.copy()[0].data
            b = datab.copy()[0].data

            N = len(dataa[0].data)
            fs = dataa[0].stats.sampling_rate
            if signal_corr:
                Cxy = signal.correlate(a, b, mode='same', method='auto')
                lag = numpy.linspace(-N // 2, N // 2, len(Cxy)) / fs

            norm = numpy.sqrt(numpy.sum((a - numpy.mean(a)) ** 2) * numpy.sum((b - numpy.mean(b)) ** 2))

            if numpy_corr:
                Cxy = numpy.correlate(a, b, 'same')
                lag = numpy.linspace(-(len(a) - 1), (len(b) - 1), max(len(a), len(b))) / fs
            ax.plot(lag, Cxy / norm, color=color, label=label, linestyle=linestyle)
            ax.legend(loc='upper right', fontsize='x-small')
            # ax.margins(0)
            # ax.xaxis_date()
            ax.set_ylabel(r'$C_{RT}$')


        start = dataa[0].stats.starttime
        plt.suptitle('LIV_%s_all_stations, Filtered' %folder)
        axs[0,0].set_title('interval vs external')
        axs[0,1].set_title('rotation vs translation')
        axs[1,0].set_xlabel('Time [s]')
        axs[1,1].set_xlabel('Time [s]')

    plt.savefig('%s/%s/LI_LIV_%s_correlation.png' %(root_savefig,folder,folder), dpi=300, bbox_inches='tight')




    #plt.show()
    #exit()

