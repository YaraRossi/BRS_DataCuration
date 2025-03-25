import os
import numpy
import pandas as pd
from scipy import signal
import matplotlib.pyplot as plt
from obspy import read, UTCDateTime

root_data = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data/miniseed_LIGO'
root_savefig = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/2_Figures/LIGO_analysis'
root_event = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data'


event_times = pd.read_csv('%s/event_times_HAN.txt' %root_event, parse_dates=True)
#Eventtime,tstart,tend,eventstart,eventend,location

signal_corr = True
numpy_corr = False

for folder in os.listdir(root_data):
    if folder.startswith('.'):
        continue

    if 'Sept15' not in folder:
        continue

    # make a folder to save all the figures
    newpath = r'%s/%s' % (root_savefig, folder)
    if not os.path.exists(newpath):
        os.makedirs(newpath)

    if folder in ['May13.2024']:
        continue

    selected_rows = event_times.loc[event_times['location'] == folder]
    do_slice = True
    if selected_rows.empty:
        do_slice = False
        raise ValueError(f"No matching rows found for folder: {folder}")

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

    fig, axs = plt.subplots(2,1, figsize=(6,6), sharex=True)
    plt.subplots_adjust(hspace=0, top=0.92, bottom=0.1)
    for ax in axs.flatten():
        ax.set_xlim(left=-5, right=15)
        ax.set_ylim(top=1, bottom=-1)
        ax.axvline(0, color='lightblue')

    for freqmin, freqmax, color in zip([0.01, 0.1, 1, 10], [0.1, 1, 10, 49],
                                       ['black', 'grey', 'lightgreen', 'green']):

        HAN_EY_DJX_proc = HAN_EY_DJX.copy()
        HAN_EY_DHX_proc = HAN_EY_DHX.copy()
        HAN_EY_DHY_proc = HAN_EY_DHY.copy()
        HAN_EY_DHZ_proc = HAN_EY_DHZ.copy()
        HAN_EX_DJY_proc = HAN_EX_DJY.copy()
        HAN_EX_DHX_proc = HAN_EX_DHX.copy()
        HAN_EX_DHY_proc = HAN_EX_DHY.copy()
        HAN_EX_DHZ_proc = HAN_EX_DHZ.copy()

        trace_list = [HAN_EY_DJX_proc, HAN_EY_DHX_proc, HAN_EY_DHY_proc, HAN_EY_DHZ_proc,
                      HAN_EX_DJY_proc, HAN_EX_DHX_proc, HAN_EX_DHY_proc, HAN_EX_DHZ_proc]
        processed_traces = []
        for trace in trace_list:
            # Select first matching row or adjust based on your logic
            if do_slice:
                starttime = UTCDateTime(str(selected_rows['eventstart'].iloc[0]))
                endtime = UTCDateTime(str(selected_rows['eventend'].iloc[0]))
                trace = trace.slice(starttime, endtime)

            trace.filter('lowpass', freq=49, zerophase=True).resample(100).detrend('simple').filter('bandpass',freqmin=freqmin,freqmax=freqmax,zerophase=True).taper(0.01)
            processed_traces.append(trace)

        (HAN_EY_DJX_proc, HAN_EY_DHX_proc, HAN_EY_DHY_proc, HAN_EY_DHZ_proc,
         HAN_EX_DJY_proc, HAN_EX_DHX_proc, HAN_EX_DHY_proc, HAN_EX_DHZ_proc) = processed_traces

        data5 = [HAN_EY_DJX_proc, HAN_EX_DJY_proc]
        data6 = [HAN_EY_DHZ_proc, HAN_EX_DHZ_proc]
        labels56 = ['EY_DJX - EY_DHZ, %s Hz' %[freqmin,freqmax], 'EX_DJY - EX_DHZ, %s Hz' %[freqmin,freqmax]]

        data7 = [HAN_EY_DJX_proc, HAN_EX_DJY_proc]
        data8 = [HAN_EY_DHY_proc, HAN_EX_DHX_proc]
        labels78 = ['EY_DJX - EY_DHY, %s Hz' %[freqmin,freqmax], 'EX_DJY - EX_DHX, %s Hz' %[freqmin,freqmax]]


        ax = axs[0]
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
            ax.legend(loc='upper right')
            # ax.margins(0)
            # ax.xaxis_date()
            ax.set_ylabel(r'$C_{RT}$')

        ax = axs[1]
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
            ax.legend(loc='upper right')
            # ax.margins(0)
            # ax.xaxis_date()
            ax.set_ylabel(r'$C_{RT}$')


    start = dataa[0].stats.starttime
    plt.suptitle('HAN_%s_all_stations, Filtered' %folder)
    axs[0].set_title('rotation vs translation')
    axs[1].set_xlabel('Lag [s]')

    plt.savefig('%s/%s/LI_HAN_%s_correlation.png' %(root_savefig,folder,folder), dpi=300, bbox_inches='tight')


    plt.show()
    #exit()

