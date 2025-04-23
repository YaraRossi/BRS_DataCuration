import obspy
from obspy import read, UTCDateTime
import matplotlib.pyplot as plt
import pandas as pd
import numpy
from scipy import signal
from numpy import pi, sin, cos


# timeseries overlap rot vs trans. for specific freq range.

root_data = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data/miniseed_LIGO'
root_event = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data'

event_times = pd.read_csv('%s/event_times_LIV.txt' %root_event, parse_dates=True)
#Eventtime,tstart,tend,eventstart,eventend,location
folder = 'Sept15.2024'
folder = 'Apr02.2024'

selected_rows = event_times.loc[event_times['location'] == folder]
starttime = UTCDateTime(str(selected_rows['eventstart'].iloc[0]))
endtime = UTCDateTime(str(selected_rows['eventend'].iloc[0]))

# location eq:
lat = 51.609
lon = -130.612
#


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

# filter, slice
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
    starttime = UTCDateTime(str(selected_rows['eventstart'].iloc[0]))
    endtime = UTCDateTime(str(selected_rows['eventend'].iloc[0]))
    trace.filter('lowpass', freq=49,
                 zerophase=True).resample(100).detrend('simple').filter('bandpass',freqmin=0.01,freqmax=0.2,zerophase=True)
    trace = trace.slice(starttime, endtime)
    processed_traces.append(trace)

(LIV_EX_DJY_proc, LIV_EX_DHX_proc, LIV_EX_DHY_proc, LIV_EX_DHZ_proc,
 LIV_IX_DJY_proc, LIV_IX_DHX_proc, LIV_IX_DHY_proc, LIV_IX_DHZ_proc,
 LIV_EY_DJX_proc, LIV_EY_DHX_proc, LIV_EY_DHY_proc, LIV_EY_DHZ_proc,
 LIV_IY_DJX_proc, LIV_IY_DHX_proc, LIV_IY_DHY_proc, LIV_IY_DHZ_proc) = processed_traces

########################################################################################################################
# making time series figure
fig,axs = plt.subplots(3,1, figsize=(11,6), sharex=True)
plt.subplots_adjust(hspace=0)

axs0 = axs[0].twinx()
axs1 = axs[1].twinx()

axs[0].plot(LIV_IX_DHZ_proc[0].times('matplotlib'), LIV_IX_DHZ_proc[0].data, color='k', label='IX_DHZ')
axs0.plot(LIV_IX_DJY_proc[0].times('matplotlib'), LIV_IX_DJY_proc[0].data, color='red', label='IX_DJY')

axs[1].plot(LIV_IY_DHZ_proc[0].times('matplotlib'), LIV_IY_DHZ_proc[0].data, color='k', label='IY_DHZ')
axs1.plot(LIV_IY_DJX_proc[0].times('matplotlib'), LIV_IY_DJX_proc[0].data, color='red', label='IY_DJX')

axs[2].plot(LIV_IY_DHZ_proc[0].times('matplotlib'), LIV_IY_DHZ_proc[0].data, color='k', label='IY_DHZ')
axs[2].plot(LIV_IX_DHZ_proc[0].times('matplotlib'), LIV_IX_DHZ_proc[0].data, color='grey', linestyle='--', label='IY_DHZ')

# Set y-limits for the primary and secondary axes to be symmetric around zero
for ax in axs:
    ymin, ymax = ax.get_ylim()
    ax.set_ylim([-max(abs(ymin), abs(ymax)), max(abs(ymin), abs(ymax))])
    ax.legend(loc='upper left')
    ax.set_ylabel('[nm/s]')
    ax.xaxis_date()
    ax.margins(0)

for ax in [axs0, axs1]:
    ymin, ymax = ax.get_ylim()
    ax.set_ylim([-max(abs(ymin), abs(ymax)), max(abs(ymin), abs(ymax))])
    ax.legend(loc='lower left')
    ax.set_ylabel('[nrad]')
    ax.xaxis_date()
    ax.margins(0)

axs[1].set_xlabel('Time [s]')

########################################################################################################################
# Backazimuth estimation

fig, axs = plt.subplots(2,1, figsize=(11,5), sharex=True)
plt.subplots_adjust(hspace=0)
angles = numpy.linspace(0,360,36)
max_C = []
max_lag = []

for theta_deg in angles:
    theta = theta_deg/180*pi
    x_data = LIV_IY_DJX_proc[0].data
    y_data = LIV_IX_DJY_proc[0].data
    z_data = LIV_IY_DHZ_proc[0].data
    new_N = y_data * cos(theta) - x_data * sin(theta)
    new_E = y_data * sin(theta) + x_data * cos(theta)
    a = new_E
    b = z_data
    N = len(x_data)
    fs = LIV_IY_DJX_proc[0].stats.sampling_rate
    Cxy = signal.correlate(a, b, mode='same', method='auto')
    norm = numpy.sqrt(numpy.sum((a - numpy.mean(a)) ** 2) * numpy.sum((b - numpy.mean(b)) ** 2))
    lag = numpy.linspace(-N // 2, N // 2, len(Cxy)) / fs

    max_C.append(max(Cxy/norm))
    max_lag.append(lag[numpy.argmax(Cxy)])


axs[0].scatter(angles, max_C)
axs[1].scatter(angles, max_lag)
    
plt.show()