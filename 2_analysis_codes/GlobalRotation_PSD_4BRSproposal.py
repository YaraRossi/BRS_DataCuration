import matplotlib.pyplot as plt
from obspy import read, UTCDateTime, Stream, Trace
from obspy.geodetics import gps2dist_azimuth
import pandas as pd
import scipy.signal as signal
import numpy
from functions import plot_traces_offset
from numpy import pi


root_savefig = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/2_Figures'
root_event = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data'

noise, event = True, False

# get data
root_data_ligo = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data/miniseed_LIGO'
folder = 'Dec05.2024'
# BRS
HAN_EY_DJX = read('%s/%s/LI_HAN_EY_DJX.ms' % (root_data_ligo, folder))
HAN_EX_DJY = read('%s/%s/LI_HAN_EX_DJY.ms' % (root_data_ligo, folder))
# BRS
LIV_EY_DJX = read('%s/%s/LI_LIV_EY_DJX.ms' %(root_data_ligo,folder))
LIV_EX_DJY = read('%s/%s/LI_LIV_EX_DJY.ms' %(root_data_ligo,folder))
LIV_IY_DJX = read('%s/%s/LI_LIV_IY_DJX.ms' %(root_data_ligo,folder))
LIV_IX_DJY = read('%s/%s/LI_LIV_IX_DJY.ms' %(root_data_ligo,folder))
# CRS Caltech
crs_pd = pd.read_csv('/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data/CRS_Caltech/CIT_CRS_DEC_05_DATA.csv',
                     names=['time','crs','X','Y','Z'])
trace = Trace()
trace.data = numpy.array(crs_pd['crs'], dtype=numpy.float32)*1e9/180*pi
trace.stats.starttime = UTCDateTime('2024-12-05T18:40:00')
trace.stats.sampling_rate = 512 #HAN_EY_DJX[0].stats.sampling_rate
# Set metadata
trace.stats.network = "LI"  # Adjust network code
trace.stats.station = "CAL"
trace.stats.location = ""
trace.stats.channel = "DJZ"  # Modify based on data
CAL_CRS = Stream(traces=[trace])

# BSPF
root_data_PFO = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data/PFO_data'
event_df = pd.read_csv("%s/event_times.txt" % root_event, parse_dates=['Eventtime', 'tstart', 'tend', 'location'])
tstart = event_df.loc[event_df['location'] == 'Cape Mendocino', 'tstart'].iloc[0]
starttime=UTCDateTime(tstart).timetuple()
BSPF_J1 = read('%s/PY.BSPF..HJ1.D.%4d.%03d' %(root_data_PFO,starttime.tm_year,starttime.tm_yday))
BSPF_J2 = read('%s/PY.BSPF..HJ2.D.%4d.%03d' %(root_data_PFO,starttime.tm_year,starttime.tm_yday))
BSPF_J3 = read('%s/PY.BSPF..HJ3.D.%4d.%03d' %(root_data_PFO,starttime.tm_year,starttime.tm_yday))
# ROMY Fürstenfeldbruck
root_data_romy = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data/ROMY'
romy = read('%s/romy_eventM6.8.mseed' % root_data_romy)
romy_N = romy.copy().select(channel='BJU')
for trace in romy_N:
    trace.data = trace.data*1e9


# filter data
freqmin, freqmax = 0.008, 100
BSPF_J1.detrend('linear').filter('bandpass', freqmin=freqmin, freqmax=freqmax, zerophase=True).integrate().filter('bandpass', freqmin=freqmin, freqmax=freqmax, zerophase=True)
BSPF_J2.detrend('linear').filter('bandpass', freqmin=freqmin, freqmax=freqmax, zerophase=True).integrate().filter('bandpass', freqmin=freqmin, freqmax=freqmax, zerophase=True)
BSPF_J3.detrend('linear').filter('bandpass', freqmin=freqmin, freqmax=freqmax, zerophase=True).integrate().filter('bandpass', freqmin=freqmin, freqmax=freqmax, zerophase=True)
HAN_EY_DJX.detrend('linear')#.filter('bandpass', freqmin=freqmin, freqmax=freqmax, zerophase=True)
HAN_EX_DJY.detrend('linear')#.filter('bandpass', freqmin=freqmin, freqmax=freqmax, zerophase=True)
LIV_EY_DJX.detrend('simple')#.filter('bandpass', freqmin=freqmin, freqmax=freqmax, zerophase=True)
LIV_EX_DJY.detrend('simple')#.filter('bandpass', freqmin=freqmin, freqmax=freqmax, zerophase=True)
LIV_IY_DJX.detrend('simple')#.filter('bandpass', freqmin=freqmin, freqmax=freqmax, zerophase=True)
LIV_IX_DJY.detrend('simple')#.filter('bandpass', freqmin=freqmin, freqmax=freqmax, zerophase=True)
CAL_CRS.detrend('simple').filter('lowpass', freq=125, zerophase=True).decimate(2)
romy_N.detrend('simple').filter('bandpass', freqmin=freqmin, freqmax=9.8, zerophase=True).integrate().filter('bandpass', freqmin=freqmin, freqmax=9.8, zerophase=True)



data_all = [BSPF_J1,BSPF_J2,BSPF_J3,HAN_EX_DJY,HAN_EY_DJX,LIV_EX_DJY,LIV_EX_DJY,LIV_IY_DJX,LIV_IX_DJY,CAL_CRS,romy_N]
label_all = ['BSPF_J1','BSPF_J2','BSPF_J3','HAN_EX_DJY','HAN_EY_DJX','LIV_EX_DJY','LIV_EY_DJX','LIV_IX_DJY','LIV_IY_DJX','CAL_CRS','romy_N']
color_all = ['darkblue', 'blue', 'deepskyblue', 'black', 'grey', 'darkgreen','seagreen', 'lightgreen', 'lime','red','rebeccapurple']
line_all = ['-','--','-.','-','--','-','--','-','--','-','-']

#plot_traces_offset(traces=data_all,label=label_all)
new_start = UTCDateTime('2024-12-05T18:40:00')
if noise:
    label_type = 'noise'
    start, end = UTCDateTime(new_start), UTCDateTime(new_start)+3*60 # only first 5 minutes to have noise before event.
    data_all = [data.slice(start, end) for data in data_all]
    time_label = end-start
    nperseg = 25000
if event:
    label_type = 'event'
    start, end = UTCDateTime(new_start), UTCDateTime(new_start)+20*60 # only first 20 minutes to accomodate lenth of CRS data.
    time_label = 'total length'
    nperseg = 50000

##############################################################################################################
lat1, lon1 = 40.374, -125.022
dist_PFO , azimuth, back_azimuth = gps2dist_azimuth(lat1, lon1, lat2=33.582325, lon2=-116.492470)
dist_HAN , azimuth, back_azimuth = gps2dist_azimuth(lat1, lon1, lat2=46.45511825164643, lon2=-119.40785307788653)
dist_LIV , azimuth, back_azimuth = gps2dist_azimuth(lat1, lon1, lat2=30.562775410180993, lon2=-90.77469251757329)
dist_CAL , azimuth, back_azimuth = gps2dist_azimuth(lat1, lon1, lat2=34.136726696087834, lon2=-118.12602285654852)
dist_ROMY , azimuth, back_azimuth = gps2dist_azimuth(lat1, lon1, lat2=48.16291747145096, lon2=11.275483826919913)
dist_all = [dist_PFO,dist_PFO,dist_PFO,dist_HAN,dist_HAN,dist_LIV,dist_LIV,dist_LIV,dist_LIV,dist_CAL,dist_ROMY]

##############################################################################################################
fig, ax = plt.subplots(1, 1, figsize=(10, 5))  # Increase width for space
for data, label, line, color, dist in zip(data_all, label_all, line_all, color_all, dist_all):
    if 'J3' in label:
        continue

    if 'CRS' in label:
        continue

    f_t, P_t = signal.welch(data[0].data, 1 / data[0].stats.delta, nperseg=nperseg)
    ax.loglog(f_t, P_t, color=color, label='%s @ %s km' %(label,int(round(dist/1000,0))), linestyle=line)

ax.set_xlim(left=0.008, right=300)
ax.set_ylim(top=1e6, bottom=1e-9)

ax.legend(loc="upper left", bbox_to_anchor=(1.05, 1), borderaxespad=0)  # Move legend further right
ax.set_ylabel(r'[$nrad^2/Hz$]')
ax.set_xlabel('Frequency [Hz]')
plt.subplots_adjust(right=0.75)  # Ensure space for legend
plt.suptitle('%s    --    starttime: %s, %s s' %(label_type, start, time_label))
plt.savefig('%s/Global_Rotation_PSD_%s_%s_4BRSproposal.png' % (root_savefig, folder, label_type), dpi=300, bbox_inches='tight')

plt.show()


