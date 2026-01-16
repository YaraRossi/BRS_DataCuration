import matplotlib.pyplot as plt
from obspy import read, UTCDateTime, Stream, Trace
from obspy.geodetics import gps2dist_azimuth
import pandas as pd
import scipy.signal as signal
import numpy
import sys
sys.path.append("../2_analysis_codes")
from functions import plot_traces_offset, get_ligo_data
from numpy import pi


root_savefig = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/2_Figures/Eigenmode'
root_event = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data'

# get data
root_data_ligo = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data/miniseed_LIGO'
folder = 'Dec05.2024'
# BRS
HAN_EY_DJX, HAN_EX_DJY, _, _, _ = get_ligo_data(folder, 'HAN')
#HAN_EY_DJX = read('%s/%s/LI_HAN_EY_DJX.ms' % (root_data_ligo, folder))
#HAN_EX_DJY = read('%s/%s/LI_HAN_EX_DJY.ms' % (root_data_ligo, folder))
# BRS
LIV_EY_DJX, LIV_EX_DJY, LIV_IY_DJX, LIV_IX_DJY, _, _, _, _ = get_ligo_data(folder, 'LIV')
#LIV_EY_DJX = read('%s/%s/LI_LIV_EY_DJX.ms' %(root_data_ligo,folder))
#LIV_EX_DJY = read('%s/%s/LI_LIV_EX_DJY.ms' %(root_data_ligo,folder))
#LIV_IY_DJX = read('%s/%s/LI_LIV_IY_DJX.ms' %(root_data_ligo,folder))
#LIV_IX_DJY = read('%s/%s/LI_LIV_IX_DJY.ms' %(root_data_ligo,folder))


# filter data
freqmin, freqmax = 0.008, 100
HAN_EY_DJX.detrend('linear')#.filter('bandpass', freqmin=freqmin, freqmax=freqmax, zerophase=True)
HAN_EX_DJY.detrend('linear')#.filter('bandpass', freqmin=freqmin, freqmax=freqmax, zerophase=True)
LIV_EY_DJX.detrend('simple')#.filter('bandpass', freqmin=freqmin, freqmax=freqmax, zerophase=True)
LIV_EX_DJY.detrend('simple')#.filter('bandpass', freqmin=freqmin, freqmax=freqmax, zerophase=True)
LIV_IY_DJX.detrend('simple')#.filter('bandpass', freqmin=freqmin, freqmax=freqmax, zerophase=True)
LIV_IX_DJY.detrend('simple')#.filter('bandpass', freqmin=freqmin, freqmax=freqmax, zerophase=True)


'''data_all_orig = [HAN_EX_DJY,HAN_EY_DJX,LIV_EX_DJY,LIV_EX_DJX,LIV_IY_DJX,LIV_IX_DJY,]
label_all = ['HAN_EX_DJY','HAN_EY_DJX','LIV_EX_DJY','LIV_EY_DJX','LIV_IY_DJX','LIV_IX_DJY']
color_event = [ 'black', 'grey', 'darkgreen','seagreen', 'lightgreen', 'lime']
color_noise = [ 'darkblue', 'lightblue', 'darkorange','orange', 'gold', 'yellow']
line_event = ['-','-','-','-','-','-']
line_noise = ['-','-','-','-','-','-']'''

data_all_orig = [HAN_EX_DJY,LIV_EX_DJY,LIV_IX_DJY]
label_all = ['HAN_EX_DJY','LIV_EX_DJY','LIV_IX_DJY']
color_event = [ 'black', 'grey', 'darkgreen']
color_noise = [ 'darkorange', 'gold', 'yellow']
line_event = ['-','-','-']
line_noise = ['--','--','--']


new_start = UTCDateTime('2024-12-05T18:40:00')
# event
label_type = 'event'
start, end = UTCDateTime(new_start), UTCDateTime(new_start)+20*60 # only first 20 minutes to accomodate lenth of CRS data.
data_all_event = data_all_orig.copy()
time_label = 'total length'
nperseg_e = 50000
# noise
label_type = 'noise'
start, end = UTCDateTime(new_start), UTCDateTime(new_start)+3*60 # only first 5 minutes to have noise before event.
data_all_noise = [data.copy().slice(start, end) for data in data_all_orig]
time_label = end-start
nperseg_n = 25000

##############################################################################################################
lat1, lon1 = 40.374, -125.022
dist_HAN , azimuth, back_azimuth = gps2dist_azimuth(lat1, lon1, lat2=46.45511825164643, lon2=-119.40785307788653)
dist_LIV , azimuth, back_azimuth = gps2dist_azimuth(lat1, lon1, lat2=30.562775410180993, lon2=-90.77469251757329)
dist_all = [dist_HAN,dist_HAN,dist_LIV,dist_LIV,dist_LIV,dist_LIV]

##############################################################################################################
##############################################################################################################
# plot for poster
fig, ax = plt.subplots(1, 1, figsize=(8, 3))  # Increase width for space
for data_all, nperseg, line_all, color_all, label_type in zip([data_all_event,data_all_noise],
                                                  [nperseg_e,nperseg_n],
                                                  [line_event,line_noise],
                                                  [color_event,color_noise],
                                                  ['event', 'noise']):
    for data, label, line, color, dist in zip(data_all, label_all, line_all, color_all, dist_all):
        if 'J3' in label:
            continue

        f_t, P_t = signal.welch(data[0].data, 1 / data[0].stats.delta, nperseg=nperseg)

        ax.loglog(f_t, P_t, color=color, label=f'{label} - {label_type}', linestyle=line)


ax.set_xlim(left=0.008, right=125)
ax.set_ylim(top=1e6, bottom=1e-7)

# --- add zoomed inset in upper right corner ---
# position = [left, bottom, width, height] in figure coordinates (0-1)
ax_inset = fig.add_axes([0.55, 0.6, 0.19, 0.26])  # adjust numbers to move/resize inset

# define the frequency range to zoom in
fmin_zoom, fmax_zoom = 21, 28
ax_inset.set_xlim(fmin_zoom, fmax_zoom)
ax_inset.set_ylim(1.5e-4, 0.35)  # choose appropriate y-range for zoomed region

for data_all, nperseg, line_all, color_all, label_type in zip([data_all_event,data_all_noise],
                                                  [nperseg_e,nperseg_n],
                                                  [line_event,line_noise],
                                                  [color_event,color_noise],
                                                  ['event', 'noise']):
    for data, label, line, color, dist in zip(data_all, label_all, line_all, color_all, dist_all):
        if 'J3' in label:
            continue
        f_t, P_t = signal.welch(data[0].data, 1 / data[0].stats.delta, nperseg=nperseg)
        ax_inset.semilogy(f_t, P_t, color=color, linestyle=line)

leg = ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1), borderaxespad=0)  # Move legend further right
for line in leg.get_lines():
    line.set_linewidth(3.0)
ax.set_ylabel(r'Amplitude [$nrad^2/Hz$]')
ax.set_xlabel('Frequency [Hz]')
ax.set_title('(c)',loc='left')
plt.subplots_adjust(right=0.75)  # Ensure space for legend
plt.savefig('%s/Eigenmode_Error_PSD_%s_4poster.png' % (root_savefig, folder), dpi=300,
            bbox_inches='tight')


plt.show()


