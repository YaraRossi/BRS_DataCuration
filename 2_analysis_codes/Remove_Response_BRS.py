import obspy
import numpy
from obspy import read, UTCDateTime, read_inventory
from obspy.clients.nrl import NRL
import scipy.signal as signal
from obspy.core.inventory import Inventory, Network, Station, Channel, Site

import matplotlib.pyplot as plt

root_data = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data/miniseed_LIGO'
root_xml = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data/metadata'
folder = '2024-04-02'

HAN_EY_DJX = read('%s/%s/LI_HAN_EY_DJX.ms' % (root_data, folder))
HAN_EX_DJY = read('%s/%s/LI_HAN_EX_DJY.ms' % (root_data, folder))

HAN_EY_DJX = HAN_EY_DJX.filter('bandpass', freqmin=0.008, freqmax = 100, zerophase=True).taper(0.02)
HAN_EX_DJY = HAN_EX_DJY.filter('bandpass', freqmin=0.008, freqmax = 100, zerophase=True).taper(0.02)

response = read_inventory('%s/LIGO_HAN_LIV.xml' %root_xml)

network = response[0]  # First network
station = network.stations[0]  # First station
channel = station.channels[0]  # First channel

# Print the response
print(channel.response)

fig, axs = plt.subplots(2,2, figsize=(9,6))
plt.subplots_adjust(hspace=0, wspace=0.3, right=0.98)

axs[0,0].plot(HAN_EY_DJX[0].times('matplotlib'), HAN_EY_DJX[0].data, color='k', label='raw')
axs[1,0].plot(HAN_EX_DJY[0].times('matplotlib'), HAN_EX_DJY[0].data, color='k', label='raw')

f_t_jx, P_t_jx = signal.welch(HAN_EY_DJX[0].data, 1 / HAN_EY_DJX[0].stats.delta, nperseg=len(HAN_EY_DJX[0].data)/10)
axs[0,1].loglog(f_t_jx, P_t_jx, color='k', label='raw')
f_t_jy, P_t_jy = signal.welch(HAN_EX_DJY[0].data, 1 / HAN_EX_DJY[0].stats.delta, nperseg=len(HAN_EX_DJY[0].data)/10)
axs[1,1].loglog(f_t_jy, P_t_jy, color='k', label='raw')

HAN_EY_DJX_rem = HAN_EY_DJX.copy()
HAN_EX_DJY_rem = HAN_EX_DJY.copy()
#pre_filt = [0.008, 0.009, 100.0, 125.0]
HAN_EY_DJX_rem[0].remove_response(inventory=response, pre_filt=None, output="DEF", water_level=60, plot=False)
HAN_EX_DJY_rem[0].remove_response(inventory=response, pre_filt=None, output="DEF", water_level=60, plot=False)

axs[0,0].plot(HAN_EY_DJX_rem[0].times('matplotlib'), HAN_EY_DJX_rem[0].data, linestyle='--', color='grey', label='remove resp.')
axs[1,0].plot(HAN_EX_DJY_rem[0].times('matplotlib'), HAN_EX_DJY_rem[0].data, linestyle='--', color='grey', label='remove resp.')

f_t_jx_rem, P_t_jx_rem = signal.welch(HAN_EY_DJX_rem[0].data, 1 / HAN_EY_DJX_rem[0].stats.delta, nperseg=len(HAN_EY_DJX_rem[0].data)/10)
axs[0,1].loglog(f_t_jx_rem, P_t_jx_rem, color='grey', label='remove resp.')
f_t_jy_rem, P_t_jy_rem = signal.welch(HAN_EX_DJY_rem[0].data, 1 / HAN_EX_DJY_rem[0].stats.delta, nperseg=len(HAN_EX_DJY_rem[0].data)/10)
axs[1,1].loglog(f_t_jy_rem, P_t_jy_rem, color='grey', label='remove resp.')


axs[0,0].plot(HAN_EY_DJX_rem[0].times('matplotlib'), HAN_EY_DJX_rem[0].data-HAN_EY_DJX[0].data, linestyle='-.', color='pink', label='difference')
axs[1,0].plot(HAN_EX_DJY_rem[0].times('matplotlib'), HAN_EX_DJY_rem[0].data-HAN_EX_DJY[0].data, linestyle='-.', color='pink', label='difference')

axs[0,1].loglog(f_t_jx_rem, P_t_jx_rem-P_t_jx, linestyle='-.', color='pink', label='difference')
axs[1,1].loglog(f_t_jy_rem, P_t_jy_rem-P_t_jy, linestyle='-.', color='pink', label='difference')

# correlation:
a = HAN_EY_DJX.copy()[0].data[50000:60000]
b = HAN_EY_DJX_rem.copy()[0].data[50000:60000]

fs = HAN_EY_DJX[0].stats.sampling_rate
norm = numpy.sqrt(numpy.sum((a - numpy.mean(a)) ** 2) * numpy.sum((b - numpy.mean(b)) ** 2))
Cxy = numpy.correlate(a, b, 'same')
lag = numpy.linspace(-(len(a) - 1), (len(b) - 1), max(len(a), len(b))) / fs

index = numpy.where(Cxy == max(Cxy))[0][0]

axs[0,0].set_title('Cxy: %s, lag: %s s' %(round(Cxy[index]/norm,3), round(lag[index],4)))

axs[0,0].xaxis_date()
axs[0,0].set_ylabel('Jx [nrad]')
axs[1,0].xaxis_date()
axs[1,0].set_ylabel('Jy [nrad]')
axs[1,0].set_xlabel('Time [s]')
axs[0,0].legend(loc = 'lower right')

axs[0,1].set_ylabel('PSD [nrad^2/Hz]')
axs[1,1].set_ylabel('PSD [nrad^2/Hz]')
axs[1,1].set_xlabel('Freq [Hz]')


plt.show()
