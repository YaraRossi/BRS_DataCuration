import obspy
from obspy import read, UTCDateTime, read_inventory
from obspy.clients.nrl import NRL
import scipy.signal as signal
from obspy.core.inventory import Inventory, Network, Station, Channel, Site

import matplotlib.pyplot as plt

root_data = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data/miniseed_LIGO'
root_xml = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data/metadata'
folder = '2024-04-02'

HAN_EX_DHX = read('%s/%s/LI_HAN_EX_DHX.ms' % (root_data, folder))
HAN_EX_DHY = read('%s/%s/LI_HAN_EX_DHY.ms' % (root_data, folder))
HAN_EX_DHZ = read('%s/%s/LI_HAN_EX_DHZ.ms' % (root_data, folder))


HAN_EX_DHX = HAN_EX_DHX.filter('bandpass', freqmin=0.08, freqmax = 100, zerophase=True).taper(0.02)
HAN_EX_DHY = HAN_EX_DHY.filter('bandpass', freqmin=0.08, freqmax = 100, zerophase=True).taper(0.02)
HAN_EX_DHZ = HAN_EX_DHZ.filter('bandpass', freqmin=0.08, freqmax = 100, zerophase=True).taper(0.02)

response = read_inventory('%s/LIGO_HAN_LIV.xml' %root_xml)

network = response[0]  # First network
station = network.stations[0]  # First station
channel = station.channels[0]  # First channel

# Print the response
print(station.channels)

fig, axs = plt.subplots(3,2, figsize=(12,6), sharex='col')
plt.subplots_adjust(hspace=0, wspace=0.3, right=0.98)


axs[0,0].plot(HAN_EX_DHX[0].times('matplotlib'), HAN_EX_DHX[0].data, color='k', label='raw')
axs[1,0].plot(HAN_EX_DHY[0].times('matplotlib'), HAN_EX_DHY[0].data, color='k', label='raw')
axs[2,0].plot(HAN_EX_DHZ[0].times('matplotlib'), HAN_EX_DHZ[0].data, color='k', label='raw')

f_t_jx, P_t_jx = signal.welch(HAN_EX_DHX[0].data, 1 / HAN_EX_DHX[0].stats.delta, nperseg=len(HAN_EX_DHX[0].data)/10)
axs[0,1].loglog(f_t_jx, P_t_jx, color='k', label='raw')
f_t_jy, P_t_jy = signal.welch(HAN_EX_DHY[0].data, 1 / HAN_EX_DHY[0].stats.delta, nperseg=len(HAN_EX_DHY[0].data)/10)
axs[1,1].loglog(f_t_jy, P_t_jy, color='k', label='raw')
f_t_jz, P_t_jz = signal.welch(HAN_EX_DHZ[0].data, 1 / HAN_EX_DHZ[0].stats.delta, nperseg=len(HAN_EX_DHZ[0].data)/10)
axs[2,1].loglog(f_t_jz, P_t_jz, color='k', label='raw')

HAN_EX_DHX_rem = HAN_EX_DHX.copy()
HAN_EX_DHY_rem = HAN_EX_DHY.copy()
HAN_EX_DHZ_rem = HAN_EX_DHZ.copy()
HAN_EX_DHX_rem[0].remove_response(inventory=response, pre_filt=None, output="DEF", water_level=60, plot=False)
HAN_EX_DHY_rem[0].remove_response(inventory=response, pre_filt=None, output="DEF", water_level=60, plot=False)
HAN_EX_DHZ_rem[0].remove_response(inventory=response, pre_filt=None, output="DEF", water_level=60, plot=False)

axs[0,0].plot(HAN_EX_DHX_rem[0].times('matplotlib'), HAN_EX_DHX_rem[0].data*1e9, linestyle='--', color='grey', label='remove resp.')
axs[1,0].plot(HAN_EX_DHY_rem[0].times('matplotlib'), HAN_EX_DHY_rem[0].data*1e9, linestyle='--', color='grey', label='remove resp.')
axs[2,0].plot(HAN_EX_DHZ_rem[0].times('matplotlib'), HAN_EX_DHZ_rem[0].data*1e9, linestyle='--', color='grey', label='remove resp.')

f_t_jx_rem, P_t_jx_rem = signal.welch(HAN_EX_DHX_rem[0].data*1e9, 1 / HAN_EX_DHX_rem[0].stats.delta, nperseg=len(HAN_EX_DHX_rem[0].data)/10)
axs[0,1].loglog(f_t_jx_rem, P_t_jx_rem, color='grey', label='remove resp.')
f_t_jy_rem, P_t_jy_rem = signal.welch(HAN_EX_DHY_rem[0].data*1e9, 1 / HAN_EX_DHY_rem[0].stats.delta, nperseg=len(HAN_EX_DHY_rem[0].data)/10)
axs[1,1].loglog(f_t_jy_rem, P_t_jy_rem, color='grey', label='remove resp.')
f_t_jz_rem, P_t_jz_rem = signal.welch(HAN_EX_DHZ_rem[0].data*1e9, 1 / HAN_EX_DHZ_rem[0].stats.delta, nperseg=len(HAN_EX_DHZ_rem[0].data)/10)
axs[2,1].loglog(f_t_jz_rem, P_t_jz_rem, color='grey', label='remove resp.')


axs[0,0].plot(HAN_EX_DHX_rem[0].times('matplotlib'), HAN_EX_DHX_rem[0].data*1e9-HAN_EX_DHX[0].data, linestyle='-.', color='pink', label='difference')
axs[1,0].plot(HAN_EX_DHY_rem[0].times('matplotlib'), HAN_EX_DHY_rem[0].data*1e9-HAN_EX_DHY[0].data, linestyle='-.', color='pink', label='difference')
axs[2,0].plot(HAN_EX_DHZ_rem[0].times('matplotlib'), HAN_EX_DHZ_rem[0].data*1e9-HAN_EX_DHZ[0].data, linestyle='-.', color='pink', label='difference')

axs[0,1].loglog(f_t_jx_rem, abs(P_t_jx_rem-P_t_jx), linestyle='-.', color='pink', label='difference')
axs[1,1].loglog(f_t_jy_rem, abs(P_t_jy_rem-P_t_jy), linestyle='-.', color='pink', label='difference')
axs[2,1].loglog(f_t_jz_rem, abs(P_t_jz_rem-P_t_jz), linestyle='-.', color='pink', label='difference')


axs[0,0].xaxis_date()
axs[0,0].set_ylabel('Hx [nm/s]')
axs[1,0].xaxis_date()
axs[1,0].set_ylabel('Hy [nm/s]')
axs[2,0].xaxis_date()
axs[2,0].set_ylabel('Hz [nm/s]')
axs[2,0].set_xlabel('Time [s]')
axs[0,0].legend(loc = 'lower right')

axs[0,1].set_ylabel('PSD [(nm/s)^2/Hz]')
axs[1,1].set_ylabel('PSD [(nm/s)^2/Hz]')
axs[2,1].set_ylabel('PSD [(nm/s)^2/Hz]')
axs[2,1].set_xlabel('Freq [Hz]')


plt.show()
