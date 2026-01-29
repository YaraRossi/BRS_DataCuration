import obspy
from obspy import read, UTCDateTime, read_inventory
from obspy.clients.nrl import NRL
import scipy.signal as signal
from obspy.core.inventory import Inventory, Network, Station, Channel, Site

import matplotlib.pyplot as plt

root_data = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data/miniseed_LIGO'
root_xml   = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data/metadata'
folder     = '2024-04-02'

# --- Statt HAN_EX_* jetzt LIV_EX_* einlesen ---
LIV_EX_DHX = read('%s/%s/LI_LIV_EX_DHX.ms' % (root_data, folder))
LIV_EX_DHY = read('%s/%s/LI_LIV_EX_DHY.ms' % (root_data, folder))
LIV_EX_DHZ = read('%s/%s/LI_LIV_EX_DHZ.ms' % (root_data, folder))

# --- Filter & Taper wie gehabt ---
LIV_EX_DHX = LIV_EX_DHX.filter('bandpass', freqmin=0.008, freqmax=100, zerophase=True).taper(0.02)
LIV_EX_DHY = LIV_EX_DHY.filter('bandpass', freqmin=0.008, freqmax=100, zerophase=True).taper(0.02)
LIV_EX_DHZ = LIV_EX_DHZ.filter('bandpass', freqmin=0.008, freqmax=100, zerophase=True).taper(0.02)

# --- Inventory bleibt gleich (enthält HAN und LIV) ---
response = read_inventory('%s/LIGO_HAN_LIV.xml' % root_xml)

# Debug-Ausgabe
network = response[0]      # First network
station = network.stations[0]
channel = station.channels[0]
print(station.channels)

# --- Plot Setup unverändert ---
fig, axs = plt.subplots(3, 2, figsize=(12,6), sharex='col')
plt.subplots_adjust(hspace=0, wspace=0.3, right=0.98)

# Zeitreihen (roh)
axs[0,0].plot(LIV_EX_DHX[0].times('matplotlib'), LIV_EX_DHX[0].data, color='k', label='raw')
axs[1,0].plot(LIV_EX_DHY[0].times('matplotlib'), LIV_EX_DHY[0].data, color='k', label='raw')
axs[2,0].plot(LIV_EX_DHZ[0].times('matplotlib'), LIV_EX_DHZ[0].data, color='k', label='raw')

# PSD (roh)
f_t_x, P_t_x = signal.welch(LIV_EX_DHX[0].data, 1 / LIV_EX_DHX[0].stats.delta, nperseg=len(LIV_EX_DHX[0].data)//10)
axs[0,1].loglog(f_t_x, P_t_x, color='k', label='raw')
f_t_y, P_t_y = signal.welch(LIV_EX_DHY[0].data, 1 / LIV_EX_DHY[0].stats.delta, nperseg=len(LIV_EX_DHY[0].data)//10)
axs[1,1].loglog(f_t_y, P_t_y, color='k', label='raw')
f_t_z, P_t_z = signal.welch(LIV_EX_DHZ[0].data, 1 / LIV_EX_DHZ[0].stats.delta, nperseg=len(LIV_EX_DHZ[0].data)//10)
axs[2,1].loglog(f_t_z, P_t_z, color='k', label='raw')

# remove_response kopieren & anwenden
LIV_EX_DHX_rem = LIV_EX_DHX.copy()
LIV_EX_DHY_rem = LIV_EX_DHY.copy()
LIV_EX_DHZ_rem = LIV_EX_DHZ.copy()

for tr_rem in (LIV_EX_DHX_rem[0], LIV_EX_DHY_rem[0], LIV_EX_DHZ_rem[0]):
    tr_rem.remove_response(inventory=response, pre_filt=None, output="DEF", water_level=60, plot=False)

# Zeitreihen (korrigiert)
axs[0,0].plot(LIV_EX_DHX_rem[0].times('matplotlib'), LIV_EX_DHX_rem[0].data * 1e9,
              linestyle='--', color='grey', label='remove resp.')
axs[1,0].plot(LIV_EX_DHY_rem[0].times('matplotlib'), LIV_EX_DHY_rem[0].data * 1e9,
              linestyle='--', color='grey', label='remove resp.')
axs[2,0].plot(LIV_EX_DHZ_rem[0].times('matplotlib'), LIV_EX_DHZ_rem[0].data * 1e9,
              linestyle='--', color='grey', label='remove resp.')

# PSD (korrigiert)
f_tx_r, P_tx_r = signal.welch(LIV_EX_DHX_rem[0].data * 1e9, 1 / LIV_EX_DHX_rem[0].stats.delta, nperseg=len(LIV_EX_DHX_rem[0].data)//10)
axs[0,1].loglog(f_tx_r, P_tx_r, color='grey', label='remove resp.')
f_ty_r, P_ty_r = signal.welch(LIV_EX_DHY_rem[0].data * 1e9, 1 / LIV_EX_DHY_rem[0].stats.delta, nperseg=len(LIV_EX_DHY_rem[0].data)//10)
axs[1,1].loglog(f_ty_r, P_ty_r, color='grey', label='remove resp.')
f_tz_r, P_tz_r = signal.welch(LIV_EX_DHZ_rem[0].data * 1e9, 1 / LIV_EX_DHZ_rem[0].stats.delta, nperseg=len(LIV_EX_DHZ_rem[0].data)//10)
axs[2,1].loglog(f_tz_r, P_tz_r, color='grey', label='remove resp.')

# Differenzen (Zeit & PSD)
axs[0,0].plot(LIV_EX_DHX_rem[0].times('matplotlib'),
              LIV_EX_DHX_rem[0].data * 1e9 - LIV_EX_DHX[0].data,
              linestyle='-.', color='pink', label='difference')
axs[1,0].plot(LIV_EX_DHY_rem[0].times('matplotlib'),
              LIV_EX_DHY_rem[0].data * 1e9 - LIV_EX_DHY[0].data,
              linestyle='-.', color='pink', label='difference')
axs[2,0].plot(LIV_EX_DHZ_rem[0].times('matplotlib'),
              LIV_EX_DHZ_rem[0].data * 1e9 - LIV_EX_DHZ[0].data,
              linestyle='-.', color='pink', label='difference')

axs[0,1].loglog(f_tx_r, abs(P_tx_r - P_t_x), linestyle='-.', color='pink', label='difference')
axs[1,1].loglog(f_ty_r, abs(P_ty_r - P_t_y), linestyle='-.', color='pink', label='difference')
axs[2,1].loglog(f_tz_r, abs(P_tz_r - P_t_z), linestyle='-.', color='pink', label='difference')

# Achsenbeschriftungen wie gehabt
axs[0,0].xaxis_date(); axs[0,0].set_ylabel('Hx [nm/s]')
axs[1,0].xaxis_date(); axs[1,0].set_ylabel('Hy [nm/s]')
axs[2,0].xaxis_date(); axs[2,0].set_ylabel('Hz [nm/s]'); axs[2,0].set_xlabel('Time [s]')
axs[0,1].set_ylabel('PSD [(nm/s)^2/Hz]')
axs[1,1].set_ylabel('PSD [(nm/s)^2/Hz]')
axs[2,1].set_ylabel('PSD [(nm/s)^2/Hz]'); axs[2,1].set_xlabel('Freq [Hz]')

axs[0,0].legend(loc='lower right')

plt.show()
