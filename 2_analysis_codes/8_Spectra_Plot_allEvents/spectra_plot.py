import numpy
from numpy import pi
import matplotlib.pyplot as plt
import matplotlib as mpl
from obspy import read, UTCDateTime, read_inventory
from obspy.geodetics.base import gps2dist_azimuth
import scipy.signal as signal
import sys


sys.path.append(
    "/Users/yararossi/Documents/Work/Towards_Quantification/1_Coding/1.11_BRS_data_curation/BRS_DataCuration/2_analysis_codes"
)

from functions import rot_vec, circular_positive_midpoint, get_ligo_data
import pandas as pd
import os
from obspy.signal.spectral_estimation import get_nlnm, get_nhnm



# roots
root_savefig = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/2_Figures/2_Raw_PSD'
root_data = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data/miniseed_LIGO'
root_xml = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data/metadata'
response = read_inventory('%s/LIGO_HAN_LIV.xml' %root_xml)
path = '/Users/yararossi/Documents/Work/Towards_Quantification/1_Coding/1.11_BRS_data_curation/BRS_DataCuration/4_Paper_Tables'
events = pd.read_csv(f'{path}/event_times_LIV.txt')

# Import Noise models
# "PSD [dB rel. (m²/s⁴)/Hz]"
per_nlnm, nlnm = get_nlnm()
per_nhnm, nhnm = get_nhnm()

# RLNM
root_RLNM = '/Users/yararossi/Documents/Work/Random_Workstuff/RLNM/'
rlnm_df = pd.read_csv("%s/ModelForFelix_freq.csv" % root_RLNM)
rlnm_freq = rlnm_df['frequency'].values
rlnm_amp = numpy.sqrt(rlnm_df['rlnm_median (rad2/Hz)'].values*rlnm_freq)
rlnm_amp_psd = rlnm_df['rlnm_median (rad2/Hz)'].values

# Make Figure:
fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(8,4))
plt.subplots_adjust(wspace=0.3)
axs[0].plot(1/per_nhnm, nhnm, color='grey', label='NHNM') # make freq out of it
axs[0].plot(1/per_nlnm, nlnm, color='k', label='NLNM')
axs[0].set_xlim(1e-3,20)

axs[1].plot(rlnm_freq, 10 * numpy.log10(rlnm_amp_psd*1e18), color='k', label='RLNM')
axs[1].set_xlim(1e-3,20)


# import data:
# BRS
for folder in os.listdir(root_data):
    if folder in ['2022-09-19_liv','2026-01-16']: # BRS data is really bad for LIV
        continue
    if folder in ['2024-07-31','2024-09-23','2024-05-28','2024-07-02_han','2025-01-20','2023-08-20','2022-01-15_han',
                  '.DS_Store','2024-05-13_liv','2024-12-05', '2022-09-19_L','2021-07-29_L',
                  '2025-06-20_lownoise_L','2025-02-16_highnoise_L','2025-03-04_han', '2025-02-15_liv', '2023-09-08_L']:
        continue
    if 'han' in folder:
        continue

    print(folder)
    data_event = events[events['location'] == folder.split('_')[0]]

    starttime = UTCDateTime(str(data_event['eventstart'].iloc[0]))
    endtime = UTCDateTime(str(data_event['eventend'].iloc[0]))

    # theoretical baz
    source_latitude = data_event[' latitude'].values.astype(float)
    source_longitude = data_event[' longitude'].values.astype(float)
    magnitude = data_event[' mag'].values.astype(float)[0]
    norm = mpl.colors.Normalize(vmin=6, vmax=9)
    cmap = mpl.cm.get_cmap('viridis')
    color = cmap(norm(magnitude))
    station_latitude = 30.562775410180993
    station_longitude = -90.77469251757329
    baz = gps2dist_azimuth(source_latitude, source_longitude, station_latitude, station_longitude)
    scale_R = -1 # because it's rayleigh waves

    # download data
    location = 'LIV'
    LIV_EY_DJX, LIV_EX_DJY, LIV_IY_DJX, LIV_IX_DJY, LIV_IX_DH, LIV_IY_DH, LIV_EX_DH, LIV_EY_DH = get_ligo_data(folder,location)
    LIV_IY_DHX, LIV_IY_DHY, LIV_IY_DHZ = LIV_IY_DH.select(channel='*HX'), LIV_IY_DH.select(channel='*HY'), LIV_IY_DH.select(channel='*HZ')

    # rotate to EN
    channel = response.select(network='LI', station='LIV', location='IX', channel='*JY')[0][0][0]
    az = channel.azimuth/180*pi
    LIV_I_DJ = numpy.array([LIV_IY_DJX[0].data, LIV_IX_DJY[0].data, LIV_IX_DJY[0].data * 0])
    LIV_I_DJ_baz = rot_vec(0, 0, az, LIV_I_DJ)
    LIV_I_DJE_orig = LIV_IY_DJX.copy()
    LIV_I_DJN_orig = LIV_IX_DJY.copy()
    LIV_I_DJE_orig[0].data = LIV_I_DJ_baz[0, :]
    LIV_I_DJN_orig[0].data = LIV_I_DJ_baz[1, :]

    LIV_I_DJE = LIV_I_DJE_orig.copy().slice(starttime, endtime).detrend('simple').taper(0.02)
    LIV_I_DJN = LIV_I_DJN_orig.copy().slice(starttime, endtime).detrend('simple').taper(0.02)

    # Rotating seismometer
    LIV_I_DH = numpy.array([LIV_IY_DHX[0].data, LIV_IY_DHY[0].data, LIV_IY_DHZ[0].data])
    LIV_I_DH_baz = rot_vec(0, 0, az, LIV_I_DH)
    LIV_I_DHE_orig = LIV_IY_DHX.copy()
    LIV_I_DHN_orig = LIV_IY_DHY.copy()
    LIV_I_DHZ_orig = LIV_IY_DHZ.copy()
    LIV_I_DHE_orig[0].data = LIV_I_DH_baz[0, :]
    LIV_I_DHN_orig[0].data = LIV_I_DH_baz[1, :]
    LIV_I_DHZ_orig[0].data = LIV_I_DH_baz[2, :]
    LIV_IY_DHE = LIV_I_DHE_orig.copy().slice(starttime, endtime).detrend('simple').taper(0.02)
    LIV_IY_DHE[0].data = numpy.gradient(LIV_IY_DHE[0].data, LIV_IY_DHE[0].stats.delta)
    LIV_IY_DHN = LIV_I_DHN_orig.copy().slice(starttime, endtime).detrend('simple').taper(0.02)
    LIV_IY_DHN[0].data = numpy.gradient(LIV_IY_DHN[0].data, LIV_IY_DHN[0].stats.delta)
    LIV_IY_DHZ = LIV_I_DHZ_orig.copy().slice(starttime, endtime).detrend('simple').taper(0.02)
    LIV_IY_DHZ[0].data = numpy.gradient(LIV_IY_DHZ[0].data, LIV_IY_DHZ[0].stats.delta)

    ########################################################################################################################
    # Plot Data, psd_db = 10 * numpy.log10(psd)
    data_t = LIV_IY_DHZ[0]
    nperseg= 1/data_t.stats.delta*100
    nfft=nperseg
    f_t, P_t = signal.welch(data_t.data, 1 / data_t.stats.delta, nperseg=nperseg, nfft=nfft)
    axs[0].semilogx(f_t, 10 * numpy.log10(P_t), linestyle='-', color=color)

    data_r = LIV_I_DJE[0]
    nperseg= 1/data_t.stats.delta*100
    nfft=nperseg
    f_r, P_r = signal.welch(data_r.data, 1 / data_r.stats.delta, nperseg=nperseg, nfft=nfft)
    axs[1].semilogx(f_r, 10 * numpy.log10(P_r), linestyle='-', color=color)


    ########################################################################################################################
    # make labels etc
    axs[0].set_xlabel('Frequency [Hz]')
    axs[1].set_xlabel('Frequency [Hz]')
    axs[0].set_ylabel(r'PSD [$(m/s/s)^2/Hz$][dB]')
    axs[1].set_ylabel(r'PSD [$(nrad)^2/Hz$][dB]')
    axs[0].set_title('a)', loc='left')
    axs[1].set_title('b)', loc='left')

    fig = axs[1].figure
    pos = axs[1].get_position()

    cax = fig.add_axes([
        pos.x1 + 0.01,  # left
        pos.y0,  # bottom
        0.02,  # width
        pos.height  # height
    ])
    sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    plt.colorbar(sm, cax=cax, label='Event magnitude')

    print('done event')

axs[0].legend(loc='lower left')
axs[1].legend(loc='lower left')
fig.savefig(f'{root_savefig}/all_PSD_{location}.png', dpi=300, bbox_inches='tight')
plt.show()