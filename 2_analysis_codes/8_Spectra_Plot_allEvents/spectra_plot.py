import numpy
from numpy import pi
import matplotlib.pyplot as plt
import matplotlib as mpl
from obspy import read, UTCDateTime, read_inventory
from obspy.signal.cross_correlation import correlate
from obspy.signal.rotate import rotate_ne_rt
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
root_savefig = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/2_Figures/Backazimuth_estimations/fixed_freq'
root_savedata = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data/correlation_baz_csv'
root_data = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data/miniseed_LIGO'
root_xml = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data/metadata'
response = read_inventory('%s/LIGO_HAN_LIV.xml' %root_xml)
path = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data'
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
axs[0].plot(1/per_nlnm, nlnm, color='grey')
axs[0].plot(1/per_nhnm, nhnm, color='grey') # make freq out of it
axs[0].set_xlim(1e-3,270)

axs[1].loglog(rlnm_freq, rlnm_amp_psd, color='grey')
axs[0].set_xlim(1e-3,270)


# import data:
# BRS
for folder in os.listdir(root_data):
    if folder in ['July31.2024','Sept23.2024','May28.2024','July02.2024_han','Jan20.2025','Aug20.2023','Jan15.2022_han',
                  '.DS_Store','May13.2024_liv','Dec05.2024', 'Sept19.2022_L','July29.2021_L',
                  'Jun20.2025_lownoise_L','Feb16.2025_highnoise_L','Mar04.2025_han', 'Feb15.2025_liv']:
        continue
    if 'han' in folder:
        continue
    data_event = events[events['location'] == folder.split('_')[0]]

    starttime = UTCDateTime(str(data_event['eventstart'].iloc[0]))
    endtime = UTCDateTime(str(data_event['eventend'].iloc[0]))

    print(folder)

    # theoretical baz
    source_latitude = data_event[' latitude'].values.astype(float)
    source_longitude = data_event[' longitude'].values.astype(float)
    station_latitude = 30.562775410180993
    station_longitude = -90.77469251757329
    baz = gps2dist_azimuth(source_latitude, source_longitude, station_latitude, station_longitude)
    scale_R = -1 # because it's rayleigh waves

    # download data
    location = 'LIV'
    LIV_EY_DJX, LIV_EX_DJY, LIV_IY_DJX, LIV_IX_DJY, LIV_IX_DH, LIV_IY_DH, LIV_EX_DH, LIV_EY_DH = get_ligo_data(folder,location)
    LIV_IY_DHX, LIV_IY_DHY, LIV_IY_DHZ = LIV_IY_DH.select(channel='DHX'), LIV_IY_DH.select(channel='DHY'), LIV_IY_DH.select(channel='DHZ')

    # rotate to EN
    channel = response.select(network='LI', station='LIV', location='IX', channel='DJY')[0][0][0]
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
    nperseg= 4048
    data_t = LIV_IY_DHZ[0]
    f_t, P_t = signal.welch(data_t.data, 1 / data_t.stats.delta, nperseg=int(len(data_t)/4))
    axs[0].semilogx(f_t, 10 * numpy.log10(P_t), linestyle='-')

    data_r = LIV_I_DJE[0]
    f_r, P_r = signal.welch(data_r.data, 1 / data_r.stats.delta, nperseg=int(len(data_t)/4))
    axs[1].loglog(f_r, P_r, linestyle='-')
    print('done event')

plt.show()