import os
from obspy import read, UTCDateTime, read_inventory
import pandas as pd
import numpy
from obspy.signal import PPSD
from obspy.imaging.cm import pqlx, viridis
import matplotlib.pyplot as plt

# roots
root_savefig = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/2_Figures/'
root_data = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data/zenodo_upload/miniseed_LIGO'
root_xml = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data/metadata'
response_H = read_inventory('%s/LIGO_HAN_LIV.xml' %root_xml)
response_L = read_inventory('%s/LIGO_HAN_LIV_lowsps.xml' %root_xml)
path = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data'

# get RLNM
root_RLNM = '/Users/yararossi/Documents/Work/Random_Workstuff/RLNM/'
rlnm_df = pd.read_csv("%s/ModelForFelix_freq.csv" % root_RLNM)
rlnm_freq = rlnm_df['frequency'].values
rlnm_amp = numpy.sqrt(rlnm_df['rlnm_median (rad2/Hz)'].values*rlnm_freq)
rlnm_amp_psd = rlnm_df['rlnm_median (rad2/Hz)'].values

freq_rlnm = rlnm_freq
rlnm_db = 10 * numpy.log10(rlnm_amp_psd*1e18) # make nrad,


# import data:
# BRS
for LIV, HAN, station in zip([True, False], [False, True], ['LIV', 'HAN']):
    j = 0
    for folder in ['2025-06-20_lownoise_L','2025-02-16_highnoise_L']:

        print(folder)
        events = pd.read_csv(f'{path}/event_times_{station}.txt')
        data_event = events[events['location'] == folder.split('_')[0]]
        starttime = UTCDateTime(str(data_event['tstart'].iloc[0]))
        eventstart = UTCDateTime(str(data_event['eventstart'].iloc[0]))

        if LIV:
            LIV_EY_DJX = read('%s/%s/LI_LIV_EY_*JX.ms' % (root_data, folder))
            LIV_EX_DJY = read('%s/%s/LI_LIV_EX_*JY.ms' % (root_data, folder))
            LIV_IY_DJX = read('%s/%s/LI_LIV_IY_*JX.ms' % (root_data, folder))
            LIV_IX_DJY = read('%s/%s/LI_LIV_IX_*JY.ms' % (root_data, folder))

            LIV_EY_DHX = read('%s/%s/LI_LIV_EY_*HX.ms' % (root_data, folder))

        if HAN:
            HAN_EY_DJX = read('%s/%s/LI_HAN_EY_*JX.ms' % (root_data, folder))
            HAN_EX_DJY = read('%s/%s/LI_HAN_EX_*JY.ms' % (root_data, folder))


        db_bins_r = (-55, 55,1)
        db_bins_t = (-185,-75,1)
        p_s_w_o = 0.4# 0.1
        p_s_o = 0.03#0.007
        ppsd_length = 2400

        if '_L' in folder:
            response = response_L
        else:
            response = response_H


        if j == 0:
            if LIV:
                r_ppsd_LIV_EY_DJX = PPSD(LIV_EY_DJX[0].stats, metadata=None, db_bins=db_bins_r,
                                period_step_octaves=p_s_o, period_smoothing_width_octaves=p_s_w_o, ppsd_length=ppsd_length)
                r_ppsd_LIV_EX_DJY = PPSD(LIV_EX_DJY[0].stats, metadata=None, db_bins=db_bins_r,
                                period_step_octaves=p_s_o, period_smoothing_width_octaves=p_s_w_o, ppsd_length=ppsd_length)
                r_ppsd_LIV_IY_DJX = PPSD(LIV_IY_DJX[0].stats, metadata=None, db_bins=db_bins_r,
                                period_step_octaves=p_s_o, period_smoothing_width_octaves=p_s_w_o, ppsd_length=ppsd_length)
                r_ppsd_LIV_IX_DJY = PPSD(LIV_IX_DJY[0].stats, metadata=None, db_bins=db_bins_r,
                                period_step_octaves=p_s_o, period_smoothing_width_octaves=p_s_w_o, ppsd_length=ppsd_length)
                t_ppsd_LIV_EY_DHX = PPSD(LIV_EY_DHX[0].stats, metadata=None, db_bins=db_bins_t,
                                period_step_octaves=p_s_o, period_smoothing_width_octaves=p_s_w_o, ppsd_length=ppsd_length)
                j = 1
            if HAN:
                r_ppsd_HAN_EY_DJX = PPSD(HAN_EY_DJX[0].stats, metadata=None, db_bins=db_bins_r,
                                         period_step_octaves=p_s_o, period_smoothing_width_octaves=p_s_w_o,
                                         ppsd_length=ppsd_length)
                r_ppsd_HAN_EX_DJY = PPSD(HAN_EX_DJY[0].stats, metadata=None, db_bins=db_bins_r,
                                         period_step_octaves=p_s_o, period_smoothing_width_octaves=p_s_w_o,
                                         ppsd_length=ppsd_length)
                j = 1

        # add stream to ppsd
        if LIV:
            r_ppsd_LIV_EY_DJX.add(LIV_EY_DJX[0])
            r_ppsd_LIV_EX_DJY.add(LIV_EX_DJY[0])
            r_ppsd_LIV_IY_DJX.add(LIV_IY_DJX[0])
            r_ppsd_LIV_IX_DJY.add(LIV_IX_DJY[0])
            t_ppsd_LIV_EY_DHX.add(LIV_EY_DHX[0])
        if HAN:
            r_ppsd_HAN_EY_DJX.add(HAN_EY_DJX[0])
            r_ppsd_HAN_EX_DJY.add(HAN_EX_DJY[0])


    cmap = pqlx
    freq = True
    if LIV:
        fig = r_ppsd_LIV_EY_DJX.plot(show=False, cmap=cmap, xaxis_frequency=freq, period_lim=(0.003, 6))
        axes = fig.axes
        axes[0].set_ylabel("Amplitude [nrad^2/Hz] [dB]")
        axes[0].plot(freq_rlnm, rlnm_db, 'grey', lw=2, label='Rotational LNM')
        axes[0].set_title("")
        axes[0].set_title('(c) LIV.EY.BJX', loc='left')
        fig.savefig('%s/PPSD/Noise_PPSD_LIV_EY_DJX.png' % (root_savefig), dpi=300)

        fig = r_ppsd_LIV_EX_DJY.plot(show=False, cmap=cmap, xaxis_frequency=freq, period_lim=(0.003, 6))
        axes = fig.axes
        axes[0].set_ylabel("Amplitude [nrad^2/Hz] [dB]")
        axes[0].plot(freq_rlnm, rlnm_db, 'grey', lw=2, label='Rotational LNM')
        axes[0].set_title("")
        axes[0].set_title('(d) LIV.EX.BJY', loc='left')
        fig.savefig('%s/PPSD/Noise_PPSD_LIV_EX_DJY.png' % (root_savefig), dpi=300)

        fig = r_ppsd_LIV_IY_DJX.plot(show=False, cmap=cmap, xaxis_frequency=freq, period_lim=(0.003, 6))
        axes = fig.axes
        axes[0].set_ylabel("Amplitude [nrad^2/Hz] [dB]")
        axes[0].plot(freq_rlnm, rlnm_db, 'grey', lw=2, label='Rotational LNM')
        axes[0].set_title("")
        axes[0].set_title('(e) LIV.IY.BJX', loc='left')
        fig.savefig('%s/PPSD/Noise_PPSD_LIV_IY_DJX.png' % (root_savefig), dpi=300)

        fig = r_ppsd_LIV_IX_DJY.plot(show=False, cmap=cmap, xaxis_frequency=freq, period_lim=(0.003, 6))
        axes = fig.axes
        axes[0].set_ylabel("Amplitude [nrad^2/Hz] [dB]")
        axes[0].plot(freq_rlnm, rlnm_db, 'grey', lw=2, label='Rotational LNM')
        axes[0].set_title("")
        axes[0].set_title('(f) LIV.IX.BJY', loc='left')
        fig.savefig('%s/PPSD/Noise_PPSD_LIV_IX_DJY.png' % (root_savefig), dpi=300)

        fig = t_ppsd_LIV_EY_DHX.plot(show=False, cmap=cmap, xaxis_frequency=freq, period_lim=(0.003, 6))
        fig.savefig('%s/PPSD/Noise_PPSD_LIV_EY_DHX.png' % (root_savefig), dpi=300)

    if HAN:
        fig = r_ppsd_HAN_EY_DJX.plot(show=False, cmap=cmap, xaxis_frequency=freq, period_lim=(0.003, 6))
        axes = fig.axes
        axes[0].set_ylabel("Amplitude [nrad^2/Hz] [dB]")
        axes[0].plot(freq_rlnm, rlnm_db, 'grey', lw=2, label='Rotational LNM')
        axes[0].set_title("")
        axes[0].set_title('(a) HAN.EY.BJX', loc='left')
        fig.savefig('%s/PPSD/Noise_PPSD_HAN_EY_DJX.png' % (root_savefig), dpi=300)

        fig = r_ppsd_HAN_EX_DJY.plot(show=False, cmap=cmap, xaxis_frequency=freq, period_lim=(0.003, 6))
        axes = fig.axes
        axes[0].set_ylabel("Amplitude [nrad^2/Hz] [dB]")
        axes[0].plot(freq_rlnm, rlnm_db, 'grey', lw=2, label='Rotational LNM')
        axes[0].set_title("")
        axes[0].set_title('(b) HAN.EX.BJY', loc='left')
        fig.savefig('%s/PPSD/Noise_PPSD_HAN_EX_DJY.png' % (root_savefig), dpi=300)