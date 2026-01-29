import os
from obspy import read, UTCDateTime, read_inventory
import pandas as pd
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


# import data:
# BRS
for LIV, HAN, station in zip([True, False], [False, True], ['LIV', 'HAN']):
    j = 0
    for folder in os.listdir(root_data):
        if folder in ['.DS_Store', '2025-02-15_liv', '2022-01-15_han', '2024-07-02_han', '2024-05-13_liv',
                      '2024-12-05', '2025-03-04_han']:
            continue
        if 'han' in folder and LIV==True:
            continue
        if 'liv' in folder and HAN==True:
            continue

        if LIV:
            if folder in ['2024-07-31','2024-09-23','2024-05-28','2024-07-02',
                          '2025-01-20','2023-08-20','2022-01-15']:
                continue


        print(folder)
        events = pd.read_csv(f'{path}/event_times_{station}.txt')
        data_event = events[events['location'] == folder.split('_')[0]]
        starttime = UTCDateTime(str(data_event['tstart'].iloc[0]))
        eventstart = UTCDateTime(str(data_event['eventstart'].iloc[0]))

        if LIV:
            LIV_EY_DJX = read('%s/%s/LI_LIV_EY_DJX.ms' % (root_data, folder))
            LIV_EX_DJY = read('%s/%s/LI_LIV_EX_DJY.ms' % (root_data, folder))
            LIV_IY_DJX = read('%s/%s/LI_LIV_IY_DJX.ms' % (root_data, folder))
            LIV_IX_DJY = read('%s/%s/LI_LIV_IX_DJY.ms' % (root_data, folder))

            LIV_EY_DHX = read('%s/%s/LI_LIV_EY_DHX.ms' % (root_data, folder))

        if HAN:
            HAN_EY_DJX = read('%s/%s/LI_HAN_EY_DJX.ms' % (root_data, folder))
            HAN_EX_DJY = read('%s/%s/LI_HAN_EX_DJY.ms' % (root_data, folder))


        db_bins_r = (-20, 70,1)
        db_bins_t = (-135,-50,1)
        p_s_w_o = 0.1# 0.1
        p_s_o = 0.007
        ppsd_length = 600

        if '_L' in folder:
            response = response_L
        else:
            response = response_H


        if j == 0:
            if LIV:
                r_ppsd_LIV_EY_DJX = PPSD(LIV_EY_DJX[0].stats, metadata=response, db_bins=db_bins_r,
                                period_step_octaves=p_s_o, period_smoothing_width_octaves=p_s_w_o, ppsd_length=ppsd_length)
                r_ppsd_LIV_EX_DJY = PPSD(LIV_EX_DJY[0].stats, metadata=response, db_bins=db_bins_r,
                                period_step_octaves=p_s_o, period_smoothing_width_octaves=p_s_w_o, ppsd_length=ppsd_length)
                r_ppsd_LIV_IY_DJX = PPSD(LIV_IY_DJX[0].stats, metadata=response, db_bins=db_bins_r,
                                period_step_octaves=p_s_o, period_smoothing_width_octaves=p_s_w_o, ppsd_length=ppsd_length)
                r_ppsd_LIV_IX_DJY = PPSD(LIV_IX_DJY[0].stats, metadata=response, db_bins=db_bins_r,
                                period_step_octaves=p_s_o, period_smoothing_width_octaves=p_s_w_o, ppsd_length=ppsd_length)
                t_ppsd_LIV_EY_DHX = PPSD(LIV_EY_DHX[0].stats, metadata=response, db_bins=db_bins_t,
                                period_step_octaves=p_s_o, period_smoothing_width_octaves=p_s_w_o, ppsd_length=ppsd_length)
                j = 1
            if HAN:
                r_ppsd_HAN_EY_DJX = PPSD(HAN_EY_DJX[0].stats, metadata=response, db_bins=db_bins_r,
                                         period_step_octaves=p_s_o, period_smoothing_width_octaves=p_s_w_o,
                                         ppsd_length=ppsd_length)
                r_ppsd_HAN_EX_DJY = PPSD(HAN_EX_DJY[0].stats, metadata=response, db_bins=db_bins_r,
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
        fig = r_ppsd_LIV_EY_DJX.plot(show=False, cmap=cmap, xaxis_frequency=freq)
        axes = fig.axes
        axes[0].set_ylabel("Amplitude [nrad^2/Hz] [dB]")
        fig.savefig('%s/PPSD/PPSD_LIV_EY_DJX.png' % (root_savefig), dpi=300)

        fig = r_ppsd_LIV_EX_DJY.plot(show=False, cmap=cmap, xaxis_frequency=freq)
        axes = fig.axes
        axes[0].set_ylabel("Amplitude [nrad^2/Hz] [dB]")
        fig.savefig('%s/PPSD/PPSD_LIV_EX_DJY.png' % (root_savefig), dpi=300)

        fig = r_ppsd_LIV_IY_DJX.plot(show=False, cmap=cmap, xaxis_frequency=freq)
        axes = fig.axes
        axes[0].set_ylabel("Amplitude [nrad^2/Hz] [dB]")
        fig.savefig('%s/PPSD/PPSD_LIV_IY_DJX.png' % (root_savefig), dpi=300)

        fig = r_ppsd_LIV_IX_DJY.plot(show=False, cmap=cmap, xaxis_frequency=freq)
        axes = fig.axes
        axes[0].set_ylabel("Amplitude [nrad^2/Hz] [dB]")
        fig.savefig('%s/PPSD/PPSD_LIV_IX_DJY.png' % (root_savefig), dpi=300)

        fig = t_ppsd_LIV_EY_DHX.plot(show=False, cmap=cmap, xaxis_frequency=freq)
        fig.savefig('%s/PPSD/PPSD_LIV_EY_DHX.png' % (root_savefig), dpi=300)

    if HAN:
        fig = r_ppsd_HAN_EY_DJX.plot(show=False, cmap=cmap, xaxis_frequency=freq)
        axes = fig.axes
        axes[0].set_ylabel("Amplitude [nrad^2/Hz] [dB]")
        fig.savefig('%s/PPSD/PPSD_HAN_EY_DJX.png' % (root_savefig), dpi=300)

        fig = r_ppsd_HAN_EX_DJY.plot(show=False, cmap=cmap, xaxis_frequency=freq)
        axes = fig.axes
        axes[0].set_ylabel("Amplitude [nrad^2/Hz] [dB]")
        fig.savefig('%s/PPSD/PPSD_HAN_EX_DJY.png' % (root_savefig), dpi=300)