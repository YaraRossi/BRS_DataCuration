import os
import pandas as pd
import numpy
import time
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.colorbar import ColorbarBase
from matplotlib.legend_handler import HandlerTuple
from numpy import sin, cos, pi
from obspy import read, UTCDateTime, read_inventory
from correlation_function import correlate2
from obspy.geodetics.base import gps2dist_azimuth
from functions import rot_vec, supylabel2
from datetime import datetime
tic = time.time()
root_data = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data/miniseed_LIGO'
root_savefig = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/2_Figures/RotationError'
root_xml = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data/metadata'
response = read_inventory('%s/LIGO_HAN_LIV.xml' %root_xml)
path = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data'

fig_tele_re, axs_tele_re = plt.subplots(1, 1, figsize=(5, 5))
plt.subplots_adjust(hspace=0, top=0.95, bottom=0.05)
fig_tele_am, axs_tele_am = plt.subplots(1, 1, figsize=(5, 5))
plt.subplots_adjust(hspace=0, top=0.95, bottom=0.05)

fig_local_re, axs_local_re = plt.subplots(1, 1, figsize=(5, 5))
plt.subplots_adjust(hspace=0, top=0.95, bottom=0.05)
fig_local_am, axs_local_am = plt.subplots(1, 1, figsize=(5, 5))
plt.subplots_adjust(hspace=0, top=0.95, bottom=0.05)

marker_type = ['*', 'o', '<', '^', '>', 'v', 'p','s','X','d','P', 'h','D', '8']
n_marker = 0
labels_tele =[]
labels_local =[]
handles_tele = []
handles_local = []

# magnitude > 5:
freqmin_b = [0.008, 0.02, 0.04, 0.06, 0.08, 0.1]
freqmax_b = [0.02, 0.04, 0.06, 0.08, 0.1, 49]
color_freq_b = ['k', 'blue', 'green', 'lightgreen', 'lightblue', 'lightgrey']
# colorbar
bounds = freqmin_b + [freqmax_b[-1]]
cmap = ListedColormap(color_freq_b)
norm = BoundaryNorm(bounds, ncolors=cmap.N, clip=True)
tick_locs = [(f0 + f1) / 2 for f0, f1 in zip(freqmin_b, freqmax_b)]
tick_labels = [f"{f0:.3f}–{f1:.3f} Hz" for f0, f1 in zip(freqmin_b, freqmax_b)]
colorbar_big = [bounds, cmap, norm, tick_locs, tick_labels]

# magnitude < 5:
freqmin_s = [0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 2, 5, 10, 20, 30, 40]
freqmax_s = [0.3, 0.5, 0.7, 0.9, 1.1, 2, 5, 10, 20, 30, 40, 49]
color_freq_s = ['grey', 'lightgrey', 'rebeccapurple', 'pink', 'darkred', 'red', 'darkorange', 'orange',
              'yellow', 'lightgreen', 'green', 'darkgreen']
# colorbar
bounds = freqmin_s + [freqmax_s[-1]]
cmap = ListedColormap(color_freq_s)
norm = BoundaryNorm(bounds, ncolors=cmap.N, clip=True)
tick_locs = [(f0 + f1) / 2 for f0, f1 in zip(freqmin_s, freqmax_s)]
tick_labels = [f"{f0:.3f}–{f1:.3f} Hz" for f0, f1 in zip(freqmin_s, freqmax_s)]
colorbar_small = [bounds, cmap, norm, tick_locs, tick_labels]

for folder in os.listdir(root_data):
    if folder in ['.DS_Store', 'Aug20.2023', 'Feb15.2025', 'Jan15.2022', 'May13.2024', 'July29.2025']:
        continue
    '''if folder not in ['Dec05.2024']:
       continue'''


    LIV, HAN = True, True

    if folder in ['July31.2024','Sept23.2024','May28.2024','July02.2024','Jan20.2025','Aug20.2023','Jan15.2022']:
        LIV = False
    if folder in ['May13.2024']:
        HAN = False
    print(folder)
    marker_event = marker_type[n_marker]

    print(n_marker)
    n_marker +=1
    h_empty = Line2D([0], [0], marker=marker_type[-1],  # e.g. 'o', '^', etc.
                     color='none', markerfacecolor=None, markeredgecolor='white', linestyle='None', markersize=4)


    if LIV:
        events = f'{path}/event_times_LIV.txt'  # Replace with your actual file path
        data = pd.read_csv(events)
        # get Backazimuth to rotate data:
        # azimuth of station:
        channel = response.select(network='LI', station='LIV', location='EX', channel='DJY')[0][0][0]
        azimuth = channel.azimuth
        # event location from event info
        data_event =  data[data['location'] == folder]
        source_latitude = data_event[' latitude'].values.astype(float)
        source_longitude = data_event[' longitude'].values.astype(float)
        starttime = UTCDateTime(str(data_event['eventstart'].iloc[0]))
        endtime = UTCDateTime(str(data_event['eventend'].iloc[0]))
        magnitude = data_event[' mag'].values.astype(float)
        station_latitude = 30.562775410180993
        station_longitude = -90.77469251757329
        dist_liv = gps2dist_azimuth(source_latitude, source_longitude, station_latitude, station_longitude)[0] / 1000
        print(magnitude)
        name = data_event[' name'].iloc[0]
        country = data_event[' country'].iloc[0]
        mags = magnitude[0]

        # theoretical backazimuth and distance
        baz = gps2dist_azimuth(source_latitude, source_longitude, station_latitude, station_longitude)
        baz_rot = (azimuth-baz[2])/180*numpy.pi
        scale_R = -1
        scale_T = 1e6

        # BRS
        LIV_EY_DJX = read('%s/%s/LI_LIV_EY_DJX.ms' %(root_data,folder))
        LIV_EX_DJY = read('%s/%s/LI_LIV_EX_DJY.ms' %(root_data,folder))
        LIV_IY_DJX = read('%s/%s/LI_LIV_IY_DJX.ms' %(root_data,folder))
        LIV_IX_DJY = read('%s/%s/LI_LIV_IX_DJY.ms' %(root_data,folder))

        # seis
        LIV_EY_DHX = read('%s/%s/LI_LIV_EY_DHX.ms' %(root_data,folder))
        LIV_EY_DHY = read('%s/%s/LI_LIV_EY_DHY.ms' %(root_data,folder))
        #LIV_EY_DHZ = read('%s/%s/LI_LIV_EY_DHZ.ms' %(root_data,folder))
        LIV_EX_DHX = read('%s/%s/LI_LIV_EX_DHX.ms' %(root_data,folder))
        LIV_EX_DHY = read('%s/%s/LI_LIV_EX_DHY.ms' %(root_data,folder))
        #LIV_EX_DHZ = read('%s/%s/LI_LIV_EX_DHZ.ms' %(root_data,folder))


        LIV_IY_DHX = read('%s/%s/LI_LIV_IY_DHX.ms' %(root_data,folder))
        LIV_IY_DHY = read('%s/%s/LI_LIV_IY_DHY.ms' %(root_data,folder))
        #LIV_IY_DHZ = read('%s/%s/LI_LIV_IY_DHZ.ms' %(root_data,folder))
        LIV_IX_DHX = read('%s/%s/LI_LIV_IX_DHX.ms' %(root_data,folder))
        LIV_IX_DHY = read('%s/%s/LI_LIV_IX_DHY.ms' %(root_data,folder))
        #LIV_IX_DHZ = read('%s/%s/LI_LIV_IX_DHZ.ms' %(root_data,folder))

        if magnitude > 5:
            freqmin = freqmin_b
            freqmax = freqmax_b
            color_freq = color_freq_b

            # colorbar
            bounds = freqmin + [freqmax[-1]]
            cmap = ListedColormap(color_freq)
            norm = BoundaryNorm(bounds, ncolors=cmap.N, clip=True)
            tick_locs = [(f0 + f1) / 2 for f0, f1 in zip(freqmin, freqmax)]
            tick_labels = [f"{f0:.3f}–{f1:.3f} Hz" for f0, f1 in zip(freqmin, freqmax)]

            colorbar_big = [bounds, cmap, norm, tick_locs, tick_labels]

            label = f'M{mags},{name},{country}'
            labels_tele.append(label)
            h = Line2D([0], [0], marker=marker_event,  # e.g. 'o', '^', etc.
                       color='none', markerfacecolor='gray', markeredgecolor='gray', linestyle='None', markersize=4)
            handles_tele.append((h, h_empty))

        elif magnitude < 5:
            freqmin = freqmin_s
            freqmax = freqmax_s
            color_freq = color_freq_s

            # colorbar
            bounds = freqmin + [freqmax[-1]]
            cmap = ListedColormap(color_freq)
            norm = BoundaryNorm(bounds, ncolors=cmap.N, clip=True)
            tick_locs = [(f0 + f1) / 2 for f0, f1 in zip(freqmin, freqmax)]
            tick_labels = [f"{f0:.3f}–{f1:.3f} Hz" for f0, f1 in zip(freqmin, freqmax)]

            colorbar_small = [bounds,cmap,norm,tick_locs,tick_labels]

            label = f'M{mags},{name},{country}'
            labels_local.append(label)
            h = Line2D([0], [0], marker=marker_event,  # e.g. 'o', '^', etc.
                       color='none', markerfacecolor='gray', markeredgecolor='gray', linestyle='None', markersize=4)
            handles_local.append((h, h_empty))



        fig_event, axs_event = plt.subplots(2, 2, figsize=(5, 5))
        plt.subplots_adjust(hspace=0, top=0.95, bottom=0.05)
        fig_event_ts_r, axs_event_ts_r = plt.subplots(2, 2, figsize=(10, 5))
        plt.subplots_adjust(hspace=0, top=0.95, bottom=0.05)
        fig_event_ts_t, axs_event_ts_t = plt.subplots(2, 2, figsize=(10, 5))
        plt.subplots_adjust(hspace=0, top=0.95, bottom=0.05)

        for freqmin, freqmax, color in zip(freqmin, freqmax, color_freq):

            for data_r, label_r, data_t, label_t, ax_event, ax_event_ts_r, ax_event_ts_t in zip([LIV_EY_DJX, LIV_EX_DJY, LIV_IY_DJX, LIV_IX_DJY],
                                                                                ['LIV_EY_DJX', 'LIV_EX_DJY', 'LIV_IY_DJX', 'LIV_IX_DJY'],
                                                                                [LIV_EY_DHY, LIV_EX_DHX, LIV_IY_DHY, LIV_IX_DHX],
                                                                                ['LIV_EY_DHY', 'LIV_EX_DHX', 'LIV_IY_DHY', 'LIV_IX_DHX'],
                                                                                axs_event.flat, axs_event_ts_r.flat, axs_event_ts_t.flat):

                data_r = data_r.filter('bandpass', freqmin=0.001, freqmax=100, zerophase=True)
                data_r = data_r.slice(starttime, endtime)
                data_r = data_r#.remove_response(inventory=response, pre_filt=None, output="DEF", water_level=60, plot=False)
                data_r_proc = data_r.copy().detrend('simple').taper(0.1).filter('bandpass', freqmin=freqmin, freqmax=freqmax,zerophase=True)

                data_t = data_t.filter('bandpass', freqmin=0.001, freqmax=100, zerophase=True)
                data_t = data_t.slice(starttime, endtime)
                data_t = data_t.remove_response(inventory=response, pre_filt=None, output="DEF", water_level=60, plot=False)
                data_t_proc = data_t.copy().detrend('simple').taper(0.1).filter('bandpass', freqmin=freqmin, freqmax=freqmax, zerophase=True)

                rot_error_acc_t = data_t_proc.copy()
                angle_rad = data_r_proc[0].data/1e9
                rot_error_acc_t[0].data = sin(angle_rad)*9.81
                rot_error_vel_t = rot_error_acc_t.copy().filter('bandpass', freqmin=freqmin, freqmax=freqmax, zerophase=True).integrate().filter('bandpass', freqmin=freqmin, freqmax=freqmax, zerophase=True)
                PREV = max(abs(rot_error_vel_t[0].data))
                PGV = max(abs(data_t_proc[0].data))
                PGRA = max(abs(angle_rad))

                if magnitude > 5:
                    axs_tele_re.scatter(PGV, PREV, marker=marker_event, color=color, label=f'f {freqmin}-{freqmax} Hz')
                    axs_tele_am.scatter(PGV, PGRA, marker=marker_event, color=color, label=f'f {freqmin}-{freqmax} Hz')
                elif magnitude < 5:
                    axs_local_re.scatter(PGV, PREV, marker=marker_event, color=color, label=f'f {freqmin}-{freqmax} Hz')
                    axs_local_am.scatter(PGV, PGRA, marker=marker_event, color=color, label=f'f {freqmin}-{freqmax} Hz')

                ax_event.scatter(PGV, PREV, marker=marker_event, color=color, label=f'f {freqmin}-{freqmax} Hz')
                ax_event_ts_r.plot(rot_error_vel_t[0].times(),rot_error_vel_t[0].data, linewidth=0.2, color=color, label=f'f {freqmin}-{freqmax} Hz')
                ax_event_ts_t.plot(data_t_proc[0].times(),data_t_proc[0].data, linewidth=0.2, color=color, label=f'f {freqmin}-{freqmax} Hz')

                ax_event.set_xscale("log")
                ax_event.set_yscale("log")
        for ax in [axs_event.flat[0], axs_event_ts_r.flat[0], axs_event_ts_t.flat[0]]:
            ax.legend(loc='upper left')
        fig_event.savefig('%s/LI_LIV_%s_rot_error.png' %(root_savefig,folder), dpi=300,bbox_inches='tight')
        fig_event_ts_r.savefig('%s/LI_LIV_%s_rot_error_ts_r.png' %(root_savefig,folder), dpi=300,bbox_inches='tight')
        fig_event_ts_t.savefig('%s/LI_LIV_%s_rot_error_ts_t.png' %(root_savefig,folder), dpi=300,bbox_inches='tight')
        #plt.show()
        #plt.close()
        #exit()

    if HAN:
        events = f'{path}/event_times_HAN.txt'  # Replace with your actual file path
        data = pd.read_csv(events)
        # get Backazimuth to rotate data:
        # azimuth of station:
        channel = response.select(network='LI', station='HAN', location='EX', channel='DJY')[0][0][0]
        azimuth = channel.azimuth
        # event location from event info
        data_event =  data[data['location'] == folder]
        source_latitude = data_event[' latitude'].values.astype(float)
        source_longitude = data_event[' longitude'].values.astype(float)
        starttime = UTCDateTime(str(data_event['eventstart'].iloc[0]))
        endtime = UTCDateTime(str(data_event['eventend'].iloc[0]))
        magnitude = data_event[' mag'].values.astype(float)

        name = data_event[' name'].iloc[0]
        country = data_event[' country'].iloc[0]
        mags = magnitude[0]


        # theoretical backazimuth and distance
        station_latitude = 46.455 #47.5
        station_longitude = -119.407 #-120.7
        dist_han = gps2dist_azimuth(source_latitude, source_longitude, station_latitude, station_longitude)[0] / 1000
        baz = gps2dist_azimuth(source_latitude, source_longitude, station_latitude, station_longitude)
        baz_rot = (azimuth - baz[2]) / 180 * numpy.pi
        scale_R = -1
        scale_T = 1e6

        # BRS
        HAN_EY_DJX = read('%s/%s/LI_HAN_EY_DJX.ms' % (root_data, folder))
        HAN_EX_DJY = read('%s/%s/LI_HAN_EX_DJY.ms' % (root_data, folder))

        # seis
        HAN_EY_DHX = read('%s/%s/LI_HAN_EY_DHX.ms' % (root_data, folder))
        HAN_EY_DHY = read('%s/%s/LI_HAN_EY_DHY.ms' % (root_data, folder))
        #HAN_EY_DHZ = read('%s/%s/LI_HAN_EY_DHZ.ms' % (root_data, folder))
        HAN_EX_DHX = read('%s/%s/LI_HAN_EX_DHX.ms' % (root_data, folder))
        HAN_EX_DHY = read('%s/%s/LI_HAN_EX_DHY.ms' % (root_data, folder))
        #HAN_EX_DHZ = read('%s/%s/LI_HAN_EX_DHZ.ms' % (root_data, folder))


        if magnitude > 5:
            freqmin = freqmin_b
            freqmax = freqmax_b
            color_freq = color_freq_b

            # colorbar
            bounds = freqmin + [freqmax[-1]]
            cmap = ListedColormap(color_freq)
            norm = BoundaryNorm(bounds, ncolors=cmap.N, clip=True)
            tick_locs = [(f0 + f1) / 2 for f0, f1 in zip(freqmin, freqmax)]
            tick_labels = [f"{f0:.3f}–{f1:.3f} Hz" for f0, f1 in zip(freqmin, freqmax)]

            colorbar_big = [bounds,cmap,norm,tick_locs,tick_labels]

            label = f'M{mags},{name},{country}'
            if label not in labels_tele:
                labels_tele.append(label)
            h = Line2D([0], [0], marker=marker_event,  # e.g. 'o', '^', etc.
                       color='none', markerfacecolor='None', markeredgecolor='gray', linestyle='None', markersize=4)
            if LIV:
                handles_tele[-1] = (handles_tele[-1][0], h)
            else:
                handles_tele.append((h_empty, h))


        elif magnitude < 5:
            freqmin = freqmin_s
            freqmax = freqmax_s
            color_freq = color_freq_s

            # colorbar
            bounds = freqmin + [freqmax[-1]]
            cmap = ListedColormap(color_freq)
            norm = BoundaryNorm(bounds, ncolors=cmap.N, clip=True)
            tick_locs = [(f0 + f1) / 2 for f0, f1 in zip(freqmin, freqmax)]
            tick_labels = [f"{f0:.3f}–{f1:.3f} Hz" for f0, f1 in zip(freqmin, freqmax)]

            colorbar_small = [bounds,cmap,norm,tick_locs,tick_labels]

            label = f'M{mags},{name},{country}'
            if label not in labels_local:
                labels_local.append(label)
            h = Line2D([0], [0], marker=marker_event,  # e.g. 'o', '^', etc.
                       color='none', markerfacecolor='None', markeredgecolor='gray', linestyle='None', markersize=4)
            if LIV:
                handles_local[-1] = (handles_local[-1][0], h)
            else:
                handles_local.append((h_empty, h))



        fig_event, axs_event = plt.subplots(1, 2, figsize=(5, 3))
        plt.subplots_adjust(hspace=0, top=0.95, bottom=0.05)
        fig_event_ts_r, axs_event_ts_r = plt.subplots(1, 2, figsize=(10, 3))
        plt.subplots_adjust(hspace=0, top=0.95, bottom=0.05)
        fig_event_ts_t, axs_event_ts_t = plt.subplots(1, 2, figsize=(10, 3))
        plt.subplots_adjust(hspace=0, top=0.95, bottom=0.05)


        for freqmin, freqmax, color in zip(freqmin, freqmax, color_freq):
            for data_r, label_r, data_t, label_t, ax_event,ax_event_ts_r, ax_event_ts_t in zip([HAN_EX_DJY, HAN_EY_DJX],
                                                                  ['HAN_EX_DJY', 'HAN_EY_DJX'],
                                                                  [HAN_EX_DHX, HAN_EY_DHY],
                                                                  ['HAN_EX_DHX', 'HAN_EY_DHY'],
                                                                  axs_event.flat,axs_event_ts_r.flat,axs_event_ts_t.flat):
                data_r = data_r.filter('bandpass', freqmin=0.001, freqmax=100, zerophase=True)
                data_r = data_r.slice(starttime, endtime)
                data_r = data_r  # .remove_response(inventory=response, pre_filt=None, output="DEF", water_level=60, plot=False)
                data_r_proc = data_r.copy().detrend('simple').taper(0.1).filter('bandpass', freqmin=freqmin,
                                                                                freqmax=freqmax, zerophase=True)

                data_t = data_t.filter('bandpass', freqmin=0.001, freqmax=100, zerophase=True)
                data_t = data_t.slice(starttime, endtime)
                data_t = data_t.remove_response(inventory=response, pre_filt=None, output="DEF", water_level=60,
                                                plot=False)
                data_t_proc = data_t.copy().detrend('simple').taper(0.1).filter('bandpass', freqmin=freqmin,
                                                                                freqmax=freqmax, zerophase=True)

                rot_error_acc_t = data_t_proc.copy()
                rot_error_acc_t[0].stats.sampling_rate = data_r_proc[0].stats.sampling_rate
                angle_rad = data_r_proc[0].data / 1e9
                rot_error_acc_t[0].data = sin(angle_rad) * 9.81
                rot_error_vel_t = rot_error_acc_t.copy().filter('bandpass', freqmin=freqmin, freqmax=freqmax,
                                                                zerophase=True).integrate().filter('bandpass',
                                                                                                   freqmin=freqmin,
                                                                                                   freqmax=freqmax,
                                                                                                   zerophase=True)
                PREV = max(abs(rot_error_vel_t[0].data))
                PGV = max(abs(data_t_proc[0].data))
                PGRA = max(abs(angle_rad))

                if magnitude > 5:
                    axs_tele_re.scatter(PGV, PREV, marker=marker_event, facecolor = 'None', edgecolor=color, label=f'f {freqmin}-{freqmax} Hz')
                    axs_tele_am.scatter(PGV, PGRA, marker=marker_event, facecolor = 'None', edgecolor=color, label=f'f {freqmin}-{freqmax} Hz')

                elif magnitude < 5:
                    axs_local_re.scatter(PGV, PREV, marker=marker_event, facecolor = 'None', edgecolor=color, label=f'f {freqmin}-{freqmax} Hz')
                    axs_local_am.scatter(PGV, PGRA, marker=marker_event, facecolor = 'None', edgecolor=color, label=f'f {freqmin}-{freqmax} Hz')

                ax_event.scatter(PGV, PREV, marker=marker_event, facecolor = 'None', edgecolor=color, label=f'f {freqmin}-{freqmax} Hz')
                ax_event_ts_r.plot(rot_error_vel_t[0].times(),rot_error_vel_t[0].data, linewidth=0.2, color=color, label=f'f {freqmin}-{freqmax} Hz')
                ax_event_ts_t.plot(data_t_proc[0].times(),data_t_proc[0].data, linewidth=0.2, color=color, label=f'f {freqmin}-{freqmax} Hz')

                ax_event.set_xscale("log")
                ax_event.set_yscale("log")
        for ax in [axs_event.flat[0], axs_event_ts_r.flat[0], axs_event_ts_t.flat[0]]:
            ax.legend(loc='upper left')
        fig_event.savefig('%s/LI_HAN_%s_rot_error.png' % (root_savefig, folder), dpi=300, bbox_inches='tight')
        fig_event_ts_r.savefig('%s/LI_HAN_%s_rot_error_ts_r.png' %(root_savefig,folder), dpi=300,bbox_inches='tight')
        fig_event_ts_t.savefig('%s/LI_HAN_%s_rot_error_ts_t.png' %(root_savefig,folder), dpi=300,bbox_inches='tight')

toc = time.time()
print((toc-tic)/60)



handles_tele = [x for x in handles_tele if not (x[0] is h_empty and x[1] is h_empty)]
handles_local = [x for x in handles_local if not (x[0] is h_empty and x[1] is h_empty)]

# HAN
h = Line2D([0], [0],marker=marker_type[-1],     # e.g. 'o', '^', etc.
    color='none', markerfacecolor=None,markeredgecolor='gray',linestyle='None',markersize=4)
handles_tele.append(h)
handles_local.append(h)
labels_tele.append('HAN Station')
labels_local.append('HAN Station')
# LIV
h = Line2D([0], [0],marker=marker_type[-1],     # e.g. 'o', '^', etc.
    color='none', markerfacecolor='gray',markeredgecolor='gray',linestyle='None',markersize=4)
handles_tele.append(h)
handles_local.append(h)
labels_tele.append('LIV Station')
labels_local.append('LIV Station')


for fig, ax, colorbar, title, ylabel, handles, labels in zip([fig_tele_re,fig_local_re, fig_tele_am, fig_local_am],
                                    [axs_tele_re, axs_local_re, axs_tele_am, axs_local_am],
                                    [colorbar_big,colorbar_small,colorbar_big,colorbar_small],
                                    ['Teleseismic events with Mw > 5', 'Local events with M < 5','Teleseismic events with Mw > 5', 'Local events with M < 5'],
                                    ['Rotation error PREV [m/s]','Rotation error PREV [m/s]','Rotation angle PGRA [rad]','Rotation angle PGRA [rad]'],
                                    [handles_tele, handles_local, handles_tele, handles_local],
                                    [labels_tele, labels_local, labels_tele, labels_local]
                                                             ):
    ax.set_title(title)
    ax.set_xscale("log")
    ax.set_yscale("log")
    #ax.set_ylimit(bottom=10^-12)
    #ax.set_xlimit(left=10^-11)
    if 'PREV' in ylabel:
        xx = numpy.logspace(numpy.log10(1e-6), numpy.log10(1e-3), 100)
        ax.plot(xx, xx, color='gray', linestyle='-',linewidth=1, zorder=0)
        ax.text(xx[-1],xx[-1],'100%', ha='right')
        ax.plot(xx, xx/10, color='gray', linestyle='--',linewidth=1, zorder=0)
        ax.text(xx[-1],xx[-1]/10,'10%', ha='right')
        ax.plot(xx, xx/100, color='gray', linestyle='-.',linewidth=1, zorder=0)
        ax.text(xx[-1],xx[-1]/100,'1%', ha='right')
    else:
        xx = numpy.logspace(numpy.log10(1e-6), numpy.log10(1e-3), 100)
        ax.plot(xx, xx, color='gray', linestyle='-',linewidth=1, zorder=0)
        ax.text(xx[-1],xx[-1],'1', ha='right')
        ax.plot(xx, xx/100, color='gray', linestyle='--',linewidth=1, zorder=0)
        ax.text(xx[-1],xx[-1]/100,'0.01', ha='right')
        ax.plot(xx, xx/10000, color='gray', linestyle='-.',linewidth=1, zorder=0)
        ax.text(xx[-1],xx[-1]/10000,'0.001', ha='right')


    ax.legend(handles,labels, handler_map={tuple: HandlerTuple(ndivide=None)}, loc='upper left',frameon=False, fontsize='small')
    ax.set_xlabel('PGV [m/s]')
    ax.set_ylabel(ylabel)

    pos = ax.get_position()

    # create a new Axes to the right: [x0, y0, width, height]
    cax = fig.add_axes([
        pos.x1 + 0.02,  # a little to the right of your main axes
        pos.y0,         # same bottom
        0.02,           # thin colorbar
        pos.height      # same height
    ])

    cb = ColorbarBase(
        cax,
        cmap=colorbar[1],#cmap,
        norm=colorbar[2],#norm,
        boundaries=colorbar[0], #bounds,
        ticks=colorbar[3],#norm,
        spacing='uniform',
        orientation='vertical'
    )
    cb.set_ticklabels(colorbar[4])#tick_labels)
    cb.set_label('Frequency band [Hz]', labelpad=10)

fig_tele_re.savefig('%s/LI_teleseismic_rot_error.png' %(root_savefig), dpi=300,bbox_inches='tight')
fig_local_re.savefig('%s/LI_local_rot_error.png' %(root_savefig), dpi=300,bbox_inches='tight')


fig_tele_am.savefig('%s/LI_teleseismic_rot_angle.png' %(root_savefig), dpi=300,bbox_inches='tight')
fig_local_am.savefig('%s/LI_local_rot_angle.png' %(root_savefig), dpi=300,bbox_inches='tight')