import os
import pandas as pd
import numpy
import matplotlib.pyplot as plt
from numpy import sin, cos, pi
from obspy import read, UTCDateTime, read_inventory
from correlation_function import correlate2
from obspy.geodetics.base import gps2dist_azimuth
from functions import rot_vec, supylabel2
from datetime import datetime

root_data = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data/miniseed_LIGO'
root_savefig = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/2_Figures/RotationError'
root_xml = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data/metadata'
response = read_inventory('%s/LIGO_HAN_LIV.xml' %root_xml)
path = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data'

fig_tele, axs_tele = plt.subplots(1, 1, figsize=(5, 5))
plt.subplots_adjust(hspace=0, top=0.95, bottom=0.05)
fig_local, axs_local = plt.subplots(1, 1, figsize=(5, 5))
plt.subplots_adjust(hspace=0, top=0.95, bottom=0.05)

for folder in os.listdir(root_data):
    if folder in ['.DS_Store', 'Aug20.2023', 'Feb15.2025']:
        continue
    '''if 'May28.2024' not in folder:
       continue'''


    LIV, HAN = True, True

    if folder in ['July31.2024','Sept23.2024','May28.2024','July02.2024','Jan20.2025','Aug20.2023','Jan15.2022']:
        LIV = False
    if folder in ['May13.2024']:
        HAN = False
    print(folder)



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
            freqmin = [0.008,0.02, 0.04, 0.06, 0.08, 0.1]
            freqmax = [0.02, 0.04, 0.06, 0.08, 0.1,  49]
            color_freq = ['k', 'blue', 'green', 'lightblue', 'lightgreen']
            #freqmin = [0.008,0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.008]
            #freqmax = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1,  49]
        elif magnitude < 5:
            freqmin = [0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 2,  5, 10, 20,30,40]
            freqmax = [0.3, 0.5, 0.7, 0.9, 1.1,   2, 5, 10, 20, 30,40,49]
            color_freq = ['grey', 'lightgrey', 'rebeccapurple', 'pink', 'darkred', 'red', 'darkorange', 'orange', 'yellow', 'lightgreen', 'green', 'darkgreen']

        fig_event, axs_event = plt.subplots(2, 2, figsize=(5, 5))
        plt.subplots_adjust(hspace=0, top=0.95, bottom=0.05)
        for freqmin, freqmax, color in zip(freqmin, freqmax, color_freq):

            for data_r, label_r, data_t, label_t, ax_event in zip([LIV_EY_DJX, LIV_EX_DJY, LIV_IY_DJX, LIV_IX_DJY],
                                                                    ['LIV_EY_DJX', 'LIV_EX_DJY', 'LIV_IY_DJX', 'LIV_IX_DJY'],
                                                                    [LIV_EY_DHY, LIV_EX_DHX, LIV_IY_DHY, LIV_IX_DHX],
                                                                    ['LIV_EY_DHY', 'LIV_EX_DHX', 'LIV_IY_DHY', 'LIV_IX_DHX'],
                                                                   axs_event.flat):

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

                if magnitude > 5:
                    axs_tele.scatter(PGV, PREV, marker='o', color=color, label=f'f {freqmin}-{freqmax} Hz')
                elif magnitude < 5:
                    axs_local.scatter(PGV, PREV, marker='o', color=color, label=f'f {freqmin}-{freqmax} Hz')

                ax_event.scatter(PGV, PREV, marker='o', color=color, label=f'f {freqmin}-{freqmax} Hz')

                ax_event.set_xscale("log")
                ax_event.set_yscale("log")
                fig_event.savefig('%s/LI_LIV_%s_rot_error.png' %(root_savefig,folder), dpi=300,bbox_inches='tight')
        #plt.show()
        #plt.close()
        #exit()

    if HAN:
        continue
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

        HAN_E_DJ = numpy.array([HAN_EY_DJX[0].data, HAN_EX_DJY[0].data, HAN_EX_DJY[0].data*0])
        HAN_E_DJ_baz = rot_vec(0,0,baz_rot, HAN_E_DJ)
        HAN_E_T = HAN_EY_DJX.copy()
        HAN_E_R = HAN_EX_DJY.copy()
        HAN_E_T[0].data = HAN_E_DJ_baz[0,:]
        HAN_E_R[0].data = HAN_E_DJ_baz[1,:]

        # seis
        HAN_EY_DHX = read('%s/%s/LI_HAN_EY_DHX.ms' % (root_data, folder))
        HAN_EY_DHY = read('%s/%s/LI_HAN_EY_DHY.ms' % (root_data, folder))
        HAN_EY_DHZ = read('%s/%s/LI_HAN_EY_DHZ.ms' % (root_data, folder))
        HAN_EX_DHX = read('%s/%s/LI_HAN_EX_DHX.ms' % (root_data, folder))
        HAN_EX_DHY = read('%s/%s/LI_HAN_EX_DHY.ms' % (root_data, folder))
        HAN_EX_DHZ = read('%s/%s/LI_HAN_EX_DHZ.ms' % (root_data, folder))


        if magnitude > 5:
            freqmin = [0.008,0.02, 0.04, 0.06, 0.08, 0.008]
            freqmax = [0.02, 0.04, 0.06, 0.08, 0.1,  49]
            #freqmin = [0.008, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.008]
            #freqmax = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 49]
        elif magnitude < 2:
            freqmin = [0.3, 0.8]
            freqmax = [0.8, 30]
        elif magnitude < 5:
            freqmin = [0.1, 0.8, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10, 20, 30, 0.008]
            freqmax = [0.8, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10, 20, 30, 49, 49]

        fig_Ey, axs_Ey = plt.subplots(len(freqmin), 1, figsize=(8.5, len(freqmin)*0.7), sharex=True)
        plt.subplots_adjust(hspace=0, top=0.95, bottom=0.05)
        fig_Ex, axs_Ex = plt.subplots(len(freqmin), 1, figsize=(8.5, len(freqmin)*0.7), sharex=True)
        plt.subplots_adjust(hspace=0, top=0.95, bottom=0.05)
        fig_ET, axs_ET = plt.subplots(len(freqmin), 1, figsize=(8.5, len(freqmin) * 0.7), sharex=True)
        plt.subplots_adjust(hspace=0, top=0.95, bottom=0.05)
        fig_ER, axs_ER = plt.subplots(len(freqmin), 1, figsize=(8.5, len(freqmin) * 0.7), sharex=True)
        plt.subplots_adjust(hspace=0, top=0.95, bottom=0.05)

        '''for axs in [axs_Ey, axs_Ex]:
            for ax in axs:
                ax.tick_params(labelleft=False, labelright=False)
                ax.tick_params(left=False, right=False)'''

        for freqmin, freqmax, ax_ex, ax_ey, ax_et, ax_er in zip(freqmin, freqmax, axs_Ex.flat, axs_Ey.flat, axs_ET.flat, axs_ER.flat):

            for data_t, label_t, ax in zip([HAN_EX_DHZ, HAN_EY_DHZ, HAN_EX_DHZ, HAN_EY_DHZ],
                                           ['HAN_EX_DHZ', 'HAN_EY_DHZ', 'HAN_EX_DHZ', 'HAN_EY_DHZ'],
                                           [ax_ex, ax_ey, ax_et, ax_er]):
                data_t = data_t.filter('bandpass', freqmin=0.008, freqmax=100, zerophase=True).resample(256)
                data_t = data_t.slice(starttime, endtime)
                data_t = data_t.remove_response(inventory=response, pre_filt=None, output="DEF", water_level=60,plot=False)
                data_proc = data_t.copy().detrend('simple').taper(0.1).filter('bandpass', freqmin=freqmin,freqmax=freqmax, zerophase=True)
                ax.plot(data_proc[0].times(), data_proc[0].data*scale_T, color='k', label=label_t)

            for data_r, label_r, ax in zip([HAN_EX_DJY, HAN_EY_DJX, HAN_E_T, HAN_E_R],
                                           ['HAN_EX_DJY', 'HAN_EY_DJX', 'HAN_E_T', 'HAN_E_R'],
                                           [ax_ex, ax_ey, ax_et, ax_er]):
                ax_twin = ax.twinx()
                data_r = data_r.filter('bandpass', freqmin=0.001, freqmax=100, zerophase=True)
                data_r = data_r.slice(starttime, endtime)
                data_r = data_r.remove_response(inventory=response, pre_filt=None, output="DEF", water_level=60,plot=False)
                data_proc = data_r.copy().detrend('simple').taper(0.1).filter('bandpass', freqmin=freqmin,freqmax=freqmax, zerophase=True)
                ax_twin.plot(data_proc[0].times(), data_proc[0].data*scale_R, color='red', label=label_r)
                ax_twin.yaxis.label.set_color('red')
                ax_twin.tick_params(axis='y', colors='red')

                # ax.legend(loc='upper right')
                ax.margins(0)
                #ax.xaxis_date()
                ax.text(0.80, 0.1, f'{freqmin}-{freqmax} Hz', transform=ax_twin.transAxes)
            start = data_proc[0].stats.starttime
        formatted_time_new = start.strftime('%Y-%m-%dT%H:%M:%S')
        for ax, title, scale in zip([axs_Ex.flat, axs_Ey.flat, axs_ET.flat, axs_ER.flat],
                                     ['EX', 'EY', 'ET', 'ER'],
                                     [scale_R,scale_R,scale_R,scale_R]):
            #ax[0].set_title(f'HAN_{folder}_{title}, M{round(float(magnitude[0]),1)}, baz = {round(baz[2],2)}°, $vel_Z$ vs $rot_H$*{scale}')
            ax[0].set_title(f'HAN @ {title} | {formatted_time_new} | M{round(float(magnitude[0]),1)} | dist={int(dist_han)} km | BAz = {round(baz[2],2)}° | $vel_Z$ vs $rot_H \cdot${scale}')
            ax[-1].set_xlabel(f'Time [s]')
        for fig, name in zip([fig_Ey, fig_Ex, fig_ET, fig_ER], ['EY', 'EX', 'ET', 'ER']):
            # fig.savefig(f'{root_savefig}/{folder}/LI_HAN_{folder}_{name}_dispersion_all_stations.png', dpi=150, bbox_inches='tight')
            fig.supylabel(f'velocity [$\mu$m/s]', ha='right', x=0.07)
            supylabel2(fig, f'angle [nrad]', color='red')
            fig.savefig('%s/LI_HAN_%s_%s_dispersion_all_stations.png' % (root_savefig, folder, name), dpi=300,bbox_inches='tight')
        #plt.show()
        # exit()
for ax in [axs_tele, axs_local]:
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.axline((0, 0), (1e-5, 1e-5),transform=ax.transAxes,
              color='gray', linestyle='-',linewidth=1, zorder=0)
    ax.axline((0, 0), (1e-6, 1e-5),transform=ax.transAxes,
              color='gray', linestyle='--',linewidth=1, zorder=0)
    ax.axline((0, 0), (1e-7, 1e-5),transform=ax.transAxes,
              color='gray', linestyle='-.',linewidth=1, zorder=0)
fig_tele.savefig('%s/LI_teleseismic_rot_error.png' %(root_savefig), dpi=300,bbox_inches='tight')
fig_local.savefig('%s/LI_local_rot_error.png' %(root_savefig), dpi=300,bbox_inches='tight')