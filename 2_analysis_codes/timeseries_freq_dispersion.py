import os
import pandas as pd
import numpy
import matplotlib.pyplot as plt
from obspy import read, UTCDateTime, read_inventory
from correlation_function import correlate2
from obspy.geodetics.base import gps2dist_azimuth
from functions import rot_vec, supylabel2
from datetime import datetime

root_data = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data/miniseed_LIGO'
root_savefig = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/2_Figures/Timeseries_freq_dispersion'
root_xml = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data/metadata'
response = read_inventory('%s/LIGO_HAN_LIV.xml' %root_xml)
path = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data'



for folder in os.listdir(root_data):
    if folder in ['.DS_Store', 'Aug20.2023', 'Feb15.2025']:
        continue
    '''if 'Apr02.2024' not in folder:
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
        print(magnitude)

        station_latitude = 30.562775410180993
        station_longitude = -90.77469251757329

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

        LIV_E_DJ = numpy.array([LIV_EY_DJX[0].data, LIV_EX_DJY[0].data, LIV_EX_DJY[0].data*0])
        LIV_E_DJ_baz = rot_vec(0,0,baz_rot, LIV_E_DJ)
        LIV_E_T = LIV_EY_DJX.copy()
        LIV_E_R = LIV_EX_DJY.copy()
        LIV_E_T[0].data = LIV_E_DJ_baz[0,:]
        LIV_E_R[0].data = LIV_E_DJ_baz[1,:]

        LIV_I_DJ = numpy.array([LIV_IY_DJX[0].data, LIV_IX_DJY[0].data, LIV_IX_DJY[0].data*0])
        LIV_I_DJ_baz = rot_vec(0,0,baz_rot, LIV_I_DJ)
        LIV_I_T = LIV_IY_DJX.copy()
        LIV_I_R = LIV_IX_DJY.copy()
        LIV_I_T[0].data = LIV_I_DJ_baz[0,:]
        LIV_I_R[0].data = LIV_I_DJ_baz[1,:]

        # seis
        LIV_EY_DHX = read('%s/%s/LI_LIV_EY_DHX.ms' %(root_data,folder))
        LIV_EY_DHY = read('%s/%s/LI_LIV_EY_DHY.ms' %(root_data,folder))
        LIV_EY_DHZ = read('%s/%s/LI_LIV_EY_DHZ.ms' %(root_data,folder))
        LIV_EX_DHX = read('%s/%s/LI_LIV_EX_DHX.ms' %(root_data,folder))
        LIV_EX_DHY = read('%s/%s/LI_LIV_EX_DHY.ms' %(root_data,folder))
        LIV_EX_DHZ = read('%s/%s/LI_LIV_EX_DHZ.ms' %(root_data,folder))


        LIV_IY_DHX = read('%s/%s/LI_LIV_IY_DHX.ms' %(root_data,folder))
        LIV_IY_DHY = read('%s/%s/LI_LIV_IY_DHY.ms' %(root_data,folder))
        LIV_IY_DHZ = read('%s/%s/LI_LIV_IY_DHZ.ms' %(root_data,folder))
        LIV_IX_DHX = read('%s/%s/LI_LIV_IX_DHX.ms' %(root_data,folder))
        LIV_IX_DHY = read('%s/%s/LI_LIV_IX_DHY.ms' %(root_data,folder))
        LIV_IX_DHZ = read('%s/%s/LI_LIV_IX_DHZ.ms' %(root_data,folder))

        if magnitude > 5:
            #freqmin = [0.008,0.02, 0.04, 0.06, 0.08, 0.008]
            #freqmax = [0.02, 0.04, 0.06, 0.08, 0.1,  49]
            freqmin = [0.008,0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.008]
            freqmax = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1,  49]
        elif magnitude < 5:
            freqmin = [0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 2,  5, 15, 25,35,0.008]
            freqmax = [0.3, 0.5, 0.7, 0.9, 1.1,   2, 5, 10, 20, 30,40,49]

        fig_Ey, axs_Ey = plt.subplots(len(freqmin), 1, figsize=(11,8), sharex=True)
        plt.subplots_adjust(hspace=0, top=0.95, bottom=0.05)
        fig_Ex, axs_Ex = plt.subplots(len(freqmin), 1, figsize=(11,8), sharex=True)
        plt.subplots_adjust(hspace=0, top=0.95, bottom=0.05)
        fig_Iy, axs_Iy = plt.subplots(len(freqmin), 1, figsize=(11,8), sharex=True)
        plt.subplots_adjust(hspace=0, top=0.95, bottom=0.05)
        fig_Ix, axs_Ix = plt.subplots(len(freqmin), 1, figsize=(11,8), sharex=True)
        plt.subplots_adjust(hspace=0, top=0.95, bottom=0.05)

        '''for axs in [axs_Ey, axs_Ex, axs_Iy, axs_Ix]:
            for ax in axs:
                ax.tick_params(labelleft=False, labelright=False)
                ax.tick_params(left=False, right=False)'''



        for freqmin, freqmax, ax_ix, ax_iy, ax_ex, ax_ey in zip(freqmin, freqmax, axs_Ix.flat, axs_Iy.flat, axs_Ex.flat, axs_Ey.flat):

            for data_t, label_t, ax in zip([LIV_IX_DHZ, LIV_IY_DHZ, LIV_EX_DHZ, LIV_EY_DHZ],
                                           ['LIV_IX_DHZ', 'LIV_IY_DHZ', 'LIV_EX_DHZ', 'LIV_EY_DHZ'],
                                           [ax_ix, ax_iy, ax_ex, ax_ey]):
                data_t = data_t.filter('bandpass', freqmin=0.008, freqmax=100, zerophase=True).resample(256)
                data_t = data_t.slice(starttime, endtime)
                data_t = data_t.remove_response(inventory=response, pre_filt=None, output="VEL", water_level=60, plot=False)
                data_proc = data_t.copy().detrend('simple').taper(0.1).filter('bandpass', freqmin=freqmin, freqmax=freqmax,zerophase=True)
                ax.plot(data_proc[0].times(), data_proc[0].data*scale_T, color='k', label=label_t)

            for data_r, label_r, ax in zip([ LIV_I_T, LIV_I_R, LIV_E_T, LIV_E_R],
                                       ['LIV_I_T','LIV_I_R','LIV_E_T','LIV_E_R'],
                                       [ax_ix, ax_iy, ax_ex, ax_ey]):
                ax_twin = ax.twinx()
                data_r = data_r.filter('bandpass', freqmin=0.001, freqmax=100, zerophase=True)
                data_r = data_r.slice(starttime, endtime)
                data_r = data_r#.remove_response(inventory=response, pre_filt=None, output="DEF", water_level=60, plot=False)
                data_proc = data_r.copy().detrend('simple').taper(0.1).filter('bandpass', freqmin=freqmin, freqmax=freqmax,zerophase=True)
                ax_twin.plot(data_proc[0].times(), data_proc[0].data*scale_R, color='red', label=label_r)
                ax_twin.yaxis.label.set_color('red')
                ax_twin.tick_params(axis='y', colors='red')

                #ax.legend(loc='upper right')
                ax.margins(0)
                #ax.xaxis_date()
                ax.text(0.80,0.1, f'{freqmin}-{freqmax} Hz', transform=ax_twin.transAxes)
        start = data_proc[0].stats.starttime
        formatted_time_new = start.strftime('%Y-%m-%dT%H:%M:%S')
        for ax, title, scale in zip([axs_Ix.flat, axs_Iy.flat, axs_Ex.flat, axs_Ey.flat],
                                    ['IT', 'IR', 'ET', 'ER'],
                                    [scale_R,scale_R,scale_R,scale_R]):
            #ax[0].set_title(f'LIV_{folder}_{title}, M{round(float(magnitude[0]),1)}, baz = {round(baz[2],2)}°, acc_Z vs rot_H*{scale}')
            ax[0].set_title(f'LIV @ {title} | {formatted_time_new} | M{round(float(magnitude[0]),1)} | BAz = {round(baz[2],2)}° | acc_Z vs rot_H*{scale}')

            ax[-1].set_xlabel(f'Time [s]')
        for fig, name in zip([fig_Ey, fig_Ex, fig_Iy, fig_Ix],['ER', 'ET', 'IR', 'IT']):
            #fig.savefig(f'{root_savefig}/{folder}/LI_LIV_{folder}_{name}_dispersion_all_stations.png', dpi=150, bbox_inches='tight')
            fig.supylabel(f'Velocity [$\mu$m]')
            supylabel2(fig, f'Angle [nrad]', color='red')
            fig.savefig('%s/LI_LIV_%s_%s_dispersion_all_stations.png' %(root_savefig,folder,name))#, dpi=10,bbox_inches='tight')
        #plt.show()
        #plt.close()
        #exit()

    if HAN:
        events = f'{path}/event_times_HAN.txt'  # Replace with your actual file path
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

        station_latitude = 47.5
        station_longitude = -120.7

        # theoretical backazimuth and distance
        # here I don't rotate, as they are too far away for the events.
        baz = gps2dist_azimuth(source_latitude, source_longitude, station_latitude, station_longitude)
        baz_rot = 0 #(azimuth - baz[2]) / 180 * numpy.pi
        scale_R = -1

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


        fig_Ey, axs_Ey = plt.subplots(12, 1, figsize=(11, 8), sharex=True)
        plt.subplots_adjust(hspace=0, top=0.95, bottom=0.05)
        fig_Ex, axs_Ex = plt.subplots(12, 1, figsize=(11, 8), sharex=True)
        plt.subplots_adjust(hspace=0, top=0.95, bottom=0.05)

        for axs in [axs_Ey, axs_Ex]:
            for ax in axs:
                ax.tick_params(labelleft=False, labelright=False)
                ax.tick_params(left=False, right=False)

        if magnitude > 5:
            freqmin = [0.008, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.008]
            freqmax = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 49]
        elif magnitude < 5:
            freqmin = [0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 2, 5, 15, 25, 35, 0.008]
            freqmax = [0.3, 0.5, 0.7, 0.9, 1.1, 2, 5, 10, 20, 30, 40, 49]

        for freqmin, freqmax, ax_ex, ax_ey in zip(freqmin, freqmax, axs_Ex.flat, axs_Ey.flat):

            for data_t, label_t, ax in zip([HAN_EX_DHZ, HAN_EY_DHZ],
                                           ['HAN_EX_DHZ', 'HAN_EY_DHZ'],
                                           [ax_ex, ax_ey]):
                data_t = data_t.filter('bandpass', freqmin=0.008, freqmax=100, zerophase=True).resample(256)
                data_t = data_t.slice(starttime, endtime)
                data_t = data_t.remove_response(inventory=response, pre_filt=None, output="DEF", water_level=60,plot=False)
                data_proc = data_t.copy().detrend('simple').taper(0.1).filter('bandpass', freqmin=freqmin,freqmax=freqmax, zerophase=True)
                ax.plot(data_proc[0].times(), data_proc[0].data, color='k', label=label_t)

            for data_r, label_r, ax in zip([HAN_EX_DJY, HAN_EY_DJX],
                                           ['HAN_EX_DJY', 'HAN_EY_DJX'],
                                           [ax_ex, ax_ey]):
                ax_twin = ax.twinx()
                data_r = data_r.filter('bandpass', freqmin=0.001, freqmax=100, zerophase=True)
                data_r = data_r.slice(starttime, endtime)
                data_r = data_r.remove_response(inventory=response, pre_filt=None, output="DEF", water_level=60,plot=False)
                data_proc = data_r.copy().detrend('simple').taper(0.1).filter('bandpass', freqmin=freqmin,freqmax=freqmax, zerophase=True)
                ax_twin.plot(data_proc[0].times(), data_proc[0].data, color='red', label=label_r)

                # ax.legend(loc='upper right')
                ax.margins(0)
                #ax.xaxis_date()
                ax.text(0.80, 0.1, f'{freqmin}-{freqmax} Hz', transform=ax_twin.transAxes)
            start = data_proc[0].stats.starttime
        formatted_time_new = start.strftime('%Y-%m-%dT%H:%M:%S')
        for ax, title, scale in zip([axs_Ex.flat, axs_Ey.flat],
                                     ['Ex', 'Ey'],
                                     [scale_R,scale_R,scale_R,scale_R]):
            #ax[0].set_title(f'HAN_{folder}_{title}, M{round(float(magnitude[0]),1)}, baz = {round(baz[2],2)}°, acc_Z vs rot_H*{scale}')
            ax[0].set_title(f'HAN @ {title} | {formatted_time_new} | M{round(float(magnitude[0]),1)} | BAz = {round(baz[2],2)}° | acc_Z vs rot_H*{scale}')
            ax[-1].set_xlabel(f'Time [s]')
        for fig, name in zip([fig_Ey, fig_Ex], ['Ey', 'Ex']):
            # fig.savefig(f'{root_savefig}/{folder}/LI_HAN_{folder}_{name}_dispersion_all_stations.png', dpi=150, bbox_inches='tight')
            fig.supylabel(f'Velocity [$\mu$m]')
            supylabel2(fig, f'Angle [nrad]', color='red')
            fig.savefig('%s/LI_HAN_%s_%s_dispersion_all_stations.png' % (root_savefig, folder, name))  # , dpi=10,bbox_inches='tight')
        #plt.show()
        # exit()