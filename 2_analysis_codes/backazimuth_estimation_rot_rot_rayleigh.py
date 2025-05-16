import numpy
import os
from numpy import pi, arctan, mean
import matplotlib.pyplot as plt
import matplotlib as mpl
from obspy import read, UTCDateTime, read_inventory
from obspy.signal.cross_correlation import correlate
from obspy.signal.rotate import rotate_ne_rt
from obspy.geodetics.base import gps2dist_azimuth
from functions import rot_vec
import pandas as pd


# roots
root_savefig = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/2_Figures/Backazimuth_estimations'
root_data = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data/miniseed_LIGO'
root_xml = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data/metadata'
response = read_inventory('%s/LIGO_HAN_LIV.xml' %root_xml)
path = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data'
events = pd.read_csv(f'{path}/event_times_LIV.txt')

# import data:
# BRS
for folder in os.listdir(root_data):
    if folder in ['.DS_Store', 'Aug20.2023', 'Feb15.2025']:
        continue
    '''if 'Sept15.2024' not in folder:
       continue'''

    if folder in ['July31.2024','Sept23.2024','May28.2024','July02.2024','Jan20.2025','Aug20.2023','Jan15.2022']:
        continue

    print(folder)
    data_event = events[events['location'] == folder]
    starttime = UTCDateTime(str(data_event['eventstart'].iloc[0]))
    endtime = UTCDateTime(str(data_event['eventend'].iloc[0]))

    # theoretical baz
    source_latitude = data_event[' latitude'].values.astype(float)
    source_longitude = data_event[' longitude'].values.astype(float)
    station_latitude = 30.562775410180993
    station_longitude = -90.77469251757329
    baz = gps2dist_azimuth(source_latitude, source_longitude, station_latitude, station_longitude)
    scale_R = -1 # because it's rayleigh waves

    # download data
    LIV_IY_DJX = read('%s/%s/LI_LIV_IY_DJX.ms' % (root_data, folder))
    LIV_IX_DJY = read('%s/%s/LI_LIV_IX_DJY.ms' % (root_data, folder))
    LIV_IY_DHZ_orig = read('%s/%s/LI_LIV_IY_DHZ.ms' % (root_data, folder))

    # rotate to EN
    channel = response.select(network='LI', station='LIV', location='IX', channel='DJY')[0][0][0]
    az = channel.azimuth/180*pi
    LIV_I_DJ = numpy.array([LIV_IY_DJX[0].data, LIV_IX_DJY[0].data, LIV_IX_DJY[0].data * 0])
    LIV_I_DJ_baz = rot_vec(0, 0, az, LIV_I_DJ)
    LIV_I_E_orig = LIV_IY_DJX.copy()
    LIV_I_N_orig = LIV_IX_DJY.copy()
    LIV_I_E_orig[0].data = LIV_I_DJ_baz[0, :]
    LIV_I_N_orig[0].data = LIV_I_DJ_baz[1, :]

    # remove response, slice, detrend, filter
    for freqmin, freqmax in zip([0.008, 0.009, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08,0.09,0.008],
                                [0.009,  0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.1]):
        #freqmin, freqmax = 0.04, 0.06
        LIV_I_E = LIV_I_E_orig.copy().remove_response(inventory=response, pre_filt=None, output="DEF", water_level=60, plot=False).slice(starttime, endtime).resample(20).detrend('simple').taper(0.02).filter('bandpass', freqmin=freqmin, freqmax=freqmax,zerophase=True)
        LIV_I_N = LIV_I_N_orig.copy().remove_response(inventory=response, pre_filt=None, output="DEF", water_level=60, plot=False).slice(starttime, endtime).resample(20).detrend('simple').taper(0.02).filter('bandpass', freqmin=freqmin, freqmax=freqmax,zerophase=True)
        LIV_IY_DHZ = LIV_IY_DHZ_orig.copy().remove_response(inventory=response, pre_filt=None, output="DEF", water_level=60, plot=False).slice(starttime, endtime).resample(20).detrend('simple').taper(0.02).filter('bandpass', freqmin=freqmin, freqmax=freqmax,zerophase=True)



        ########################################################################################################################

        sampling_rate = int(LIV_I_N[0].stats.sampling_rate)
        sec = int(1/freqmin*5)  # window length for correlation
        num_windows = len(LIV_I_N[0]) // (int(sampling_rate * sec))

        # estimate the Backazimuth for each time window
        step = 5
        backas = numpy.linspace(0, 360 - step, int(360 / step))
        baz_rr = []
        ind=None
        righthandrule = -1
        for i_win in range(0, num_windows):
            '''corrbazz = correlate(LIV_IY_DHZ[0][sampling_rate * sec * i_win : sampling_rate * sec * (i_win + 1)],
                                 rotate_ne_rt(LIV_I_N[0].data, LIV_I_E[0].data, backas[i_deg])
                                 [1][sampling_rate * sec * i_win : sampling_rate * sec * (i_win + 1)]*scale_R*righthandrule,0)'''
            tanbaz = -arctan(LIV_I_N[0].data[sampling_rate * sec * i_win : sampling_rate * sec * (i_win + 1)]/LIV_I_E[0].data[sampling_rate * sec * i_win : sampling_rate * sec * (i_win + 1)])
            baz_rr.append(mean(tanbaz)/pi*180+180)

        # plotting data:
        time = numpy.linspace(0, len(LIV_I_N[0].data)/sampling_rate,len(LIV_I_N[0].data)) ## Time-array
        # vertical velocity
        fig = plt.figure(figsize=(12, 8))
        plt.subplot2grid((3, 30), (0, 0), colspan=29)
        plt.plot(time, LIV_IY_DHZ[0].data, 'k', label='vertical velocity')
        plt.xlim(0, time[-1])
        plt.ylabel('transv. acc. \n[$nm/s$]')
        plt.legend()

        '''# add P- and S-wave arrivals
        plt.axvline(arriv_p);plt.annotate('P-arrival', xy=(arriv_p+20,numpy.max(LIV_IY_DHZ[0].data)),xycoords='data');
        plt.axvline(arriv_s);plt.annotate('S-arrival', xy=(arriv_s+20,numpy.max(LIV_IY_DHZ[0].data)),xycoords='data');'''

        # transverse rotation
        r,t = rotate_ne_rt(LIV_I_N[0].data, LIV_I_E[0].data, baz[2])
        plt.subplot2grid((3, 30), (1, 0), colspan=29)
        plt.plot(time, r, linestyle='--', color='lightpink',label='radial angle')
        plt.plot(time, t, 'red',label='transversal angle')
        plt.xlim(0, time[-1])
        plt.ylabel('vert. rot. rate \n[nrad]')
        plt.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
        plt.legend(loc='lower left')
        '''plt.axvline(arriv_p)
        plt.axvline(arriv_s)'''


        # backazimuth estimation plot
        plt.subplot2grid((3, 30), (2, 0), colspan=29)
        plt.plot(numpy.arange(sec/2., sec * num_windows, sec), baz_rr, '.k')
        plt.xlim(0, time[-1])
        plt.ylim(0, 360)
        plt.ylabel(u'estimated \nbackazimuth [°]')
        plt.xlabel('time [s]')

        # plot theoretical Backazimuth for comparison
        xx = numpy.arange(0, sec * num_windows + 1, sec)
        tba = numpy.ones(len(xx)) * baz[2]
        plt.plot(xx, tba, c='.5', lw=1.5)
        plt.text(2600, 106, u'Theor. BAz = '+str(round(baz[2],2))+'°', color='k',fontsize=12,backgroundcolor='.9')


        # add colorbar
        cax = plt.subplot2grid((3, 30), (2, 29))
        norm = mpl.colors.Normalize(vmin=-1, vmax=1)
        cb1 = mpl.colorbar.ColorbarBase(cax, cmap=plt.cm.RdYlGn_r, norm=norm, orientation='vertical')

        plt.suptitle(f'f={freqmin}-{freqmax}')

        fig.savefig(f'{root_savefig}/BAZ_rot_rot_{folder}_{freqmin}_{freqmax}.png', dpi=300, bbox_inches='tight')

    plt.show()