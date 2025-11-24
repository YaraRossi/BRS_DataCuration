import numpy
from numpy import pi
import matplotlib.pyplot as plt
import matplotlib as mpl
from obspy import read, UTCDateTime, read_inventory
from obspy.signal.cross_correlation import correlate
from obspy.signal.rotate import rotate_ne_rt
from obspy.geodetics.base import gps2dist_azimuth
from functions import rot_vec, circular_positive_midpoint, get_ligo_data
import pandas as pd
import os


# roots
root_savefig = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/2_Figures/Backazimuth_estimations/fixed_freq'
root_savedata = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data/correlation_baz_csv'
root_data = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data/miniseed_LIGO'
root_xml = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data/metadata'
response = read_inventory('%s/LIGO_HAN_LIV.xml' %root_xml)
path = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data'
events = pd.read_csv(f'{path}/event_times_LIV.txt')
event_all = pd.read_csv('../4_Paper_Tables/events_2019_2025_updated.txt')

# import data:
# BRS
folder = 'Jun28.2024' #'Sept15.2024'
for folder in os.listdir(root_data):
    if folder in ['July31.2024','Sept23.2024','May28.2024','July02.2024','Jan20.2025','Aug20.2023','Jan15.2022','.DS_Store','May13.2024']:
        continue
    data_event = events[events['location'] == folder]
    data_event_all = event_all[event_all['location'] == folder]
    if data_event.empty:
        data_event = data_event_all
        if data_event.empty:
            print(folder, 'event still missing from eventfile)')
            continue

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
    _, _, LIV_IY_DHZ_orig = LIV_IY_DH.select(channel='DHX'), LIV_IY_DH.select(channel='DHY'), LIV_IY_DH.select(channel='DHZ')
    '''LIV_IY_DJX = read('%s/%s/LI_LIV_IY_DJX.ms' % (root_data, folder))
    LIV_IX_DJY = read('%s/%s/LI_LIV_IX_DJY.ms' % (root_data, folder))
    LIV_IY_DHZ_orig = read('%s/%s/LI_LIV_IY_DHZ.ms' % (root_data, folder))'''

    # rotate to EN
    channel = response.select(network='LI', station='LIV', location='IX', channel='DJY')[0][0][0]
    az = channel.azimuth/180*pi
    LIV_I_DJ = numpy.array([LIV_IY_DJX[0].data, LIV_IX_DJY[0].data, LIV_IX_DJY[0].data * 0])
    LIV_I_DJ_baz = rot_vec(0, 0, az, LIV_I_DJ)
    LIV_I_E_orig = LIV_IY_DJX.copy()
    LIV_I_N_orig = LIV_IX_DJY.copy()
    LIV_I_E_orig[0].data = LIV_I_DJ_baz[0, :]
    LIV_I_N_orig[0].data = LIV_I_DJ_baz[1, :]

    mid_freq = 0.025
    step_freq =0.002
    # remove response, slice, detrend, filter
    for i in range(1,2):

        freqmin, freqmax = 0.008,0.2#round(mid_freq-step_freq*i,3), round(mid_freq+step_freq*i,3)
        print(freqmin, freqmax)
        LIV_I_E = LIV_I_E_orig.copy().slice(starttime, endtime).filter('lowpass',freq=10,zerophase=True).detrend('simple').taper(0.02).resample(20).filter('bandpass', freqmin=freqmin, freqmax=freqmax,zerophase=True)
        LIV_I_N = LIV_I_N_orig.copy().slice(starttime, endtime).filter('lowpass',freq=10,zerophase=True).detrend('simple').taper(0.02).resample(20).filter('bandpass', freqmin=freqmin, freqmax=freqmax,zerophase=True)
        LIV_IY_DHZ = LIV_IY_DHZ_orig.copy().slice(starttime, endtime).filter('lowpass',freq=10,zerophase=True).detrend('simple').taper(0.02).resample(20).filter('bandpass', freqmin=freqmin, freqmax=freqmax,zerophase=True)



        ########################################################################################################################

        sampling_rate = int(LIV_I_N[0].stats.sampling_rate)
        sec = 90  # window length for correlation
        overlap = 0.5
        win_samples = int(sampling_rate * sec)
        hop_samples = int(win_samples * (1.0 - overlap))
        if hop_samples <= 0:
            raise ValueError("overlap too large -> hop_samples <= 0")

        data_len = len(LIV_I_N[0])  # total samples
        # number of windows that fully fit (change to ceil if you want last partial window)
        num_windows = 0 if data_len < win_samples else ((data_len - win_samples) // hop_samples + 1)

        # estimate the Backazimuth for each time window
        step = 5
        righthandrule = -1
        backas = numpy.linspace(0, 360 - step, int(360 / step))
        corrbaz = []
        #ind=None
        for i_deg in range(len(backas)):
            # compute rotated trace once per azimuth
            rotated = rotate_ne_rt(LIV_I_N[0].data, LIV_I_E[0].data, backas[i_deg])[1] * scale_R * righthandrule

            for i_win in range(0, num_windows):
                start = i_win * hop_samples
                end = start + win_samples
                # grab the windowed slices
                z_slice = LIV_IY_DHZ[0][start:end]
                rotated_slice = rotated[start:end]

                corrbazz = correlate(z_slice, rotated_slice, 0)
                corrbaz.append(corrbazz[0])

        corrbaz = numpy.asarray(corrbaz)
        corrbaz = corrbaz.reshape(len(backas), num_windows)

        hop_seconds = hop_samples / sampling_rate
        time_calc = numpy.array([((i * hop_samples + win_samples / 2) / sampling_rate) for i in range(num_windows)])
        # calculate the meshgrid ov crosscorrelation for plotting
        X, Y = numpy.meshgrid(time_calc, backas)


        baz_maxcorr = [] # backazimuth of max corr value
        max_corr_value = []
        for l1 in range(0, num_windows):
            maxcor_r = backas[corrbaz[:, l1].argmax()]
            max_corr_value_r = corrbaz[:, l1].max()
            baz_maxcorr.append(maxcor_r)
            max_corr_value.append(max_corr_value_r)
        baz_maxcorr = numpy.asarray(baz_maxcorr)
        max_corr_value = numpy.asarray(max_corr_value)

        # getting middle of positive values
        baz_midpoints = []
        corr_midpoints = []
        corrbaz_numpy = numpy.array(corrbaz)
        for v in corrbaz_numpy.T:
            mid_azimuth = circular_positive_midpoint(v)
            azimuth_midpoint = backas[mid_azimuth]
            corr_midpoint = v[mid_azimuth]
            baz_midpoints.append(azimuth_midpoint)
            corr_midpoints.append(corr_midpoint)



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
        time_calc = numpy.array([((i * hop_samples + win_samples / 2) / sampling_rate) for i in range(num_windows)])
        # old time_calc = numpy.arange(sec/2., sec * num_windows, sec)
        plt.subplot2grid((3, 30), (2, 0), colspan=29)
        im = plt.pcolor(X, Y, corrbaz, cmap=plt.cm.RdYlGn_r, vmin=-1,vmax=1)
        plt.plot(time_calc, baz_maxcorr, '.k')
        plt.plot(time_calc, baz_midpoints, '*b')
        plt.xlim(0, time[-1])
        plt.ylim(0, 360)
        plt.ylabel(u'estimated \nbackazimuth [°]')
        plt.xlabel('time [s]')

        # plot theoretical Backazimuth for comparison
        xx = [0, time[-1]]
        tba = [baz[2], baz[2]]
        plt.plot(xx, tba, c='.5', lw=1.5)
        plt.text(2600, 106, u'Theor. BAz = '+str(round(baz[2],2))+'°', color='k',fontsize=12,backgroundcolor='.9')


        # add colorbar
        cax = plt.subplot2grid((3, 30), (2, 29))
        norm = mpl.colors.Normalize(vmin=-1, vmax=1)
        cb1 = mpl.colorbar.ColorbarBase(cax, cmap=plt.cm.RdYlGn_r, norm=norm, orientation='vertical')

        plt.suptitle(f'f={freqmin}-{freqmax}')

        fig.savefig(f'{root_savefig}/BAZ_{folder}_{freqmin}_{freqmax}.png', dpi=300, bbox_inches='tight')

        # Save correlation and BAz data

        time = time_calc # numpy.array([((i * hop_samples + win_samples / 2) / sampling_rate) for i in range(num_windows)])
        corr_value_max = max_corr_value
        corr_value_midpoint = corr_midpoints
        baz_value_max = baz_maxcorr
        baz_value_midpoint = baz_midpoints
        true_baz = time*0+round(baz[2],2)
        # max values of rot and trans
        max_rot = max(t)
        max_trans = max(LIV_IY_DHZ[0].data)
        data = {'time': time,
                'corr_value_max': corr_value_max,
                'baz_value_max': baz_value_max,
                'corr_value_midpoint': corr_value_midpoint,
                'baz_value_midpoint': baz_value_midpoint,
                'true_baz': true_baz,
                'max_rot': max_rot,
                'max_trans':max_trans
                }
        df = pd.DataFrame(data)
        df.to_csv(f'{root_savedata}/{folder}_{freqmin}_{freqmax}_corr_baz.csv')

plt.show()