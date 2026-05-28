import numpy
from numpy import pi
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import MultipleLocator
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
#event_all = pd.read_csv('../4_Paper_Tables/events_2019_2025_updated.txt')

# import data:
# BRS
folder = '2024-06-28' #'2024-09-15'
for folder in os.listdir(root_data):
    if 'han' in folder:
        continue
    if folder in ['2024-07-31','2024-09-23','2024-05-28','2024-07-02','2025-01-20','2023-08-20',
                  '2022-01-15','.DS_Store','2024-05-13_liv','2024-12-05','2025-06-20_lownoise_L',
                  '2025-02-16_highnoise_L', '2025-02-15_liv','2022-09-19_liv']:
        continue

    print(folder)
    data_event = events[events['location'] == folder.split('_')[0]]


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
    location = 'LIV'
    LIV_EY_DJX, LIV_EX_DJY, LIV_IY_DJX, LIV_IX_DJY, LIV_IX_DH, LIV_IY_DH, LIV_EX_DH, LIV_EY_DH = get_ligo_data(folder,location)
    _, _, LIV_EY_DHZ_orig = LIV_EY_DH.select(channel='*HX'), LIV_EY_DH.select(channel='*HY'), LIV_EY_DH.select(channel='*HZ')

    # rotate to EN
    channel = response.select(network='LI', station='LIV', location='EX', channel='*JY')[0][0][0]
    az = channel.azimuth/180*pi
    LIV_E_DJ = numpy.array([LIV_EY_DJX[0].data, LIV_EX_DJY[0].data, LIV_EX_DJY[0].data * 0])
    LIV_E_DJ_baz = rot_vec(0, 0, az, LIV_E_DJ)
    LIV_E_E_orig = LIV_EY_DJX.copy()
    LIV_E_N_orig = LIV_EX_DJY.copy()
    LIV_E_E_orig[0].data = LIV_E_DJ_baz[0, :]
    LIV_E_N_orig[0].data = LIV_E_DJ_baz[1, :]

    mid_freq = 0.025
    step_freq =0.002
    # remove response, slice, detrend, filter
    for i in range(1,2):

        freqmin, freqmax = 0.01,0.1#round(mid_freq-step_freq*i,3), round(mid_freq+step_freq*i,3)
        print(freqmin, freqmax)
        LIV_E_E = LIV_E_E_orig.copy().slice(starttime, endtime).filter('lowpass',freq=4.5,zerophase=True).detrend('simple').taper(0.1).resample(10).filter('bandpass', freqmin=freqmin, freqmax=freqmax,zerophase=True)
        LIV_E_N = LIV_E_N_orig.copy().slice(starttime, endtime).filter('lowpass',freq=4.5,zerophase=True).detrend('simple').taper(0.1).resample(10).filter('bandpass', freqmin=freqmin, freqmax=freqmax,zerophase=True)
        LIV_EY_DHZ = LIV_EY_DHZ_orig.copy().slice(starttime, endtime).filter('lowpass',freq=4.5,zerophase=True).detrend('simple').taper(0.1).resample(10).filter('bandpass', freqmin=freqmin, freqmax=freqmax,zerophase=True)

        ########################################################################################################################

        sampling_rate = int(LIV_E_N[0].stats.sampling_rate)
        sec = 1/freqmin*2  # window length for correlation
        overlap = 0.5
        win_samples = int(sampling_rate * sec)
        hop_samples = int(win_samples * (1.0 - overlap))
        if hop_samples <= 0:
            raise ValueError("overlap too large -> hop_samples <= 0")

        data_len = len(LIV_E_N[0])  # total samples
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
            rotated = rotate_ne_rt(LIV_E_N[0].data, LIV_E_E[0].data, backas[i_deg])[1] * scale_R * righthandrule

            for i_win in range(0, num_windows):
                start = i_win * hop_samples
                end = start + win_samples
                # grab the windowed slices
                z_slice = LIV_EY_DHZ[0][start:end]
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
        time = numpy.linspace(0, len(LIV_E_N[0].data)/sampling_rate,len(LIV_E_N[0].data)) ## Time-array
        # vertical velocity
        fig = plt.figure(figsize=(9, 5))
        ax1 = plt.subplot2grid((2, 50), (0, 0), colspan=49)
        ax1.plot(time, LIV_EY_DHZ[0].data, 'k', label='vertical velocity')
        ax1.set_xlim(0, time[-1])
        ax1.set_ylabel('[$nm/s$]')
        ax1.legend(loc='lower left')

        '''# add P- and S-wave arrivals
        plt.axvline(arriv_p);plt.annotate('P-arrival', xy=(arriv_p+20,numpy.max(LIV_EY_DHZ[0].data)),xycoords='data');
        plt.axvline(arriv_s);plt.annotate('S-arrival', xy=(arriv_s+20,numpy.max(LIV_EY_DHZ[0].data)),xycoords='data');'''

        # transverse rotation
        r,t = rotate_ne_rt(LIV_E_N[0].data, LIV_E_E[0].data, baz[2])
        #plt.subplot2grid((3, 30), (1, 0), colspan=29)
        ax1_twin = ax1.twinx()
        ax1_twin.plot(time, r, linestyle='--', color='lightpink',label='radial angle')
        ax1_twin.plot(time, t, 'red',label='transversal angle')
        ax1_twin.set_xlim(0, time[-1])
        ax1_twin.set_ylabel('[nrad]', color='red')
        ax1_twin.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
        ax1_twin.legend(loc='lower right')
        '''plt.axvline(arriv_p)
        plt.axvline(arriv_s)'''


        # backazimuth estimation plot
        time_calc = numpy.array([((i * hop_samples + win_samples / 2) / sampling_rate) for i in range(num_windows)])
        # old time_calc = numpy.arange(sec/2., sec * num_windows, sec)
        plt.subplot2grid((2, 50), (1, 0), colspan=49)
        im = plt.pcolor(X, Y, corrbaz, cmap=plt.cm.RdYlGn_r, vmin=-1,vmax=1)
        plt.plot(time_calc, baz_maxcorr, '.k')
        plt.plot(time_calc, baz_midpoints, '*b')
        plt.xlim(0, time[-1])
        plt.ylim(0, 360)
        plt.ylabel(u'estimated \nbackazimuth [°]')
        plt.xlabel('time [s]')
        ax = plt.gca()
        ax.yaxis.set_major_locator(MultipleLocator(60))
        ax.yaxis.set_minor_locator(MultipleLocator(10))

        # plot theoretical Backazimuth for comparison
        xx = [0, time[-1]]
        tba = [baz[2], baz[2]]
        plt.plot(xx, tba, c='.5', lw=1.5)
        plt.text(0.1, 0.1, u'Theor. BAz = '+str(round(baz[2],2))+'°',
                 color='k',fontsize=12,backgroundcolor='.9',transform=plt.gca().transAxes)


        # add colorbar
        cax = plt.subplot2grid((2, 50), (1, 49))
        norm = mpl.colors.Normalize(vmin=-1, vmax=1)
        cb1 = mpl.colorbar.ColorbarBase(cax, cmap=plt.cm.RdYlGn_r, norm=norm, orientation='vertical')

        plt.suptitle(f'{location} external | {folder} | f={freqmin}-{freqmax}')

        fig.savefig(f'{root_savefig}/BAZ_{folder}_{freqmin}_{freqmax}_LIV_external.png', dpi=300, bbox_inches='tight')

        # Save correlation and BAz data

        time = time_calc # numpy.array([((i * hop_samples + win_samples / 2) / sampling_rate) for i in range(num_windows)])
        corr_value_max = max_corr_value
        corr_value_midpoint = corr_midpoints
        baz_value_max = baz_maxcorr
        baz_value_midpoint = baz_midpoints
        true_baz = time*0+round(baz[2],2)
        # max values of rot and trans
        max_rot = max(t)
        max_trans = max(LIV_EY_DHZ[0].data)
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
        df.to_csv(f'{root_savedata}/{folder}_{freqmin}_{freqmax}_LIV_external_corr_baz.csv')
        #plt.show()

#plt.show()