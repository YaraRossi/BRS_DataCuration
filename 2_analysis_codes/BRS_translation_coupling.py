import numpy
from debugpy.launcher import channel
from obspy import read, UTCDateTime, read_inventory

import matplotlib.pyplot as plt



def get_ligo_data(folder, location):
    root_data = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data/miniseed_LIGO'
    root_xml = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data/metadata'

    if location == 'LIV':
        ## 1. READ DATA
        # BRS
        LIV_EY_DJX = read('%s/%s/LI_LIV_EY_DJX.ms' %(root_data,folder))
        LIV_EX_DJY = read('%s/%s/LI_LIV_EX_DJY.ms' %(root_data,folder))
        LIV_IY_DJX = read('%s/%s/LI_LIV_IY_DJX.ms' %(root_data,folder))
        LIV_IX_DJY = read('%s/%s/LI_LIV_IX_DJY.ms' %(root_data,folder))

        # seis
        LIV_EY_DH = read('%s/%s/LI_LIV_EY_DHX.ms' %(root_data,folder))
        LIV_EY_DH += read('%s/%s/LI_LIV_EY_DHY.ms' %(root_data,folder))
        LIV_EY_DH += read('%s/%s/LI_LIV_EY_DHZ.ms' %(root_data,folder))
        LIV_EX_DH = read('%s/%s/LI_LIV_EX_DHX.ms' %(root_data,folder))
        LIV_EX_DH += read('%s/%s/LI_LIV_EX_DHY.ms' %(root_data,folder))
        LIV_EX_DH += read('%s/%s/LI_LIV_EX_DHZ.ms' %(root_data,folder))

        LIV_IY_DH = read('%s/%s/LI_LIV_IY_DHX.ms' %(root_data,folder))
        LIV_IY_DH += read('%s/%s/LI_LIV_IY_DHY.ms' %(root_data,folder))
        LIV_IY_DH += read('%s/%s/LI_LIV_IY_DHZ.ms' %(root_data,folder))
        LIV_IX_DH = read('%s/%s/LI_LIV_IX_DHX.ms' %(root_data,folder))
        LIV_IX_DH += read('%s/%s/LI_LIV_IX_DHY.ms' %(root_data,folder))
        LIV_IX_DH += read('%s/%s/LI_LIV_IX_DHZ.ms' %(root_data,folder))

        response = read_inventory('%s/LIGO_HAN_LIV.xml' % root_xml)

        ## 3. REMOVE RESPONSE
        # BRS
        LIV_EY_DJX.remove_response(inventory=response, pre_filt=None, output="DEF", water_level=60, plot=False)
        LIV_EX_DJY.remove_response(inventory=response, pre_filt=None, output="DEF", water_level=60, plot=False)
        LIV_IY_DJX.remove_response(inventory=response, pre_filt=None, output="DEF", water_level=60, plot=False)
        LIV_IX_DJY.remove_response(inventory=response, pre_filt=None, output="DEF", water_level=60, plot=False)

        # SEIS
        LIV_EY_DH.remove_response(inventory=response, pre_filt=None, output="DEF", water_level=60, plot=False)
        LIV_EX_DH.remove_response(inventory=response, pre_filt=None, output="DEF", water_level=60, plot=False)
        LIV_IY_DH.remove_response(inventory=response, pre_filt=None, output="DEF", water_level=60, plot=False)
        LIV_IX_DH.remove_response(inventory=response, pre_filt=None, output="DEF", water_level=60, plot=False)

        return LIV_EY_DJX, LIV_EX_DJY, LIV_IY_DJX, LIV_IX_DJY, LIV_IX_DH, LIV_IY_DH, LIV_EX_DH, LIV_EY_DH

    if location == 'HAN':
        # BRS
        HAN_EY_DJX = read('%s/%s/LI_HAN_EY_DJX.ms' % (root_data, folder))
        HAN_EX_DJY = read('%s/%s/LI_HAN_EX_DJY.ms' % (root_data, folder))

        # seis
        HAN_EY_DH = read('%s/%s/LI_HAN_EY_DHX.ms' % (root_data, folder))
        HAN_EY_DH += read('%s/%s/LI_HAN_EY_DHY.ms' % (root_data, folder))
        HAN_EY_DH += read('%s/%s/LI_HAN_EY_DHZ.ms' % (root_data, folder))
        HAN_EX_DH = read('%s/%s/LI_HAN_EX_DHX.ms' % (root_data, folder))
        HAN_EX_DH += read('%s/%s/LI_HAN_EX_DHY.ms' % (root_data, folder))
        HAN_EX_DH += read('%s/%s/LI_HAN_EX_DHZ.ms' % (root_data, folder))

        HAN_IY_DH = read('%s/%s/LI_HAN_IY_DHX.ms' % (root_data, folder))
        HAN_IY_DH += read('%s/%s/LI_HAN_IY_DHY.ms' % (root_data, folder))
        HAN_IY_DH += read('%s/%s/LI_HAN_IY_DHZ.ms' % (root_data, folder))
        HAN_VA_DH = read('%s/%s/LI_HAN_VA_DHX.ms' % (root_data, folder))
        HAN_VA_DH += read('%s/%s/LI_HAN_VA_DHY.ms' % (root_data, folder))

        response = read_inventory('%s/LIGO_HAN_LIV.xml' % root_xml)

        ## 3. REMOVE RESPONSE
        # BRS
        HAN_EY_DJX.remove_response(inventory=response, pre_filt=None, output="DEF", water_level=60, plot=False)
        HAN_EX_DJY.remove_response(inventory=response, pre_filt=None, output="DEF", water_level=60, plot=False)

        # SEIS
        HAN_EY_DH.remove_response(inventory=response, pre_filt=None, output="DEF", water_level=60, plot=False)
        HAN_EX_DH.remove_response(inventory=response, pre_filt=None, output="DEF", water_level=60, plot=False)
        HAN_IY_DH.remove_response(inventory=response, pre_filt=None, output="DEF", water_level=60, plot=False)
        HAN_VA_DH.remove_response(inventory=response, pre_filt=None, output="DEF", water_level=60, plot=False)

        return HAN_EY_DJX, HAN_EX_DJY, HAN_EY_DH, HAN_EX_DH, HAN_IY_DH, HAN_VA_DH



folder = 'Jun28.2024'
location = 'LIV'
root_savefig = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/2_Figures/Tcouple'

LIV_EY_DJX, LIV_EX_DJY, LIV_IY_DJX, LIV_IX_DJY, LIV_IX_DH, LIV_IY_DH, LIV_EX_DH, LIV_EY_DH = get_ligo_data(folder, location)

LIV_EY_DJX = LIV_EY_DJX.detrend('simple').filter('bandpass', freqmin=0.008, freqmax=100, zerophase=True).taper(0.05)
LIV_EX_DJY = LIV_EX_DJY.detrend('simple').filter('bandpass', freqmin=0.008, freqmax=100, zerophase=True).taper(0.05)
LIV_IY_DJX = LIV_IY_DJX.detrend('simple').filter('bandpass', freqmin=0.008, freqmax=100, zerophase=True).taper(0.05)
LIV_IX_DJY = LIV_IX_DJY.detrend('simple').filter('bandpass', freqmin=0.008, freqmax=100, zerophase=True).taper(0.05)

LIV_EX_DH_disp = LIV_EX_DH.detrend('simple').filter('bandpass', freqmin=0.008, freqmax=100, zerophase=True).integrate().filter('bandpass', freqmin=0.008, freqmax=100, zerophase=True).taper(0.05)
LIV_EY_DH_disp = LIV_EY_DH.detrend('simple').filter('bandpass', freqmin=0.008, freqmax=100, zerophase=True).integrate().filter('bandpass', freqmin=0.008, freqmax=100, zerophase=True).taper(0.05)
LIV_IX_DH_disp = LIV_IX_DH.detrend('simple').filter('bandpass', freqmin=0.008, freqmax=100, zerophase=True).integrate().filter('bandpass', freqmin=0.008, freqmax=100, zerophase=True).taper(0.05)
LIV_IY_DH_disp = LIV_IY_DH.detrend('simple').filter('bandpass', freqmin=0.008, freqmax=100, zerophase=True).integrate().filter('bandpass', freqmin=0.008, freqmax=100, zerophase=True).taper(0.05)

for disp, angle, chh, cha, station in zip([LIV_EX_DH_disp,LIV_EY_DH_disp,LIV_IX_DH_disp,LIV_IY_DH_disp],
                                           [LIV_EY_DJX, LIV_EX_DJY, LIV_IY_DJX, LIV_IX_DJY],
                                           ['DHY','DHX','DHY','DHX'],
                                           ['DJX','DJY','DJX','DJY'],
                                           ['EY','EX','IY','IX']):

    color = ['black', 'red', 'pink', 'lightpink']
    label = [chh, cha, 'T coupling max', 'T coupling min']
    # Create the figure
    fig = plt.figure(figsize=(10, 8))
    ax_top = fig.add_subplot(2, 1, 1)
    ax_bot_l = fig.add_subplot(2, 2, 3)
    ax_bot_r = fig.add_subplot(2, 2, 4)

    ax_top.plot(disp.select(channel=chh)[0].times(), disp.select(channel=chh)[0].data, label=chh, color='black')
    ax_twin = ax_top.twinx()
    ax_twin.plot(angle[0].times(), angle[0].data, label=cha, color='red')
    ax_twin.plot(disp.select(channel=chh)[0].times(), disp.select(channel=chh)[0].data*1e-5, label=label[2], color='pink')
    ax_twin.plot(disp.select(channel=chh)[0].times(), disp.select(channel=chh)[0].data*1e-6, label=label[3], color='lightpink')

    # fft's
    for data, ax, scale, color, label, channel in zip([disp, angle, disp, disp],
                                [ax_bot_l, ax_bot_r, ax_bot_r, ax_bot_r],
                                [1,1,1e-5,1e-6], color, label, [chh, cha, chh, chh]):
        N = len(data[0].data)
        fft_signal = numpy.fft.fft(data.select(channel=channel)[0].data)
        freqs = numpy.fft.fftfreq(N, d=data.select(channel=channel)[0].stats.delta)  # Frequency vector
        # Take the positive frequencies only
        positive_freqs = freqs[:N // 2]
        magnitude = numpy.abs(fft_signal[:N // 2]) / N * scale  # Normalize
        ax.loglog(positive_freqs, magnitude, color=color, label=label)

    ax_top.set_ylabel('disp [m]')
    ax_twin.set_ylabel('angle [nrad]')
    ax_bot_l.set_ylabel('disp [m]')
    ax_bot_r.set_ylabel('angle [nrad]')

    ax_top.set_xlabel('time [s]')
    ax_bot_l.set_xlabel('freq [Hz]')
    ax_bot_r.set_xlabel('freq [Hz]')

    ax_top.legend(loc='upper right')
    ax_twin.legend(loc='lower right')
    ax_bot_r.legend(loc='lower left')
    ax_bot_l.legend(loc='lower left')

    plt.suptitle(f'{folder} | {location} | {station}')
    # Adjust spacing to prevent overlap
    fig.tight_layout()
    fig.savefig(f'{root_savefig}/{location}_Tcouple_{folder}_{station}_{cha}.png', bbox_inches='tight', dpi=300)

plt.show()