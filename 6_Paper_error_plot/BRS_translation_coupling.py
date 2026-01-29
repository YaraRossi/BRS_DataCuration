import numpy
from debugpy.launcher import channel
from obspy import read, UTCDateTime, read_inventory
import sys
sys.path.append("../2_analysis_codes")

from functions import get_ligo_data
import matplotlib.pyplot as plt



folder = '2024-06-28'
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

    color = ['black', 'grey', 'deeppink', 'lightpink']
    label = [f'translation {chh}', f'rotation {cha}', 'rotation from\ntranslation coupling\n1e-5 rad/m', 'rotation from\ntranslation coupling\n1e-6 rad/m']
    # Create the figure
    fig = plt.figure(figsize=(10, 8))
    ax_top = fig.add_subplot(2, 1, 1)
    ax_bot_l = fig.add_subplot(2, 2, 3)
    ax_bot_r = fig.add_subplot(2, 2, 4)

    fig_pap, ax_pap = plt.subplots(1,1,figsize=(4, 4))

    ax_top.plot(disp.select(channel=chh)[0].times(), disp.select(channel=chh)[0].data, label=chh, color='black')
    ax_twin = ax_top.twinx()
    ax_twin.plot(angle[0].times(), angle[0].data, label=cha, color='red')
    ax_twin.plot(disp.select(channel=chh)[0].times(), disp.select(channel=chh)[0].data*1e-5, label=label[2], color='deeppink')
    ax_twin.plot(disp.select(channel=chh)[0].times(), disp.select(channel=chh)[0].data*1e-6, label=label[3], color='lightpink')

    # fft's
    for data, ax, scale, color, label, channel in zip([disp, angle, disp, disp],
                                [ax_bot_l, ax_bot_r, ax_bot_r, ax_bot_r],
                                [1,1*1e-9,1e-5,1e-6], color, label, [chh, cha, chh, chh]):
        N = len(data[0].data)
        fft_signal = numpy.fft.fft(data.select(channel=channel)[0].data)
        freqs = numpy.fft.fftfreq(N, d=data.select(channel=channel)[0].stats.delta)  # Frequency vector
        # Take the positive frequencies only
        positive_freqs = freqs[:N // 2]
        magnitude = numpy.abs(fft_signal[:N // 2]) / N * scale  # Normalize
        ax.loglog(positive_freqs, magnitude, color=color, label=label)

        ax_pap.loglog(positive_freqs, magnitude, color=color, label=label)


    ax_top.set_ylabel('disp [m]')

    ax_twin.set_ylabel('angle [nrad]')
    ax_bot_l.set_ylabel('disp [m]')
    ax_bot_r.set_ylabel('angle [rad]') # check scale if in nrad or rad
    ax_pap.set_ylabel('angle [rad] or disp [m]')

    ax_top.set_xlabel('time [s]')
    ax_bot_l.set_xlabel('freq [Hz]')
    ax_bot_r.set_xlabel('freq [Hz]')
    ax_pap.set_xlabel('freq [Hz]')

    ax_top.legend(loc='upper right')
    ax_twin.legend(loc='lower right')
    ax_bot_r.legend(loc='lower left')
    ax_bot_l.legend(loc='lower left')
    ax_pap.legend(loc='upper right')

    fig.suptitle(f'{folder} | {location} | {station}')
    # Adjust spacing to prevent overlap
    fig.tight_layout()
    fig.savefig(f'{root_savefig}/{location}_Tcouple_{folder}_{station}_{cha}.png', bbox_inches='tight', dpi=300)

    ax_pap.set_title('(a)',loc='left')
    fig_pap.tight_layout()
    ax_pap.set_ylim(top = 1e-5, bottom=1e-14)
    ax_pap.set_xlim(left=0.003, right=100)
    fig_pap.savefig(f'{root_savefig}/{location}_Tcouple_{folder}_{station}_{cha}_4paper.png', bbox_inches='tight', dpi=300)

plt.show()