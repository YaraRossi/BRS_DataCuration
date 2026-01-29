import numpy
from IPython.core.pylabtools import figsize
from debugpy.launcher import channel
import matplotlib.dates as mdates
from obspy import read, UTCDateTime, read_inventory
from functions import get_ligo_data
import os
import matplotlib.pyplot as plt
import scipy.signal as signal


root_savefig = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/2_Figures/2_Raw_PSD'
root_data = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data/miniseed_LIGO'

for folder in os.listdir(root_data):
    if folder in ['.DS_Store']:
        continue
    if folder not in ['2026-01-16','2023-09-08_L','2022-09-19_liv','2022-09-19_han_L']:
        continue

    LIV, HAN = True, True

    if folder in ['2024-07-31','2024-09-23','2024-05-28','2024-07-02','2025-01-20','2023-08-20','2022-01-15',
                  '2021-07-29_L']:
        LIV = False
    if folder in ['2024-05-13', '2025-02-15']:
        HAN = False
    if 'han' in folder:
        LIV = False
    if 'liv' in folder:
        HAN = False

    for location, runornot in zip(['LIV','HAN'],[LIV,HAN]):
        if not runornot:
            continue

        if location == 'LIV':
            LIV_EY_DJX, LIV_EX_DJY, LIV_IY_DJX, LIV_IX_DJY, LIV_IX_DH, LIV_IY_DH, LIV_EX_DH, LIV_EY_DH = get_ligo_data(folder, location)
            streams = LIV_EY_DJX, LIV_EX_DJY, LIV_IY_DJX, LIV_IX_DJY, LIV_IX_DH, LIV_IY_DH, LIV_EX_DH, LIV_EY_DH
            name = ['LIV_EY_DJX', 'LIV_EX_DJY', 'LIV_IY_DJX', 'LIV_IX_DJY', 'LIV_IX_DH', 'LIV_IY_DH', 'LIV_EX_DH', 'LIV_EY_DH']
            nrows=8
        elif location == 'HAN':
            HAN_EY_DJX, HAN_EX_DJY, HAN_EY_DH, HAN_EX_DH, HAN_VA_DH = get_ligo_data(folder, location)
            streams = HAN_EY_DJX, HAN_EX_DJY, HAN_EY_DH, HAN_EX_DH, HAN_VA_DH
            name = ['HAN_EY_DJX', 'HAN_EX_DJY', 'HAN_EY_DH', 'HAN_EX_DH', 'HAN_VA_DH']
            nrows=5


        starttime = streams[0][0].stats.starttime
        endtime = streams[0][0].stats.endtime
        zoom_start = starttime+(endtime-starttime)/4
        zoom_end = zoom_start+(endtime-starttime)/4/2
        d = starttime.date

        freqmin, freqmax = 0.01, 20


        fig, axs = plt.subplots(nrows,2, figsize=(8,nrows), sharex='col')
        plt.subplots_adjust(hspace=0)

        for i, seed in enumerate(streams):
            if 'VA' in name[i]:
                if seed == None:
                    continue

            # raw data
            ax = axs[i,0]
            if 'J' in name[i]:
                color = 'red'
                f_t, P_t = signal.welch(seed[0].data, 1 / seed[0].stats.delta, nperseg=int(len(seed[0].data)/2))
                ax.loglog(f_t, P_t, color=color, label=name[i][4:])
            else:
                color = 'k'
                if 'VA' in name[i]:
                    f_t, P_t = signal.welch(seed.select(channel=f'DHX')[0].data, 1 / seed.select(channel=f'DHX')[0].stats.delta, nperseg=int(len(seed[0].data)/2))
                    ax.loglog(f_t, P_t, color=color, label=f'{name[i][4:]}X')
                else:
                    f_t, P_t = signal.welch(seed.select(channel=f'DH{name[i][5]}')[0].data, 1 / seed.select(channel=f'DH{name[i][5]}')[0].stats.delta, nperseg=int(len(seed[0].data)/2))
                    ax.loglog(f_t, P_t, color=color, label=f'{name[i][4:]}{name[i][5]}')
            ax.legend(loc='lower right')
            ax.margins(0)


            # zoomed and bandpass filtered
            ax = axs[i, 1]
            seed_f = seed.detrend('demean').taper(0.01).copy().filter('bandpass', freqmin=freqmin, freqmax=freqmax)
            seed_f_s = seed_f.slice(zoom_start, zoom_end)
            if 'J' in name[i]:
                color = 'red'
                f_t, P_t = signal.welch(seed_f_s[0].data, 1 / seed_f_s[0].stats.delta, nperseg=int(len(seed_f_s[0].data)/2))
                ax.loglog(f_t, P_t, color=color, label=name[i][4:])
            else:
                color = 'k'
                if 'VA' in name[i]:
                    f_t, P_t = signal.welch(seed_f_s.select(channel=f'DHX')[0].data, 1 / seed_f_s.select(channel=f'DHX')[0].stats.delta, nperseg=int(len(seed_f_s[0].data)/2))
                    ax.loglog(f_t, P_t, color=color, label=f'{name[i][4:]}X')
                else:
                    f_t, P_t = signal.welch(seed_f_s.select(channel=f'DH{name[i][5]}')[0].data, 1 / seed_f_s.select(channel=f'DH{name[i][5]}')[0].stats.delta, nperseg=int(len(seed_f_s[0].data)/2))
                    ax.loglog(f_t,P_t, color=color, label=f'{name[i][4:]}{name[i][5]}')
            #ax.legend(loc='lower right')
            ax.margins(0)


        axs[0,0].set_title('Raw data')
        axs[0,1].set_title(f'Zoom @ {freqmin}-{freqmax}Hz')

        # setting x label format
        axs[-1,0].tick_params(axis='x', rotation=30)
        axs[-1,1].tick_params(axis='x', rotation=30)
        axs[-1,0].set_xlabel(f'Frequency Hz')
        axs[-1,1].set_xlabel(f'Frequency Hz')

        fig.savefig(f'{root_savefig}/PSD_{folder}_{location}.png', dpi=300, bbox_inches='tight')

        #plt.show()
