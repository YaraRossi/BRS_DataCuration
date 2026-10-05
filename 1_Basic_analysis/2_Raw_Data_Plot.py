import numpy
from IPython.core.pylabtools import figsize
from debugpy.launcher import channel
import matplotlib.dates as mdates
from obspy import read, UTCDateTime, read_inventory
from functions import get_ligo_data, get_UW_data
import os
import matplotlib.pyplot as plt

root_savefig = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/2_Figures/1_Raw_Filter'
root_data = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data/zenodo_upload/miniseed_LIGO'

def plot_all_traces(ax, stream, color, label):
    """Plot every trace in a stream individually (no merging), so gaps stay visible."""
    for j, tr in enumerate(stream):
        ax.plot(tr.times('matplotlib'), tr.data, color=color, label=label if j == 0 else None)

for folder in os.listdir(root_data):
    if folder in ['.DS_Store']:
        continue
    if folder not in ['2022-12-27_han_highnoise_L']:
        continue

    LIV, HAN = True, True

    if folder in ['2024-07-31','2024-09-23','2024-05-28','2024-07-02','2025-01-20','2023-08-20','2022-01-15','2021-07-29_L']:
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
            UW_NIKE_HH, UW_NIKE_EN, UW_LIGO_HH = get_UW_data(folder, location)
            streams = HAN_EY_DJX, HAN_EX_DJY, HAN_EY_DH, HAN_EX_DH, HAN_VA_DH, UW_NIKE_HH, UW_NIKE_EN, UW_LIGO_HH
            name = ['HAN_EY_DJX', 'HAN_EX_DJY', 'HAN_EY_DH', 'HAN_EX_DH', 'HAN_VA_DH', 'UW_NIKE_HH', 'UW_NIKE_EN', 'UW_LIGO_HH']
            nrows=8


        starttime = streams[0][0].stats.starttime
        endtime = streams[0][0].stats.endtime
        fmt = "%Y-%m-%dT%H:%M:%S"
        
        # Print start and end times for each station/variable
        print(f"\n=== Event: {folder}, Location: {location} ===")
        for seed, stream_name in zip(streams, name):
            if seed is not None:
                print(f'{stream_name}: {seed[0].stats.starttime.strftime(fmt)} - {seed[0].stats.endtime.strftime(fmt)}')
            else:
                print(f'{stream_name}: No data available')
        print()
        
        zoom_start = starttime+(endtime-starttime)/4
        zoom_end = zoom_start+(endtime-starttime)/4/2
        d = starttime.date

        freqmin, freqmax = 0.01, 20


        fig, axs = plt.subplots(nrows,2, figsize=(8,nrows), sharex='col')
        plt.subplots_adjust(hspace=0)

        for i, seed in enumerate(streams):
            # Skip if data is None (for missing HAN_VA_DH or UW data)
            if seed == None:
                continue

            # raw data
            ax = axs[i,0]
            if 'J' in name[i]:
                color = 'red'
                plot_all_traces(ax, seed, color, name[i][4:])
            elif 'VA' in name[i]:
                color = 'k'
                plot_all_traces(ax, seed.select(channel='*HX'), color, f'{name[i][4:]}X')
            elif 'UW' in name[i]:
                color = 'blue'
                plot_all_traces(ax, seed.select(channel='*E'), color, f'{name[i]}E')
            else:
                color = 'k'
                plot_all_traces(ax, seed.select(channel=f'*H{name[i][5]}'), color, f'{name[i][4:]}{name[i][5]}')
            ax.legend(loc='lower right')
            ax.margins(0)


            # zoomed and bandpass filtered
            ax = axs[i, 1]
            seed_f = seed.detrend('demean').taper(0.01).copy().filter('bandpass', freqmin=freqmin, freqmax=freqmax)
            seed_f_s = seed_f.slice(zoom_start, zoom_end)
            if 'J' in name[i]:
                color = 'red'
                plot_all_traces(ax, seed_f_s, color, name[i][4:])
            elif 'VA' in name[i]:
                color = 'k'
                plot_all_traces(ax, seed_f_s.select(channel='*HX'), color, f'{name[i][4:]}X')
            elif 'UW' in name[i]:
                color = 'blue'
                plot_all_traces(ax, seed_f_s.select(channel='*E'), color, f'{name[i]}E')
            else:
                color = 'k'
                plot_all_traces(ax, seed_f_s.select(channel=f'*H{name[i][5]}'), color, f'{name[i][4:]}{name[i][5]}')
            #ax.legend(loc='lower right')
            ax.margins(0)


        axs[0,0].set_title('Raw data')
        axs[0,1].set_title(f'Zoom @ {freqmin}-{freqmax}Hz')

        # setting x label format
        time_fmt = mdates.DateFormatter('%H:%M')
        axs[-1,0].xaxis.set_major_formatter(time_fmt)
        axs[-1,1].xaxis.set_major_formatter(time_fmt)
        axs[-1,0].tick_params(axis='x', rotation=30)
        axs[-1,1].tick_params(axis='x', rotation=30)
        axs[-1,0].set_xlabel(f'Time of Day {d}')
        axs[-1,1].set_xlabel(f'time of Day {d}')

        fig.savefig(f'{root_savefig}/{folder}_{location}.png', dpi=300, bbox_inches='tight')

        plt.show()
