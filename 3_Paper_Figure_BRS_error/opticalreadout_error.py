import numpy
from debugpy.launcher import channel
from obspy import read, UTCDateTime, read_inventory
import sys
sys.path.append("../2_analysis_codes")

from functions import get_ligo_data
import matplotlib.pyplot as plt

folder = '2024-06-28'
location = 'LIV'
root_savefig = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/2_Figures/optical_readout_error'

LIV_EY_DJX, LIV_EX_DJY, LIV_IY_DJX, LIV_IX_DJY, LIV_IX_DH, LIV_IY_DH, LIV_EX_DH, LIV_EY_DH = get_ligo_data(folder, location='LIV')
HAN_EY_DJX, HAN_EX_DJY, HAN_EY_DH, HAN_EX_DH, HAN_VA_DH = get_ligo_data(folder, location='HAN')

LIV_EY_DJX = LIV_EY_DJX.detrend('simple').taper(0.001)
LIV_EX_DJY = LIV_EX_DJY.detrend('simple').taper(0.001)
LIV_IY_DJX = LIV_IY_DJX.detrend('simple').taper(0.001)
LIV_IX_DJY = LIV_IX_DJY.detrend('simple').taper(0.001)

HAN_EY_DJX = HAN_EY_DJX.detrend('simple').taper(0.001)
HAN_EX_DJY = HAN_EX_DJY.detrend('simple').taper(0.001)

starttime =LIV_IY_DJX[0].stats.starttime+100
endtime = starttime+10

fig, axs = plt.subplots(6, 1, figsize=(4, 4))

for ax, angle, ch, loc, station in zip( axs.flat,
                                         [LIV_EY_DJX, LIV_EX_DJY, LIV_IY_DJX, LIV_IX_DJY, HAN_EY_DJX, HAN_EX_DJY],
                                         ['DJX','DJY','DJX','DJY','DJX','DJY'],
                                         ['EY','EX','IY','IX','IY','IX'],
                                         ['LIV','LIV','LIV','LIV','HAN','HAN']):
    angle = angle.slice(starttime,endtime).detrend('simple')
    color = 'red'

    ax.plot(angle[0].times(), angle[0].data, label=f'{station}.{loc}.{ch}', color=color)
    ax.legend(loc='lower right')

axs[3].set_ylabel('Amplitude [nrad]')

ax.set_xlabel('Time [s]')


axs[0].set_title('(b)',loc='left')
fig.tight_layout()
plt.subplots_adjust(hspace=0)
fig.savefig(f'{root_savefig}/optical_readout_error_4paper.png', bbox_inches='tight', dpi=300)

plt.show()