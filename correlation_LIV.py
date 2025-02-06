import matplotlib.pyplot as plt
import obspy
from obspy import read, UTCDateTime
from correlation_function import yaracorrelate2, yaracorrelate

root_data_HAN = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data/HAN_APR_02'
root_data_LIV = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data/LIV_APR_02'
root_savefig = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/2_Figures/2024Apr02'


start = UTCDateTime('2024-0403T00:15:00')
end = UTCDateTime('2024-0403T00:45:00')

# BRS
LIV_ES_DJ1 = read('%s/LIV_ES_DJ1_Apr_02_2024_23-58-12-000_256.csv.ms' %root_data_LIV)
LIV_EW_DJ2 = read('%s/LIV_EW_DJ2_Apr_02_2024_23-58-12-000_256.csv.ms' %root_data_LIV)
LIV_IS_DJ1 = read('%s/LIV_IS_DJ1_Apr_02_2024_23-58-12-000_256.csv.ms' %root_data_LIV)
LIV_IW_DJ2 = read('%s/LIV_IW_DJ2_Apr_02_2024_23-58-12-000_256.csv.ms' %root_data_LIV)

# seis
LIV_ES_DH1 = read('%s/LIV_ES_DH1_Apr_02_2024_23-58-12-000_512.csv.ms' %root_data_LIV)
LIV_ES_DH2 = read('%s/LIV_ES_DH2_Apr_02_2024_23-58-12-000_512.csv.ms' %root_data_LIV)
LIV_ES_DH3 = read('%s/LIV_ES_DH3_Apr_02_2024_23-58-12-000_512.csv.ms' %root_data_LIV)
LIV_EW_DH1 = read('%s/LIV_EW_DH1_Apr_02_2024_23-58-12-000_512.csv.ms' %root_data_LIV)
LIV_EW_DH2 = read('%s/LIV_EW_DH2_Apr_02_2024_23-58-12-000_512.csv.ms' %root_data_LIV)
LIV_EW_DH3 = read('%s/LIV_EW_DH3_Apr_02_2024_23-58-12-000_512.csv.ms' %root_data_LIV)

LIV_IS_DH1 = read('%s/LIV_IS_DH1_Apr_02_2024_23-58-12-000_512.csv.ms' %root_data_LIV)
LIV_IS_DH2 = read('%s/LIV_IS_DH2_Apr_02_2024_23-58-12-000_512.csv.ms' %root_data_LIV)
LIV_IS_DH3 = read('%s/LIV_IS_DH3_Apr_02_2024_23-58-12-000_512.csv.ms' %root_data_LIV)
LIV_IW_DH1 = read('%s/LIV_IW_DH1_Apr_02_2024_23-58-12-000_512.csv.ms' %root_data_LIV)
LIV_IW_DH2 = read('%s/LIV_IW_DH2_Apr_02_2024_23-58-12-000_512.csv.ms' %root_data_LIV)
LIV_IW_DH3 = read('%s/LIV_IW_DH3_Apr_02_2024_23-58-12-000_512.csv.ms' %root_data_LIV)

fig, axs = plt.subplots(8,1, figsize=(11,16))
plt.subplots_adjust(hspace=0)
labels = ['LIV_ES_DJ1 - LIV_IS_DJ1', 'LIV_EW_DJ2 - LIV_IW_DJ2',
          'LIV_ES_DH1 - LIV_IS_DH1', 'LIV_ES_DH2 - LIV_IS_DH2', 'LIV_ES_DH3 - LIV_IS_DH3',
          'LIV_EW_DH1 - LIV_IW_DH1', 'LIV_EW_DH2 - LIV_IW_DH2', 'LIV_EW_DH3 - LIV_IW_DH3']

for data1, data2, ax, label in zip([LIV_ES_DJ1, LIV_EW_DJ2,
                                    LIV_ES_DH1, LIV_ES_DH2, LIV_ES_DH3,
                                    LIV_EW_DH1, LIV_EW_DH2, LIV_EW_DH3],
                                   [LIV_IS_DJ1, LIV_IW_DJ2,
                                    LIV_IS_DH1, LIV_IS_DH2, LIV_IS_DH3,
                                    LIV_IW_DH1, LIV_IW_DH2, LIV_IW_DH3], axs, labels):
    if 'DJ' in label:
        color='red'
    elif 'DH' in label:
        color='black'

    data1_proc = data1.copy().detrend('simple').filter('highpass', freq=0.008, zerophase=True).slice(starttime=start, endtime=end)
    data2_proc = data2.copy().detrend('simple').filter('highpass', freq=0.008, zerophase=True).slice(starttime=start, endtime=end)

    '''output = yaracorrelate(data1_proc[0].data, data2_proc[0].data, method='same')
    ax.plot(output, color=color, label=label)'''

    ax.axvline(0, color='grey')

    lag, corr, lag_part, corr_part = yaracorrelate2(data1_proc[0].data, data2_proc[0].data, method='same', sps=data2_proc[0].stats.sampling_rate, onlypos=False)
    ax.plot(lag, corr, color=color, label=label)
    ax.legend(loc='upper right')
    ax.margins(0)
    print('done with: %s' % label)

axs[0].set_title('LIV_April02_crosscorrelation, .detrend(simple).filter(highpass, freq=0.008, zerophase=True)')
axs[-1].set_xlabel('Lag [s]')
plt.savefig('%s/LIV_April02_correlation_Bodywave_all_LIGO_stations.png' % root_savefig, dpi=300, bbox_inches='tight')

plt.show()