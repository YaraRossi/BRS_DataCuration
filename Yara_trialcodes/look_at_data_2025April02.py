import matplotlib.pyplot as plt
import obspy
from obspy import read, UTCDateTime

root_data_HAN = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data/HAN_APR_02'
root_data_LIV = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data/LIV_APR_02'
root_savefig = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/2_Figures/2024Apr02'

LIV = False
HAN = True

start = UTCDateTime('2024-0403T00:15:00')
end = UTCDateTime('2024-0403T02:15:00')

if LIV:
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



    fig, axs = plt.subplots(16,1, figsize=(11,16), sharex=True)
    plt.subplots_adjust(hspace=0, top=0.98, bottom=0.02)
    labels = ['LIV_ES_DJ1', 'LIV_ES_DH1', 'LIV_ES_DH2', 'LIV_ES_DH3',
              'LIV_EW_DJ2', 'LIV_EW_DH1', 'LIV_EW_DH2', 'LIV_EW_DH3',
              'LIV_IS_DJ1', 'LIV_IS_DH1', 'LIV_IS_DH2', 'LIV_IS_DH3',
              'LIV_IW_DJ2', 'LIV_IW_DH1', 'LIV_IW_DH2', 'LIV_IW_DH3']

    for data, ax, label in zip([LIV_ES_DJ1, LIV_ES_DH1, LIV_ES_DH2, LIV_ES_DH3,
                                LIV_EW_DJ2, LIV_EW_DH1, LIV_EW_DH2, LIV_EW_DH3,
                                LIV_IS_DJ1, LIV_IS_DH1, LIV_IS_DH2, LIV_IS_DH3,
                                LIV_IW_DJ2, LIV_IW_DH1, LIV_IW_DH2, LIV_IW_DH3], axs, labels):
        if 'DJ' in label:
            color='red'
        elif 'DH' in label:
            color='black'

        data_proc = data.copy().detrend('simple').filter('highpass', freq=0.008, zerophase=True).slice(starttime=start, endtime=end)
        ax.plot(data_proc[0].times('matplotlib'), data_proc[0].data, color=color, label=label)
        ax.legend(loc='upper right')
        ax.margins(0)
        ax.xaxis_date()

    axs[0].set_title('LIV_April02_all_stations, .detrend(simple).filter(highpass, freq=0.008, zerophase=True)')
    axs[-1].set_xlabel('Time since %s UTC [s]' %start)
    plt.savefig('%s/LIV_April02_all_stations.png' % root_savefig, dpi=300, bbox_inches='tight')
    plt.show()

########################################################################################################################
################################################ LIGO Hanford Site #####################################################
########################################################################################################################
if HAN:
    # BRS
    HAN_EN_DJ1 = read('%s/HAN_EN_DJ1_Apr_02_2024_23-58-12-000_256.csv.ms' %root_data_HAN)
    HAN_ES_DJ2 = read('%s/HAN_ES_DJ2_Apr_02_2024_23-58-12-000_256.csv.ms' %root_data_HAN)
    # seis
    HAN_EN_DH1 = read('%s/HAN_EN_DH1_Apr_02_2024_23-58-12-000_512.csv.ms' %root_data_HAN)
    HAN_EN_DH2 = read('%s/HAN_EN_DH2_Apr_02_2024_23-58-12-000_512.csv.ms' %root_data_HAN)
    HAN_EN_DH3 = read('%s/HAN_EN_DH3_Apr_02_2024_23-58-12-000_512.csv.ms' %root_data_HAN)
    HAN_ES_DH1 = read('%s/HAN_ES_DH1_Apr_02_2024_23-58-12-000_512.csv.ms' %root_data_HAN)
    HAN_ES_DH2 = read('%s/HAN_ES_DH2_Apr_02_2024_23-58-12-000_512.csv.ms' %root_data_HAN)
    HAN_ES_DH3 = read('%s/HAN_ES_DH3_Apr_02_2024_23-58-12-000_512.csv.ms' %root_data_HAN)
    HAN_IS_DH1 = read('%s/HAN_IS_DH1_Apr_02_2024_23-58-12-000_512.csv.ms' %root_data_HAN)
    HAN_IS_DH2 = read('%s/HAN_IS_DH2_Apr_02_2024_23-58-12-000_512.csv.ms' %root_data_HAN)
    HAN_IS_DH3 = read('%s/HAN_IS_DH3_Apr_02_2024_23-58-12-000_512.csv.ms' %root_data_HAN)
    # vault
    HAN_VA_HHE = read('%s/LIGO_VA_HHE_Apr_02_2024_23-58-12-000_256.csv.ms' %root_data_HAN)
    HAN_VA_HHN = read('%s/LIGO_VA_HHN_Apr_02_2024_23-58-12-000_256.csv.ms' %root_data_HAN)
    # NIKE
    NIKE_UW_ENE = read('%s/NIKE.UW..ENE.2024.093.23.58.12.008-2024.094.02.58.11.998.rm.miniseed' %root_data_HAN)
    NIKE_UW_ENN = read('%s/NIKE.UW..ENN.2024.093.23.58.12.008-2024.094.02.58.11.998.rm.miniseed' %root_data_HAN)
    NIKE_UW_ENZ = read('%s/NIKE.UW..ENZ.2024.093.23.58.12.008-2024.094.02.58.11.998.rm.miniseed' %root_data_HAN)
    NIKE_UW_HHE = read('%s/NIKE.UW..HHE.2024.093.23.58.12.008-2024.094.02.58.11.998.rm.miniseed' %root_data_HAN)
    NIKE_UW_HHN = read('%s/NIKE.UW..HHN.2024.093.23.58.12.008-2024.094.02.58.11.998.rm.miniseed' %root_data_HAN)
    NIKE_UW_HHZ = read('%s/NIKE.UW..HHZ.2024.093.23.58.12.008-2024.094.02.58.11.998.rm.miniseed' %root_data_HAN)


    fig, axs = plt.subplots(19,1, figsize=(11,19), sharex=True)
    plt.subplots_adjust(hspace=0, top=0.98, bottom=0.02)
    labels = ['HAN_EN_DJ1', 'HAN_EN_DH1', 'HAN_EN_DH2', 'HAN_EN_DH3',
              'HAN_ES_DJ2', 'HAN_ES_DH1', 'HAN_ES_DH2', 'HAN_ES_DH3',
                            'HAN_IS_DH1', 'HAN_IS_DH2', 'HAN_IS_DH3',
                            'HAN_VA_HHE', 'HAN_VA_HHN',
                            'NIKE_UW_ENE','NIKE_UW_ENN','NIKE_UW_ENZ',
                            'NIKE_UW_HHE','NIKE_UW_HHN','NIKE_UW_HHZ'
                            ]

    for data, ax, label in zip([HAN_EN_DJ1, HAN_EN_DH1, HAN_EN_DH2, HAN_EN_DH3,
                                HAN_ES_DJ2, HAN_ES_DH1, HAN_ES_DH2, HAN_ES_DH3,
                                            HAN_IS_DH1, HAN_IS_DH2, HAN_IS_DH3,
                                            HAN_VA_HHE, HAN_VA_HHN,
                                            NIKE_UW_ENE,NIKE_UW_ENN,NIKE_UW_ENZ,
                                            NIKE_UW_HHE,NIKE_UW_HHN,NIKE_UW_HHZ], axs, labels):
        if 'DJ' in label:
            color='red'
        elif 'DH' in label:
            color='black'
        elif 'VA' in label:
            color='black'
        elif 'UW' in label:
            color='grey'

        data_proc = data.copy().detrend('simple').filter('highpass', freq=0.008, zerophase=True).slice(starttime=start, endtime=end)
        ax.plot(data_proc[0].times(), data_proc[0].data, color=color, label=label)
        ax.legend(loc='upper right')
        ax.margins(0)

    axs[0].set_title('HAN_April02_all_stations, .detrend(simple).filter(highpass, freq=0.008, zerophase=True)')
    axs[-1].set_xlabel('Time since %s UTC [s]' %start)
    plt.savefig('%s/HAN_April02_all_stations.png' % root_savefig, dpi=300, bbox_inches='tight')
    plt.show()





'''print(HAN_EN_DJ1[0].stats)
HAN_EN_DJ1_proc = HAN_EN_DJ1.copy().detrend('simple').filter('highpass', freq=0.1, zerophase=True)
HAN_EN_DJ1_proc.plot(method='full')'''