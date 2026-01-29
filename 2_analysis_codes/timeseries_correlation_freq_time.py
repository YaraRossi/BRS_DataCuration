import os

import matplotlib.pyplot as plt
from obspy import read, UTCDateTime, read_inventory
from correlation_function import correlate2

root_data = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data/miniseed_LIGO'
root_savefig = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/2_Figures/LIGO_analysis'
root_xml = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data/metadata'
response = read_inventory('%s/LIGO_HAN_LIV.xml' %root_xml)

for folder in os.listdir(root_data):
    if folder in ['.DS_Store', '2023-08-20', '2025-02-15']:
        continue
    #if '2024-09-15' not in folder:
    #   continue

    # make a folder to save all the figures
    newpath = r'%s/%s' % (root_savefig, folder)
    if not os.path.exists(newpath):
        os.makedirs(newpath)

    LIV, HAN = True, False

    if folder in ['2024-07-31','2024-09-23','2024-05-28','2024-07-02','2025-01-20','2023-08-20','2022-01-15']:
        LIV = False
    if folder in ['2024-05-13']:
        HAN = False
    print(folder)

    if LIV:
        # BRS
        LIV_EY_DJX = read('%s/%s/LI_LIV_EY_DJX.ms' %(root_data,folder))
        LIV_EX_DJY = read('%s/%s/LI_LIV_EX_DJY.ms' %(root_data,folder))
        LIV_IY_DJX = read('%s/%s/LI_LIV_IY_DJX.ms' %(root_data,folder))
        LIV_IX_DJY = read('%s/%s/LI_LIV_IX_DJY.ms' %(root_data,folder))

        # seis
        LIV_EY_DHX = read('%s/%s/LI_LIV_EY_DHX.ms' %(root_data,folder))
        LIV_EY_DHY = read('%s/%s/LI_LIV_EY_DHY.ms' %(root_data,folder))
        LIV_EY_DHZ = read('%s/%s/LI_LIV_EY_DHZ.ms' %(root_data,folder))
        LIV_EX_DHX = read('%s/%s/LI_LIV_EX_DHX.ms' %(root_data,folder))
        LIV_EX_DHY = read('%s/%s/LI_LIV_EX_DHY.ms' %(root_data,folder))
        LIV_EX_DHZ = read('%s/%s/LI_LIV_EX_DHZ.ms' %(root_data,folder))

        LIV_IY_DHX = read('%s/%s/LI_LIV_IY_DHX.ms' %(root_data,folder))
        LIV_IY_DHY = read('%s/%s/LI_LIV_IY_DHY.ms' %(root_data,folder))
        LIV_IY_DHZ = read('%s/%s/LI_LIV_IY_DHZ.ms' %(root_data,folder))
        LIV_IX_DHX = read('%s/%s/LI_LIV_IX_DHX.ms' %(root_data,folder))
        LIV_IX_DHY = read('%s/%s/LI_LIV_IX_DHY.ms' %(root_data,folder))
        LIV_IX_DHZ = read('%s/%s/LI_LIV_IX_DHZ.ms' %(root_data,folder))

        fig_Ey, axs_Ey = plt.subplots(11, 1, figsize=(11,8), sharex=True)
        plt.subplots_adjust(hspace=0, top=0.95, bottom=0.05)
        fig_Ex, axs_Ex = plt.subplots(11, 1, figsize=(11,8), sharex=True)
        plt.subplots_adjust(hspace=0, top=0.95, bottom=0.05)
        fig_Iy, axs_Iy = plt.subplots(11, 1, figsize=(11,8), sharex=True)
        plt.subplots_adjust(hspace=0, top=0.95, bottom=0.05)
        fig_Ix, axs_Ix = plt.subplots(11, 1, figsize=(11,8), sharex=True)
        plt.subplots_adjust(hspace=0, top=0.95, bottom=0.05)

        for axs in [axs_Ey, axs_Ex, axs_Iy, axs_Ix]:
            for ax in axs:
                ax.tick_params(labelleft=False, labelright=False)
                ax.tick_params(left=False, right=False)

        for freqmin, freqmax, ax_ix, ax_iy, ax_ex, ax_ey in zip([0.008, 0.009, 0.01, 0.02, 0.03, 0.05, 0.06, 0.07, 0.08,0.09,0.008],
                                                                [0.009,  0.01, 0.02, 0.03,  0.04, 0.06, 0.07, 0.08, 0.09, 0.1, 49],
                                                                axs_Ix.flat, axs_Iy.flat, axs_Ex.flat, axs_Ey.flat):

            for data_t, label_t, ax in zip([LIV_IX_DHX, LIV_IY_DHY, LIV_EX_DHX, LIV_EY_DHY],
                                           ['LIV_IX_DHX', 'LIV_IY_DHY', 'LIV_EX_DHX', 'LIV_EY_DHY'],
                                           [ax_ix, ax_iy, ax_ex, ax_ey]):
                data_t = data_t.filter('bandpass', freqmin=0.008, freqmax=49, zerophase=True).resample(256)
                data_t = data_t.remove_response(inventory=response, pre_filt=None, output="DEF", water_level=60, plot=False)
                data_proc = data_t.copy().detrend('simple').taper(0.1).filter('bandpass', freqmin=freqmin, freqmax=freqmax,zerophase=True)
                ax.plot(data_proc[0].times('matplotlib'), data_proc[0].data, color='k', label=label_t)

            for data_r, label_r, ax in zip([ LIV_IX_DJY, LIV_IY_DJX, LIV_EX_DJY, LIV_EY_DJX],
                                       ['LIV_IX_DJY','LIV_IY_DJX','LIV_EX_DJY','LIV_EY_DJX'],
                                       [ax_ix, ax_iy, ax_ex, ax_ey]):
                ax_twin = ax.twinx()
                data_r = data_r.filter('bandpass', freqmin=0.001, freqmax=100, zerophase=True)
                data_r = data_r.remove_response(inventory=response, pre_filt=None, output="DEF", water_level=60, plot=False)
                data_proc = data_r.copy().detrend('simple').taper(0.1).filter('bandpass', freqmin=freqmin, freqmax=freqmax,zerophase=True)
                ax_twin.plot(data_proc[0].times('matplotlib'), data_proc[0].data, color='red', label=label_r)

                #ax.legend(loc='upper right')
                ax.margins(0)
                ax.xaxis_date()
                ax.text(0.80,0.1, f'{freqmin}-{freqmax} Hz', transform=ax_twin.transAxes)
            start = data_proc[0].stats.starttime
        for ax, title in zip([axs_Ix.flat, axs_Iy.flat, axs_Ex.flat, axs_Ey.flat],
                             ['Ix', 'Iy', 'Ex', 'Ey']):
            ax[0].set_title(f'LIV_{folder}_{title}')
            ax[-1].set_xlabel(f'Time since {start} UTC [s]')
        for fig, name in zip([fig_Ey, fig_Ex, fig_Iy, fig_Ix],['Ey', 'Ex', 'IY', 'Ix']):
            #fig.savefig(f'{root_savefig}/{folder}/LI_LIV_{folder}_{name}_dispersion_all_stations.png', dpi=150, bbox_inches='tight')
            fig.savefig('%s/%s/LI_LIV_%s_%s_dispersion_all_stations.png' %(root_savefig,folder,folder,name))#, dpi=10,bbox_inches='tight')
        plt.show()
        #exit()

    if HAN:
        # BRS
        HAN_EY_DJX = read('%s/%s/LI_HAN_EY_DJX.ms' % (root_data, folder))
        HAN_EX_DJY = read('%s/%s/LI_HAN_EX_DJY.ms' % (root_data, folder))

        # seis
        HAN_EY_DHX = read('%s/%s/LI_HAN_EY_DHX.ms' % (root_data, folder))
        HAN_EY_DHY = read('%s/%s/LI_HAN_EY_DHY.ms' % (root_data, folder))
        HAN_EY_DHZ = read('%s/%s/LI_HAN_EY_DHZ.ms' % (root_data, folder))
        HAN_EX_DHX = read('%s/%s/LI_HAN_EX_DHX.ms' % (root_data, folder))
        HAN_EX_DHY = read('%s/%s/LI_HAN_EX_DHY.ms' % (root_data, folder))
        HAN_EX_DHZ = read('%s/%s/LI_HAN_EX_DHZ.ms' % (root_data, folder))

        HAN_IY_DHX = read('%s/%s/LI_HAN_IY_DHX.ms' % (root_data, folder))
        HAN_IY_DHY = read('%s/%s/LI_HAN_IY_DHY.ms' % (root_data, folder))
        HAN_IY_DHZ = read('%s/%s/LI_HAN_IY_DHZ.ms' % (root_data, folder))
        HAN_VA_DHX = read('%s/%s/LI_HAN_VA_DHX.ms' % (root_data, folder))
        HAN_VA_DHY = read('%s/%s/LI_HAN_VA_DHY.ms' % (root_data, folder))
        #HAN_IX_DHZ = read('%s/%s/LI_HAN_IX_DHZ.ms' % (root_data, folder))

        channels = [HAN_EY_DJX, HAN_EY_DHX, HAN_EY_DHY, HAN_EY_DHZ,
                                    HAN_EX_DJY, HAN_EX_DHX, HAN_EX_DHY, HAN_EX_DHZ,
                                     HAN_IY_DHX, HAN_IY_DHY, HAN_IY_DHZ,
                                     HAN_VA_DHX, HAN_VA_DHY #, HAN_IX_DHZ
                    ]
        for freqmin, freqmax in zip([0.008, 0.008, 0.01, 0.1, 1, 10], [49, 0.05, 0.1, 1, 10, 49]):
            fig, axs = plt.subplots(len(channels), 1, figsize=(11, len(channels)), sharex=True)
            plt.subplots_adjust(hspace=0, top=0.98, bottom=0.02)
            labels = ['HAN_EY_DJX', 'HAN_EY_DHX', 'HAN_EY_DHY', 'HAN_EY_DHZ',
                      'HAN_EX_DJY', 'HAN_EX_DHX', 'HAN_EX_DHY', 'HAN_EX_DHZ',
                      'HAN_IY_DHX', 'HAN_IY_DHY', 'HAN_IY_DHZ',
                      'HAN_VA_DHX', 'HAN_VA_DHY' #, 'HAN_IX_DHZ'
                      ]

            for data, ax, label in zip(channels, axs, labels):
                if 'DJ' in label:
                    color = 'red'
                    ax.set_ylabel('[nrad]')
                elif 'DH' in label:
                    color = 'black'
                    ax.set_ylabel('[nm/s]')

                data_proc = data.copy().detrend('simple').filter('bandpass',freqmin=freqmin,freqmax=freqmax, zerophase=True).taper(0.01)
                ax.plot(data_proc[0].times('matplotlib'), data_proc[0].data, color=color, label=label)
                ax.legend(loc='upper right')
                ax.margins(0)
                ax.xaxis_date()
            start = data_proc[0].stats.starttime
            axs[0].set_title('HAN_%s_all_stations: bandpass, freq=%s - %s Hz' % (folder,freqmin,freqmax))
            axs[-1].set_xlabel('Time since %s UTC [s]' % start)
            plt.savefig('%s/%s/LI_HAN_%s_freq%s_%s_all_stations.png' % (root_savefig, folder, folder,freqmin,freqmax), dpi=300, bbox_inches='tight')
            #plt.show()

