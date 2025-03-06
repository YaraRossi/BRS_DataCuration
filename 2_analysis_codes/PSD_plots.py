import os
from scipy import signal
import matplotlib.pyplot as plt
from obspy import read, UTCDateTime
from correlation_function import correlate2

root_data = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data/miniseed_LIGO'
root_savefig = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/2_Figures/LIGO_analysis'

for folder in os.listdir(root_data):
    if folder == '.DS_Store':
        continue

    # make a folder to save all the figures
    newpath = r'%s/%s' % (root_savefig, folder)
    if not os.path.exists(newpath):
        os.makedirs(newpath)

    LIV, HAN = True, True

    if folder in ['July31.2024','Sept23.2024','May28.2024','July02.2024','Jan20.2025']:
        LIV = False
    if folder in ['May13.2024']:
        HAN = False

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

        fig, axs = plt.subplots(5,1, figsize=(7,8), sharex=True)
        plt.subplots_adjust(hspace=0) #, top=0.98, bottom=0.02)
        plt.grid(True)
        labels = ['LIV_EY_DJX', 'LIV_EY_DHX', 'LIV_EY_DHY', 'LIV_EY_DHZ',
                  'LIV_EX_DJY', 'LIV_EX_DHX', 'LIV_EX_DHY', 'LIV_EX_DHZ',
                  'LIV_IY_DJX', 'LIV_IY_DHX', 'LIV_IY_DHY', 'LIV_IY_DHZ',
                  'LIV_IX_DJY', 'LIV_IX_DHX', 'LIV_IX_DHY', 'LIV_IX_DHZ']

        ax_jx, ax_jy, ax_hx, ax_hy, ax_hz = axs[0], axs[1], axs[2], axs[3], axs[4]
        axs_plot = [ax_jx, ax_hx, ax_hy, ax_hz,
                    ax_jy, ax_hx, ax_hy, ax_hz,
                    ax_jx, ax_hx, ax_hy, ax_hz,
                    ax_jy, ax_hx, ax_hy, ax_hz,]

        for data, ax, label in zip([LIV_EY_DJX, LIV_EY_DHX, LIV_EY_DHY, LIV_EY_DHZ,
                                    LIV_EX_DJY, LIV_EX_DHX, LIV_EX_DHY, LIV_EX_DHZ,
                                    LIV_IY_DJX, LIV_IY_DHX, LIV_IY_DHY, LIV_IY_DHZ,
                                    LIV_IX_DJY, LIV_IX_DHX, LIV_IX_DHY, LIV_IX_DHZ], axs_plot, labels):
            if 'DJ' in label:
                ax.set_ylabel('[nrad]')
            elif 'DH' in label:
                ax.set_ylabel('[nm/s]')

            if 'DJ' in label:
                if '_I' in label:
                    color='red'
                elif '_E' in label:
                    color='indianred'

            if 'DH' in label:
                if '_I' in label:
                    color='black'
                elif '_E' in label:
                    color='grey'
            if 'X_' in label:
                linestyle = '-'
            elif 'Y_' in label:
                linestyle = '--'


            length = int(len(data[0].data) / 5)
            f_t, P_t = signal.welch(data[0].data, 1 / data[0].stats.delta, nperseg=length)
            ax.loglog(f_t, P_t, color=color, label=label, linestyle=linestyle)
            ax.legend(loc='lower left')
            ax.grid(which='major', axis='both')
            ax.margins(0)
        start = data[0].stats.starttime
        axs[0].set_title('LIV_%s_PSD' %folder)
        axs[-1].set_xlabel('Frequency [Hz]')
        plt.savefig('%s/%s/LI_LIV_%s_PSD.png' %(root_savefig,folder,folder), dpi=300, bbox_inches='tight')
        #plt.show()
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

        fig, axs = plt.subplots(5, 1, figsize=(7, 8), sharex=True)
        plt.subplots_adjust(hspace=0)  # , top=0.98, bottom=0.02)
        plt.grid(True)
        labels = ['HAN_EY_DJX', 'HAN_EY_DHX', 'HAN_EY_DHY', 'HAN_EY_DHZ',
                  'HAN_EX_DJY', 'HAN_EX_DHX', 'HAN_EX_DHY', 'HAN_EX_DHZ',
                                'HAN_IY_DHX', 'HAN_IY_DHY', 'HAN_IY_DHZ',
                                'HAN_VA_DHX', 'HAN_VA_DHY' #, 'HAN_IX_DHZ'
                  ]
        ax_jx, ax_jy, ax_hx, ax_hy, ax_hz = axs[0], axs[1], axs[2], axs[3], axs[4]
        axs_plot = [ax_jx, ax_hx, ax_hy, ax_hz,
                    ax_jy, ax_hx, ax_hy, ax_hz,
                           ax_hx, ax_hy, ax_hz,
                           ax_hx, ax_hy, ]

        for data, ax, label in zip(channels, axs_plot, labels):
            if 'DJ' in label:
                ax.set_ylabel('[nrad]')
            elif 'DH' in label:
                ax.set_ylabel('[nm/s]')

            if 'DJ' in label:
                if '_I' in label:
                    color='red'
                elif '_E' in label:
                    color='indianred'

            if 'DH' in label:
                if '_I' in label:
                    color='black'
                elif '_E' in label:
                    color='grey'
            if 'X_' in label:
                linestyle = '-'
            elif 'Y_' in label:
                linestyle = '--'

            length = int(len(data[0].data) / 5)
            f_t, P_t = signal.welch(data[0].data, 1 / data[0].stats.delta, nperseg=length)
            ax.loglog(f_t, P_t, color=color, label=label, linestyle=linestyle)
            ax.legend(loc='lower left')
            ax.grid(which='major', axis='both')
            ax.margins(0)
        axs[0].set_title('HAN_%s_PSD' %folder)
        axs[-1].set_xlabel('Frequency [Hz]')
        plt.savefig('%s/%s/LI_HAN_%s_PSD.png' %(root_savefig,folder,folder), dpi=300, bbox_inches='tight')
        #plt.show()

        #exit()

