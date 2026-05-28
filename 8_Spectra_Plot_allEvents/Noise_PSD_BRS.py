import os
from obspy import read, UTCDateTime, read_inventory
import pandas as pd
import numpy
from obspy.signal import PPSD
from obspy.imaging.cm import pqlx, viridis
import sys
from pathlib import Path

sys.path.append('/Users/yararossi/Documents/Work/Towards_Quantification/1_Coding/1.11_BRS_data_curation/BRS_DataCuration/2_analysis_codes/')
from functions import get_ligo_data, getBRS_SN
import matplotlib.pyplot as plt

# roots
root_savefig = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/2_Figures/'
root_data = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data/zenodo_upload/miniseed_LIGO'
root_xml = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data/metadata'
response_H = read_inventory('%s/LIGO_HAN_LIV.xml' %root_xml)
response_L = read_inventory('%s/LIGO_HAN_LIV_lowsps.xml' %root_xml)
path = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data'

# get RLNM
root_RLNM = '/Users/yararossi/Documents/Work/Random_Workstuff/RLNM/'
rlnm_df = pd.read_csv("%s/ModelForFelix_freq.csv" % root_RLNM)
rlnm_freq = rlnm_df['frequency'].values
rlnm_amp = numpy.sqrt(rlnm_df['rlnm_median (rad2/Hz)'].values*rlnm_freq)
rlnm_amp_psd = rlnm_df['rlnm_median (rad2/Hz)'].values

freq_rlnm = rlnm_freq
rlnm_db = 10 * numpy.log10(rlnm_amp_psd*1e18) # make nrad,

# get BRS theoretical noise level
freq_BRSSN, noise_BRSSN_asd = getBRS_SN()
noise_BRSSN_psd = (noise_BRSSN_asd) ** 2
noise_BRSSN_db = 10 * numpy.log10(noise_BRSSN_psd*1e18)

fig_liv_I_djx, ax_liv_I_djx = plt.subplots(figsize=(6, 4))
fig_liv_I_djy, ax_liv_I_djy = plt.subplots(figsize=(6, 4))
fig_liv_E_djx, ax_liv_E_djx = plt.subplots(figsize=(6, 4))
fig_liv_E_djy, ax_liv_E_djy = plt.subplots(figsize=(6, 4))
fig_han_E_djx, ax_han_E_djx = plt.subplots(figsize=(6, 4))
fig_han_E_djy, ax_han_E_djy = plt.subplots(figsize=(6, 4))


# import data:
# BRS

for folder_liv, folder_han, label in zip(['2025-06-20_lownoise_L','2025-02-16_highnoise_L'],
                                         ['2024-12-08_highnoise_han_L','2024-06-20_lownoise_han_L'],
                                         ['low noise','high noise']):

    print(folder_liv)
    LIV_EY_DJX, LIV_EX_DJY, LIV_IY_DJX, LIV_IX_DJY, LIV_IX_DH, LIV_IY_DH, LIV_EX_DH, LIV_EY_DH = get_ligo_data(folder_liv, location='LIV')
    streams = LIV_EY_DJX, LIV_EX_DJY, LIV_IY_DJX, LIV_IX_DJY, LIV_IX_DH, LIV_IY_DH, LIV_EX_DH, LIV_EY_DH
    name = ['LIV_EY_DJX', 'LIV_EX_DJY', 'LIV_IY_DJX', 'LIV_IX_DJY', 'LIV_IX_DH', 'LIV_IY_DH', 'LIV_EX_DH', 'LIV_EY_DH']
    nrows=8

    HAN_EY_DJX, HAN_EX_DJY, HAN_EY_DH, HAN_EX_DH, HAN_VA_DH = get_ligo_data(folder_han, location='HAN')
    streams = HAN_EY_DJX, HAN_EX_DJY, HAN_EY_DH, HAN_EX_DH, HAN_VA_DH
    name = ['HAN_EY_DJX', 'HAN_EX_DJY', 'HAN_EY_DH', 'HAN_EX_DH', 'HAN_VA_DH']
    nrows=5

    db_bins_r = (-55, 55, 1)
    db_bins_t = (-185, -75, 1)
    p_s_w_o = 1  # 0.4
    p_s_o = 0.03  # 0.007
    ppsd_length = 3600

    response = response_L



    r_ppsd_LIV_EY_DJX = PPSD(LIV_EY_DJX[0].stats, metadata=response, db_bins=db_bins_r,
                             period_step_octaves=p_s_o, period_smoothing_width_octaves=p_s_w_o,
                             ppsd_length=ppsd_length)
    r_ppsd_LIV_EX_DJY = PPSD(LIV_EX_DJY[0].stats, metadata=response, db_bins=db_bins_r,
                             period_step_octaves=p_s_o, period_smoothing_width_octaves=p_s_w_o,
                             ppsd_length=ppsd_length)
    r_ppsd_LIV_IY_DJX = PPSD(LIV_IY_DJX[0].stats, metadata=response, db_bins=db_bins_r,
                             period_step_octaves=p_s_o, period_smoothing_width_octaves=p_s_w_o,
                             ppsd_length=ppsd_length)
    r_ppsd_LIV_IX_DJY = PPSD(LIV_IX_DJY[0].stats, metadata=response, db_bins=db_bins_r,
                             period_step_octaves=p_s_o, period_smoothing_width_octaves=p_s_w_o,
                             ppsd_length=ppsd_length)


    r_ppsd_HAN_EY_DJX = PPSD(HAN_EY_DJX[0].stats, metadata=response, db_bins=db_bins_r,
                             period_step_octaves=p_s_o, period_smoothing_width_octaves=p_s_w_o,
                             ppsd_length=ppsd_length)
    r_ppsd_HAN_EX_DJY = PPSD(HAN_EX_DJY[0].stats, metadata=response, db_bins=db_bins_r,
                             period_step_octaves=p_s_o, period_smoothing_width_octaves=p_s_w_o,
                             ppsd_length=ppsd_length)

    r_ppsd_LIV_EY_DJX.add(LIV_EY_DJX[0])
    r_ppsd_LIV_EX_DJY.add(LIV_EX_DJY[0])
    r_ppsd_LIV_IY_DJX.add(LIV_IY_DJX[0])
    r_ppsd_LIV_IX_DJY.add(LIV_IX_DJY[0])

    r_ppsd_HAN_EY_DJX.add(HAN_EY_DJX[0])
    r_ppsd_HAN_EX_DJY.add(HAN_EX_DJY[0])

    cmap = pqlx
    freq = True
    fmin = 0.003      # because period_lim max = 6 s  -> min freq
    fmax = 5    # because period_lim min = 0.003 s -> max freq

    fig, ax = fig_liv_E_djx, ax_liv_E_djx
    ppsd = r_ppsd_LIV_EY_DJX
    period, p5 = ppsd.get_percentile(percentile=5)
    _, p50 = ppsd.get_percentile(percentile=50)
    _, p95 = ppsd.get_percentile(percentile=95)
    freq = 1.0 / period
    mask = (freq >= fmin) & (freq <= fmax)
    mask_brs = (freq_BRSSN >= fmin) & (freq_BRSSN <= fmax)
    mask_rlnm = (freq_rlnm >= fmin) & (freq_rlnm <= fmax)
    if 'low' in label:
        ax.semilogx(freq_rlnm[mask_rlnm], rlnm_db[mask_rlnm], color="grey", linewidth=2, label="RLNM", zorder=2)
        ax.semilogx(freq_BRSSN[mask_brs], noise_BRSSN_db[mask_brs], color="black", linewidth=2, linestyle='--', label="BRS noise Model", zorder=2)
    ax.fill_between(freq[mask], p5[mask,], p95[mask], alpha=0.25, label="5–95%", zorder=0)
    ax.semilogx(freq[mask], p50[mask], linewidth=2.0, label=f'median {label}', zorder=5)
    ax.set_xlabel("Frequency [Hz]")
    ax.set_ylabel("Amplitude [nrad$^2$/Hz] [dB]")  # keep your label
    ax.set_title('(c) LIV.EY.BJX', loc='left')
    ax.grid(True, which="both", linewidth=0.3)
    ax.legend(loc='lower right')

    fig, ax = fig_liv_E_djy, ax_liv_E_djy
    ppsd = r_ppsd_LIV_EX_DJY
    period, p5 = ppsd.get_percentile(percentile=5)
    _, p50 = ppsd.get_percentile(percentile=50)
    _, p95 = ppsd.get_percentile(percentile=95)
    freq = 1.0 / period
    mask = (freq >= fmin) & (freq <= fmax)
    mask_brs = (freq_BRSSN >= fmin) & (freq_BRSSN <= fmax)
    mask_rlnm = (freq_rlnm >= fmin) & (freq_rlnm <= fmax)
    if 'low' in label:
        ax.semilogx(freq_rlnm[mask_rlnm], rlnm_db[mask_rlnm], color="grey", linewidth=2, label="RLNM", zorder=2)
        ax.semilogx(freq_BRSSN[mask_brs], noise_BRSSN_db[mask_brs], color="black", linewidth=2, linestyle='--',label="BRS noise Model", zorder=2)
    ax.fill_between(freq[mask], p5[mask,], p95[mask], alpha=0.25, label="5–95%", zorder=0)
    ax.semilogx(freq[mask], p50[mask], linewidth=2.0, label=f'median {label}', zorder=5)
    ax.set_xlabel("Frequency [Hz]")
    ax.set_ylabel("Amplitude [nrad$^2$/Hz] [dB]")  # keep your label
    ax.set_title('(d) LIV.EX.BJY', loc='left')
    ax.grid(True, which="both", linewidth=0.3)
    ax.legend(loc='lower right')

    fig, ax = fig_liv_I_djx, ax_liv_I_djx
    ppsd = r_ppsd_LIV_IY_DJX
    period, p5 = ppsd.get_percentile(percentile=5)
    _, p50 = ppsd.get_percentile(percentile=50)
    _, p95 = ppsd.get_percentile(percentile=95)
    freq = 1.0 / period
    mask = (freq >= fmin) & (freq <= fmax)
    mask_brs = (freq_BRSSN >= fmin) & (freq_BRSSN <= fmax)
    mask_rlnm = (freq_rlnm >= fmin) & (freq_rlnm <= fmax)
    if 'low' in label:
        ax.semilogx(freq_rlnm[mask_rlnm], rlnm_db[mask_rlnm], color="grey", linewidth=2, label="RLNM", zorder=2)
        ax.semilogx(freq_BRSSN[mask_brs], noise_BRSSN_db[mask_brs], color="black", linewidth=2, linestyle='--', label="BRS noise Model", zorder=2)
    ax.fill_between(freq[mask], p5[mask,], p95[mask], alpha=0.25, label="5–95%", zorder=0)
    ax.semilogx(freq[mask], p50[mask], linewidth=2.0, label=f'median {label}', zorder=5)
    ax.set_xlabel("Frequency [Hz]")
    ax.set_ylabel("Amplitude [nrad$^2$/Hz] [dB]")  # keep your label
    ax.set_title('(e) LIV.IY.BJX', loc='left')
    ax.grid(True, which="both", linewidth=0.3)
    ax.legend(loc='lower right')

    fig, ax = fig_liv_I_djy, ax_liv_I_djy
    ppsd = r_ppsd_LIV_IX_DJY
    period, p5 = ppsd.get_percentile(percentile=5)
    _, p50 = ppsd.get_percentile(percentile=50)
    _, p95 = ppsd.get_percentile(percentile=95)
    freq = 1.0 / period
    mask = (freq >= fmin) & (freq <= fmax)
    mask_brs = (freq_BRSSN >= fmin) & (freq_BRSSN <= fmax)
    mask_rlnm = (freq_rlnm >= fmin) & (freq_rlnm <= fmax)
    if 'low' in label:
        ax.semilogx(freq_rlnm[mask_rlnm], rlnm_db[mask_rlnm], color="grey", linewidth=2, label="RLNM", zorder=2)
        ax.semilogx(freq_BRSSN[mask_brs], noise_BRSSN_db[mask_brs], color="black", linewidth=2, linestyle='--',label="BRS noise Model", zorder=2)
    ax.fill_between(freq[mask], p5[mask,], p95[mask], alpha=0.25, label="5–95%", zorder=0)
    ax.semilogx(freq[mask], p50[mask], linewidth=2.0, label=f'median {label}', zorder=5)
    ax.set_xlabel("Frequency [Hz]")
    ax.set_ylabel("Amplitude [nrad$^2$/Hz] [dB]")  # keep your label
    ax.set_title('(f) LIV.IX.BJY', loc='left')
    ax.grid(True, which="both", linewidth=0.3)
    ax.legend(loc='lower right')




    fig, ax = fig_han_E_djx, ax_han_E_djx
    ppsd = r_ppsd_HAN_EY_DJX
    period, p5 = ppsd.get_percentile(percentile=5)
    _, p50 = ppsd.get_percentile(percentile=50)
    _, p95 = ppsd.get_percentile(percentile=95)
    freq = 1.0 / period
    mask = (freq >= fmin) & (freq <= fmax)
    mask_brs = (freq_BRSSN >= fmin) & (freq_BRSSN <= fmax)
    mask_rlnm = (freq_rlnm >= fmin) & (freq_rlnm <= fmax)
    if 'low' in label:
        ax.semilogx(freq_rlnm[mask_rlnm], rlnm_db[mask_rlnm], color="grey", linewidth=2, label="RLNM", zorder=2)
        ax.semilogx(freq_BRSSN[mask_brs], noise_BRSSN_db[mask_brs], color="black", linewidth=2, linestyle='--',label="BRS noise Model", zorder=2)
    ax.fill_between(freq[mask], p5[mask,], p95[mask], alpha=0.25, label="5–95%", zorder=0)
    ax.semilogx(freq[mask], p50[mask], linewidth=2.0, label=f'median {label}', zorder=5)
    ax.set_xlabel("Frequency [Hz]")
    ax.set_ylabel("Amplitude [nrad$^2$/Hz] [dB]")  # keep your label
    ax.set_title('(a) HAN.EY.BJX', loc='left')
    ax.grid(True, which="both", linewidth=0.3)
    ax.legend(loc='lower right')

    fig, ax = fig_han_E_djy, ax_han_E_djy
    ppsd = r_ppsd_HAN_EX_DJY
    period, p5 = ppsd.get_percentile(percentile=5)
    _, p50 = ppsd.get_percentile(percentile=50)
    _, p95 = ppsd.get_percentile(percentile=95)
    freq = 1.0 / period
    mask = (freq >= fmin) & (freq <= fmax)
    mask_brs = (freq_BRSSN >= fmin) & (freq_BRSSN <= fmax)
    mask_rlnm = (freq_rlnm >= fmin) & (freq_rlnm <= fmax)
    if 'low' in label:
        ax.semilogx(freq_rlnm[mask_rlnm], rlnm_db[mask_rlnm], color="grey", linewidth=2, label="RLNM", zorder=2)
        ax.semilogx(freq_BRSSN[mask_brs], noise_BRSSN_db[mask_brs], color="black", linewidth=2, linestyle='--', label="BRS noise Model", zorder=2)
    ax.fill_between(freq[mask], p5[mask,], p95[mask], alpha=0.25, label="5–95%", zorder=0)
    ax.semilogx(freq[mask], p50[mask], linewidth=2.0, label=f'median {label}', zorder=5)
    ax.set_xlabel("Frequency [Hz]")
    ax.set_ylabel("Amplitude [nrad$^2$/Hz] [dB]")  # keep your label
    ax.set_title('(b) HAN.EX.BJY', loc='left')
    ax.grid(True, which="both", linewidth=0.3)
    ax.legend(loc='lower right')



fig_liv_E_djx.savefig('%s/PPSD/Noise_PPSDper_LIV_EY_DJX.png' % (root_savefig), dpi=300, bbox_inches='tight')
fig_liv_E_djy.savefig('%s/PPSD/Noise_PPSDper_LIV_EX_DJY.png' % (root_savefig), dpi=300, bbox_inches='tight')

fig_liv_I_djx.savefig('%s/PPSD/Noise_PPSDper_LIV_IY_DJX.png' % (root_savefig), dpi=300, bbox_inches='tight')
fig_liv_I_djy.savefig('%s/PPSD/Noise_PPSDper_LIV_IX_DJY.png' % (root_savefig), dpi=300, bbox_inches='tight')

fig_han_E_djx.savefig('%s/PPSD/Noise_PPSDper_HAN_EY_DJX.png' % (root_savefig), dpi=300, bbox_inches='tight')
fig_han_E_djy.savefig('%s/PPSD/Noise_PPSDper_HAN_EX_DJY.png' % (root_savefig), dpi=300, bbox_inches='tight')

plt.show()