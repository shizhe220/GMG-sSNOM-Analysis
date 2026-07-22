"""
rp(q,w) overlay combining all 3 independent q_p cross-checks for graphene_4x1_manual
lp1 (edge), wn=860-1000cm-1: hankel (single-wave 860-890, two-wave 900-1000), FFT
(single-peak fit 860-890, 2-peak fit 900-1000), and real-space peak2-peak3/2 spacing.

Source of truth for all values is
figures/rp_dispersion/graphene_4x1_manual_lp1_peak_spacing_crosscheck_v1/GMG#3_4x1um_GMGregion_leftedge_v2(07222026).csv
(see the sibling _README.txt in that folder for what each column means). The FFT series
plotted here always uses fft_recommended_qp_cm-1/_fwhm_cm-1 -- a per-row pre-selected
channel (amp for 860-980, phase for 991/1000, where amp is essentially unresolved) so
this script never has to hand-pick amp vs phase itself.
"""
import os
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from snippet.nearfield import rp_schematic
from NearFieldOptics import Materials as M

plt.rcParams.update({
    'font.size': 12, 'font.family': 'Arial',
    'axes.linewidth': 1.5, 'xtick.major.width': 1.5, 'ytick.major.width': 1.5,
    'xtick.direction': 'in', 'ytick.direction': 'in', 'xtick.top': True,
    'ytick.right': True
})

SAVE_DIR = 'figures/rp_dispersion/graphene_4x1_manual_lp1_peak_spacing_crosscheck_v1/two_wave_hankel_900to1000'
CSV_PATH = 'figures/rp_dispersion/graphene_4x1_manual_lp1_peak_spacing_crosscheck_v1/GMG#3_4x1um_GMGregion_leftedge_v2(07222026).csv'
os.makedirs(SAVE_DIR, exist_ok=True)

df = pd.read_csv(CSV_PATH)

# --- rp background ---
phonon_a = [(972, 4, 820, 4)]
phonon_b = [(1004, 2, 958, 2)]
MoO3_a = M.AnisotropicMaterial(eps_infinity=[4, 5.2, 2.4], phonon_params=[phonon_a, phonon_a, phonon_b])
materials_4x1_a = {'SiO2': M.SiO2_300nm, 'MoO3': MoO3_a}
SUBS = M.Si
stack_4x1 = [('MoO3', 4.2), 'graphene', ('SiO2', 300)]
GAMMA = 20
EF_REF = 0.69
FREQS = np.linspace(850, 1010, 300)
QS = np.linspace(1e4, 9e5, 400)

fig, ax = plt.subplots(figsize=(10, 9))
rp_schematic(stack_4x1, materials_4x1_a, subs=SUBS, freqs=FREQS, qs=QS,
             Ef=EF_REF, gamma=GAMMA, figAx=(fig, ax), vmin=0, vmax=8, show_Ef=True, ytickspacing=20)

# --- hankel: single-wave (860-890) vs two-wave (900-1000) ---
single_mask = df['freq_cm-1'] < 900
two_mask = df['freq_cm-1'] >= 900
ax.plot(df.loc[single_mask, 'hankel_qp_cm-1'], df.loc[single_mask, 'freq_cm-1'],
        's', ms=9, mfc='none', mec='#1c7293', mew=1.8, label='hankel (single-wave, 860-890)', zorder=6)
ax.plot(df.loc[two_mask, 'hankel_2wave_qp_cm-1'], df.loc[two_mask, 'freq_cm-1'],
        's', ms=9, color='#1c7293', mec='white', mew=1, label='hankel (two-wave, 900-1000)', zorder=6)

# --- real-space peak2-peak3/2 (860-1000, corrected) ---
ax.plot(df['peak2to3_spacing_qp_cm-1'], df['freq_cm-1'],
        'D', ms=8, color='#4daf4a', mec='white', mew=1, label='real-space peak2-peak3 / 2 (860-1000, corrected)', zorder=6)

# --- FFT: pre-selected "recommended" channel per row (amp for 860-980, phase for 991/1000) ---
# single-peak fit (860-890, only ever resolves "2q") vs 2-peak fit (900-1000, resolves q and 2q)
single_fft = df[(df['freq_cm-1'] < 900) & (df['fft_recommended_channel'] == 'amp')]
ax.errorbar(single_fft['fft_recommended_qp_cm-1'], single_fft['freq_cm-1'], xerr=single_fft['fft_recommended_fwhm_cm-1'] / 2,
            marker='v', color='#762a83', ms=8, mec='white', mew=1, ls='none', capsize=3, elinewidth=1.3,
            label='fft (single-peak fit, 860-890)', zorder=7)

amp_2peak = df[(df['freq_cm-1'] >= 900) & (df['fft_recommended_channel'] == 'amp')]
ax.errorbar(amp_2peak['fft_recommended_qp_cm-1'], amp_2peak['freq_cm-1'], xerr=amp_2peak['fft_recommended_fwhm_cm-1'] / 2,
            marker='v', color='#762a83', ms=8, mec='white', mew=1, ls='none', capsize=3, elinewidth=1.3,
            label='fft (2-peak fit, amp ch., 900-980)', zorder=7)

phase_rows = df[df['fft_recommended_channel'] == 'phase']
ax.errorbar(phase_rows['fft_recommended_qp_cm-1'], phase_rows['freq_cm-1'], xerr=phase_rows['fft_recommended_fwhm_cm-1'] / 2,
            marker='v', color='#762a83', ms=8, mec='#f7d94c', mew=2, ls='none', capsize=3, elinewidth=1.3,
            label='fft (2-peak fit, phase ch., 991/1000 -- amp too noisy)', zorder=8)

ax.set_xlim(0, 4.5e5)
ax.set_title(f'graphene_4x1_manual lp1 (edge), MoO3(a), E_F={EF_REF}eV (fixed ref), gamma={GAMMA}cm-1\n'
             f'wn=860-1000cm-1: hankel + FFT + real-space peak2-3 spacing', fontsize=12, fontweight='bold')
ax.legend(frameon=False, fontsize=9, loc='upper left')
fig.tight_layout()
out_path = f'{SAVE_DIR}/rp_hankel_fft_peak23_860to1000_Ef{EF_REF}eV.png'
fig.savefig(out_path, dpi=200, bbox_inches='tight')
plt.close(fig)
print(f'Saved {out_path}')
