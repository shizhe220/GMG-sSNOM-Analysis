"""
Combined 2x3 subplot grid of graphene_4x1_manual lp2's "clean" hankel+FFT overlay
(plot_rp_hankel_fft_lp2_clean.py) across an E_F scan, for direct side-by-side
comparison -- easier to judge which E_F fits best than flipping between 5 separate
full-size figures. E_F=0.64eV added (not in the original 5-value request) purely to
fill the 2x3 grid evenly.
"""
import os
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from snippet.nearfield import rp_schematic
from NearFieldOptics import Materials as M

plt.rcParams.update({
    'font.size': 11, 'font.family': 'Arial',
    'axes.linewidth': 1.3, 'xtick.major.width': 1.3, 'ytick.major.width': 1.3,
    'xtick.direction': 'in', 'ytick.direction': 'in', 'xtick.top': True,
    'ytick.right': True
})

SAVE_DIR = 'figures/rp_dispersion/graphene_4x1_manual_lp2_hankel_fft_compare_v1'
os.makedirs(SAVE_DIR, exist_ok=True)

phonon_a = [(972, 4, 820, 4)]
phonon_b = [(1004, 2, 958, 2)]
MoO3_a = M.AnisotropicMaterial(eps_infinity=[4, 5.2, 2.4], phonon_params=[phonon_a, phonon_a, phonon_b])
materials_4x1_a = {'SiO2': M.SiO2_300nm, 'MoO3': MoO3_a}
SUBS = M.Si
stack_4x1 = [('MoO3', 4.2), 'graphene', ('SiO2', 300)]
GAMMA = 20
EF_LIST = [0.60, 0.63, 0.64, 0.65, 0.67, 0.69]
FREQS = np.linspace(850, 1010, 300)
QS = np.linspace(1e4, 9e5, 400)

single_wn = [860, 870, 880, 890]
single_q =  [94900, 104000, 109000, 111000]

two_wn = [900,    911,    920,    930,    941,    950,    960,    970,    980,    991,    1000]
two_q =  [150000, 155000, 150000, 202000, 191000, 191000, 206000, 265000, 189000, 206000, 257000]

fft_wn =       [860,   870,    880,    890,    900,    911,    920,    930,    941,    950,    960,    970,    980,    991,    1000]
fft_2q2 =      [97310, 112616, 118051, 124067, 154612, 159336, 152389, 210067, 198871, 220156, 265793, 196091, 251742, 276041, 323365]
fft_2q2_fwhm = [112919, 68990,  87756,  88043,  132695, 42166,  47543,  57410,  55862,  73908,  50734,  23579,  63950,  74801,  91331]
# ^ corrected 2026-07-23 -- see plot_rp_hankel_fft_lp2_clean.py for the double-halving bug fix.

peak23_wn = [860,    870,    880,    890,    900,    911,    920,    930,    941,    950,    960,    970,    980,    991,    1000]
peak23_q =  [113400, 108900, 103300, 113400, 165700, 196700, 170400, 104700, 214200, 191600, 197700, 223000, 111900, 127100, 169500]

fig, axes = plt.subplots(2, 3, figsize=(16, 10), sharex=True, sharey=True)

for ax, EF_REF in zip(axes.flat, EF_LIST):
    rp_schematic(stack_4x1, materials_4x1_a, subs=SUBS, freqs=FREQS, qs=QS,
                 Ef=EF_REF, gamma=GAMMA, figAx=(fig, ax), vmin=0, vmax=8, show_Ef=False,
                 ytickspacing=20, fs=10)
    ax.plot(single_q, single_wn, 's', ms=6, color='#1c7293', mec='white', mew=0.8, zorder=6)
    ax.plot(two_q, two_wn, 's', ms=6, color='#d6604d', mec='white', mew=0.8, zorder=6)
    ax.errorbar(fft_2q2, fft_wn, xerr=[f / 2 for f in fft_2q2_fwhm],
                marker='v', color='#762a83', ms=6, mec='white', mew=0.8, ls='none',
                capsize=2, elinewidth=1.0, zorder=7)
    ax.plot(peak23_q, peak23_wn, 'D', ms=6, color='#4daf4a', mec='white', mew=0.8, zorder=7)
    ax.set_xlim(0, 4.5e5)
    ax.set_title(f'$E_F$ = {EF_REF:.2f} eV', fontsize=13, fontweight='bold')

for ax in axes[-1]:
    ax.set_xlabel(r'q ($10^5$ cm$^{-1}$)', fontweight='bold')
for ax in axes[:, 0]:
    ax.set_ylabel(r'$\omega$ (cm$^{-1}$)', fontweight='bold')

handles = [
    plt.Line2D([], [], marker='s', ls='none', ms=8, color='#1c7293', mec='white', mew=1, label='hankel (single-wave, 860-890)'),
    plt.Line2D([], [], marker='s', ls='none', ms=8, color='#d6604d', mec='white', mew=1, label='hankel (two-wave, 900-1000)'),
    plt.Line2D([], [], marker='v', ls='none', ms=8, color='#762a83', mec='white', mew=1, label='fft ("2q/2", last peak, halved)'),
    plt.Line2D([], [], marker='D', ls='none', ms=8, color='#4daf4a', mec='white', mew=1, label='real-space peak2-peak3 / 2'),
]
fig.legend(handles=handles, loc='lower center', ncol=4, frameon=False, fontsize=11, bbox_to_anchor=(0.5, -0.02))
fig.suptitle('graphene_4x1_manual lp2 (interior), MoO3(a), gamma=20cm-1 -- E_F scan\n'
             'hankel (single-wave 860-890 + two-wave 900-1000) + FFT (2q/2) + peak2-3 spacing',
             fontsize=14, fontweight='bold')
fig.tight_layout(rect=[0, 0.03, 1, 0.93])
out_path = f'{SAVE_DIR}/graphene_4x1_manual_lp2_hankel_fft_clean_Ef_grid.png'
fig.savefig(out_path, dpi=200, bbox_inches='tight')
plt.close(fig)
print(f'Saved {out_path}')
