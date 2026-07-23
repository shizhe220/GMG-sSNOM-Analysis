"""
rp(q,w) "clean" overlay for graphene_4x1_manual lp2 (interior) -- the working-answer
version, as opposed to plot_rp_hankel_fft_lp2_compare.py (kept as-is, a historical
record of the q-vs-2q/2 comparison that led here).

Decisions made after that comparison (2026-07-23):
- FFT: only "2q/2" (last peak, halved) is plotted -- the "q" (first peak, direct) series
  had much larger FWHM error bars and was judged too unreliable to keep using.
- Hankel: single-wave for 860-890 only (all genuine, no bound-saturation issues there),
  two-wave for 900-1000 (now all genuine after 920/960 were re-tuned). The 7 wn where
  single-wave hits fit_cavity_prefactor_compare's default lam_bound_factor=(0.7,1.2)
  bound (911,930,941,950,960,980,991 -- see log entry) are simply not shown at all here,
  since two-wave already covers that same 900-1000 range with a genuine result --
  no need to flag an artifact when a good alternative for the same point exists.

Added 2026-07-23 (same day, follow-up): real-space peak2-peak3 spacing, all 15 wn, from
peak23_diagnostic_graphene_4x1_lp2.py after the user's manual peak-position review
(several rounds -- multiple wn had default peak-finding pick a valley or a shoulder
instead of the genuine local max; see that script's MANUAL_PEAK_XNM and the log entry
for the full correction history).
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
    'font.size': 12, 'font.family': 'Arial',
    'axes.linewidth': 1.5, 'xtick.major.width': 1.5, 'ytick.major.width': 1.5,
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
EF_LIST = [0.60, 0.63, 0.65, 0.67, 0.69]  # E_F scan, per user request 2026-07-23
FREQS = np.linspace(850, 1010, 300)
QS = np.linspace(1e4, 9e5, 400)

# hankel single-wave, 860-890 (genuine, no bound-saturation there)
single_wn = [860, 870, 880, 890]
single_q =  [94900, 104000, 109000, 111000]

# hankel two-wave "q" term, 900-1000 (920/960 re-tuned 2026-07-23, all genuine now)
two_wn = [900,    911,    920,    930,    941,    950,    960,    970,    980,    991,    1000]
two_q =  [150000, 155000, 150000, 202000, 191000, 191000, 206000, 265000, 189000, 206000, 257000]

# fft "2q/2" only, all 15 wn
fft_wn =       [860,   870,    880,    890,    900,    911,    920,    930,    941,    950,    960,    970,    980,    991,    1000]
fft_2q2 =      [97310, 112616, 118051, 124067, 154612, 159336, 152389, 210067, 198871, 220156, 265793, 196091, 251742, 276041, 323365]
fft_2q2_fwhm = [112919, 68990,  87756,  88043,  132695, 42166,  47543,  57410,  55862,  73908,  50734,  23579,  63950,  74801,  91331]
# ^ corrected 2026-07-23: earlier values here were exactly half the true FWHM (verified
# against the freshly-computed shared CSV, scripts/export_shared_csv_graphene_4x1_lp2.py --
# an ad hoc script earlier this session applied an extra, erroneous /2 on top of the
# errorbar's own xerr=fwhm/2 halving below, quietly halving the error bars again).

# real-space peak2-peak3 spacing / 2, all 15 wn (from peak23_diagnostic_graphene_4x1_lp2.py,
# after the user's manual peak-position review -- see that script for method/overrides)
peak23_wn = [860,    870,    880,    890,    900,    911,    920,    930,    941,    950,    960,    970,    980,    991,    1000]
peak23_q =  [113400, 108900, 103300, 113400, 165700, 196700, 170400, 104700, 214200, 191600, 197700, 223000, 111900, 127100, 169500]

for EF_REF in EF_LIST:
    fig, ax = plt.subplots(figsize=(9, 8))
    rp_schematic(stack_4x1, materials_4x1_a, subs=SUBS, freqs=FREQS, qs=QS,
                 Ef=EF_REF, gamma=GAMMA, figAx=(fig, ax), vmin=0, vmax=8, show_Ef=True, ytickspacing=20)

    ax.plot(single_q, single_wn, 's', ms=9, color='#1c7293', mec='white', mew=1,
            label='hankel (single-wave, 860-890)', zorder=6)
    ax.plot(two_q, two_wn, 's', ms=9, color='#d6604d', mec='white', mew=1,
            label='hankel (two-wave, 900-1000)', zorder=6)
    ax.errorbar(fft_2q2, fft_wn, xerr=[f / 2 for f in fft_2q2_fwhm],
                marker='v', color='#762a83', ms=8, mec='white', mew=1, ls='none', capsize=3, elinewidth=1.3,
                label='fft ("2q/2", last peak, halved)', zorder=7)
    ax.plot(peak23_q, peak23_wn, 'D', ms=8, color='#4daf4a', mec='white', mew=1,
            label='real-space peak2-peak3 / 2', zorder=7)

    ax.set_xlim(0, 4.5e5)
    ax.set_title(f'graphene_4x1_manual lp2 (interior), MoO3(a), E_F={EF_REF:.2f}eV, gamma={GAMMA}cm-1\n'
                 f'wn=860-1000cm-1: hankel + FFT (2q/2 only) + peak2-3 spacing', fontsize=12, fontweight='bold')
    ax.legend(frameon=False, fontsize=9.5, loc='upper left')
    fig.tight_layout()
    out_path = f'{SAVE_DIR}/graphene_4x1_manual_lp2_hankel_fft_clean_Ef{EF_REF:.2f}eV.png'
    fig.savefig(out_path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f'Saved {out_path}')
