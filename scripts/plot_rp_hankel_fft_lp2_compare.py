"""
rp(q,w) comparison overlay for graphene_4x1_manual lp2 (interior) -- still in the
comparison/exploration stage (2026-07-23), so this deliberately shows BOTH conventions
side by side rather than pre-committing to one, unlike lp1's final combined plot:

- Hankel: single-wave (fits assuming the whole oscillation IS the tip round-trip "2q"
  wave -- i.e. its own q_p output already represents "2q/2" under that assumption,
  labeled as such here) vs. two-wave's own directly-fitted "q" (edge-launched) term,
  for the SAME wn (900-1000) so they can be compared head to head.
- FFT: "q" (first/smaller peak, direct, no /2) vs. "2q/2" (last/larger peak, halved) --
  both plotted with FWHM/2 error bars wherever 2 peaks were resolved, instead of
  collapsing to one "established" choice (lp1's "always 2q" convention doesn't
  automatically carry over here since lp2 hasn't been validated the same way yet).

No FFT/peak2-3 cross-check has been done for lp2 real-space peak spacing yet (separate,
later task per user). Values read from fitting_pipeline_graphene_4x1_lp2.ipynb's own
tuned params, re-run here with plot=False for full-precision numbers (see log entry).
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
EF_REF = 0.69
FREQS = np.linspace(850, 1010, 300)
QS = np.linspace(1e4, 9e5, 400)

# --- Hankel: single-wave ("2q/2" convention), all 15 wn; bound-saturated ones flagged ---
# 911 and 960 fixed 2026-07-23 (were bound-saturated, now genuine -- see log entry);
# 930,941,950,980,991 remain bound-saturated -- confirmed (wide lam_bound_factor search)
# that single-wave has no genuine local minimum near the two-wave/FFT consensus q there,
# i.e. not a fixable tuning issue, the model itself is inadequate at those 5 wn.
single_wn =    [860,   870,    880,    890,    900,    911,    920,    930,    941,    950,    960,    970,    980,    991,   1000]
single_q =     [94900, 104000, 109000, 111000, 122000, 148400, 146000, 131000, 131000, 187000, 169800, 179000, 131000, 175000, 288000]
single_boundsat = {930, 941, 950, 980, 991}

# --- Hankel: two-wave "q" term, 900-1000 (920/960 re-tuned 2026-07-23, no longer degenerate) ---
two_wn =     [900,    911,    920,    930,    941,    950,    960,    970,    980,    991,    1000]
two_q =      [150000, 155000, 150000, 202000, 191000, 191000, 206000, 265000, 189000, 206000, 257000]
two_flagged = set()

# --- FFT: "2q/2" (all 15) and "q direct" (911-1000, where 2 peaks resolved) ---
fft_wn =       [860,   870,    880,    890,    900,    911,    920,    930,    941,    950,    960,    970,    980,    991,    1000]
fft_2q2 =      [97310, 112616, 118051, 124067, 154612, 159336, 152389, 210067, 198871, 220156, 265793, 196091, 251742, 276041, 323365]
fft_2q2_fwhm = [112919, 68990,  87756,  88043,  132695, 42166,  47543,  57410,  55862,  73908,  50734,  23579,  63950,  74801,  91331]
fft_direct_wn =    [911,    920,    930,    941,    950,    960,    970,    980,    991,    1000]
fft_direct =       [140449, 116763, 210891, 172854, 173727, 262351, 271256, 190145, 242116, 273502]
fft_direct_fwhm =  [90859,  125804, 98883,  157792, 113258, 232240, 81862,  80340,  138089, 166433]
# ^ both fwhm arrays corrected 2026-07-23: earlier values here were exactly half the true
# FWHM in both conventions (verified against a from-scratch recompute, see
# scripts/export_shared_csv_graphene_4x1_lp2.py and plot_rp_hankel_fft_lp2_clean.py) --
# an ad hoc script earlier this session applied an extra, erroneous /2 on top of the
# errorbar's own xerr=fwhm/2 halving below, quietly halving the error bars again.

fig, ax = plt.subplots(figsize=(10, 9))
rp_schematic(stack_4x1, materials_4x1_a, subs=SUBS, freqs=FREQS, qs=QS,
             Ef=EF_REF, gamma=GAMMA, figAx=(fig, ax), vmin=0, vmax=8, show_Ef=True, ytickspacing=20)

# hankel single-wave (2q/2), flagged points open with x-hatch look via marker edge only
ok_wn = [wn for wn in single_wn if wn not in single_boundsat]
ok_q = [q for wn, q in zip(single_wn, single_q) if wn not in single_boundsat]
ax.plot(ok_q, ok_wn, 's', ms=8, color='#1c7293', mec='white', mew=1,
        label='hankel single-wave ("2q/2")', zorder=6)
bad_wn = [wn for wn in single_wn if wn in single_boundsat]
bad_q = [q for wn, q in zip(single_wn, single_q) if wn in single_boundsat]
ax.plot(bad_q, bad_wn, 's', ms=10, mfc='none', mec='#1c7293', mew=2,
        label='hankel single-wave, bound-saturated (not genuine)', zorder=7)

# hankel two-wave (q), flagged degenerate ones open
two_ok_wn = [wn for wn in two_wn if wn not in two_flagged]
two_ok_q = [q for wn, q in zip(two_wn, two_q) if wn not in two_flagged]
ax.plot(two_ok_q, two_ok_wn, 's', ms=8, color='#d6604d', mec='white', mew=1,
        label='hankel two-wave ("q" term)', zorder=6)
two_bad_wn = [wn for wn in two_wn if wn in two_flagged]
two_bad_q = [q for wn, q in zip(two_wn, two_q) if wn in two_flagged]
if two_bad_wn:
    ax.plot(two_bad_q, two_bad_wn, 's', ms=10, mfc='none', mec='#d6604d', mew=2,
            label='hankel two-wave, degenerate (Q_2q→inf)', zorder=7)

# fft 2q/2
ax.errorbar(fft_2q2, fft_wn, xerr=[f/2 for f in fft_2q2_fwhm],
            marker='v', color='#762a83', ms=8, mec='white', mew=1, ls='none', capsize=3, elinewidth=1.3,
            label='fft "2q/2" (last peak, halved)', zorder=8)
# fft q direct
ax.errorbar(fft_direct, fft_direct_wn, xerr=[f/2 for f in fft_direct_fwhm],
            marker='^', mfc='none', mec='#e08214', mew=1.8, ms=9, ls='none', capsize=3, elinewidth=1.3,
            label='fft "q" (first peak, direct)', zorder=8)

ax.set_xlim(0, 4.5e5)
ax.set_title(f'graphene_4x1_manual lp2 (interior), MoO3(a), E_F={EF_REF}eV (ref), gamma={GAMMA}cm-1\n'
             f'wn=860-1000cm-1: hankel + FFT, q vs 2q/2 shown side by side (comparison stage)', fontsize=12, fontweight='bold')
ax.legend(frameon=False, fontsize=8.5, loc='upper left')
fig.tight_layout()
out_path = f'{SAVE_DIR}/graphene_4x1_manual_lp2_hankel_fft_qvs2q_compare_Ef{EF_REF}eV.png'
fig.savefig(out_path, dpi=200, bbox_inches='tight')
plt.close(fig)
print(f'Saved {out_path}')
