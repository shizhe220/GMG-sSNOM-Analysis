"""
rp(q,w) overlay for graphene_4x1_manual lp2 (interior), hankel only -- single-wave
(860-890cm-1) + two-wave (900-1000cm-1), same convention as lp1's
plot_rp_hankel_fft_peak23_crosscheck.py. FFT and real-space peak2-peak3 spacing are NOT
included yet (lp2 hasn't had those done -- per user request, 2026-07-23).

Values read live from fitting_pipeline_graphene_4x1_lp2.ipynb's own executed outputs
(not re-fit here) -- last synced 2026-07-23 after the user re-tuned 920/960's two-wave
cells (previously degenerate, now genuine -- see log entry history for that fix).

IMPORTANT: many of lp2's single-wave Hankel fits are NOT genuine local minima -- they
silently saturate against fit_cavity_prefactor_compare's default lam_bound_factor=(0.7,1.2)
bound around lam0_guess. Checked systematically (see log entry): of the 15 wn, 7
(911,930,941,950,960,980,991) hit this bound exactly. For 900-1000, the two-wave model
(intended fix) is used instead of the bound-prone single-wave one, matching lp1's
convention. All flagged points are marked with a hatched/open marker below rather than
silently plotted as trustworthy.
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

# --- rp background (same MoO3-a / stack as lp1, same sample) ---
phonon_a = [(972, 4, 820, 4)]
phonon_b = [(1004, 2, 958, 2)]
MoO3_a = M.AnisotropicMaterial(eps_infinity=[4, 5.2, 2.4], phonon_params=[phonon_a, phonon_a, phonon_b])
materials_4x1_a = {'SiO2': M.SiO2_300nm, 'MoO3': MoO3_a}
SUBS = M.Si
stack_4x1 = [('MoO3', 4.2), 'graphene', ('SiO2', 300)]
GAMMA = 20
EF_REF = 0.69  # same reference value used for lp1, for a like-for-like comparison
FREQS = np.linspace(850, 1010, 300)
QS = np.linspace(1e4, 9e5, 400)

# single-wave Hankel, 860-890cm-1 -- all genuine (non-bound-saturated) convergence
single_wn = [860, 870, 880, 890]
single_q = [94900, 104000, 109000, 111000]

# two-wave Hankel "q" term, 900-1000cm-1 -- re-tuned 2026-07-23 (920/960 no longer
# degenerate after user re-tuned xr_range_2wave/q_guess_2wave for those two)
two_wn =     [900,    911,    920,    930,    941,    950,    960,    970,    980,    991,    1000]
two_q =      [150000, 155000, 150000, 202000, 191000, 191000, 206000, 265000, 189000, 206000, 257000]
two_flagged = set()

fig, ax = plt.subplots(figsize=(9, 8))
rp_schematic(stack_4x1, materials_4x1_a, subs=SUBS, freqs=FREQS, qs=QS,
             Ef=EF_REF, gamma=GAMMA, figAx=(fig, ax), vmin=0, vmax=8, show_Ef=True, ytickspacing=20)

ax.plot(single_q, single_wn, 's', ms=9, color='#1c7293', mec='white', mew=1,
        label='hankel (single-wave, 860-890)', zorder=6)

two_ok_wn = [wn for wn in two_wn if wn not in two_flagged]
two_ok_q = [q for wn, q in zip(two_wn, two_q) if wn not in two_flagged]
ax.plot(two_ok_q, two_ok_wn, 's', ms=9, color='#d6604d', mec='white', mew=1,
        label='hankel (two-wave, 900-1000)', zorder=6)

two_bad_wn = [wn for wn in two_wn if wn in two_flagged]
two_bad_q = [q for wn, q in zip(two_wn, two_q) if wn in two_flagged]
if two_bad_wn:
    ax.plot(two_bad_q, two_bad_wn, 's', ms=11, mfc='none', mec='#d6604d', mew=2.2,
            label='hankel (two-wave, degenerate -- Q_2q→inf, low confidence)', zorder=7)

ax.set_xlim(0, 4.5e5)
ax.set_title(f'graphene_4x1_manual lp2 (interior), MoO3(a), E_F={EF_REF}eV (ref, shared w/ lp1), gamma={GAMMA}cm-1\n'
             f'wn=860-1000cm-1: hankel only (no FFT/peak2-3 yet)', fontsize=12, fontweight='bold')
ax.legend(frameon=False, fontsize=9, loc='upper left')
fig.tight_layout()
out_path = f'{SAVE_DIR}/graphene_4x1_manual_lp2_hankel_only_Ef{EF_REF}eV.png'
fig.savefig(out_path, dpi=200, bbox_inches='tight')
plt.close(fig)
print(f'Saved {out_path}')
