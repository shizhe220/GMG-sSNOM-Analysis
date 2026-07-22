"""
FFT-only rp cross-check, full wn range (860-1000cm-1), graphene_4x1_manual lp1.

Extends the 860-911cm-1 peak-spacing crosscheck (peak_spacing_crosscheck_graphene_4x1_lp1.py)
now that the user has tuned Lorentzian FFT params (q_guess/search_window/fit_window) for
every wn up through 1000cm-1, not just 900/911. Per the established physics-based labeling
(2026-07-16 log entry): 860-890cm-1 still only resolve a single peak = "2q" (fft_q = peaks[-1]/2);
900cm-1 and above now resolve TWO peaks, [q, 2q] order -- the established "2q, /2-corrected"
source is still peaks[-1]/2 (same purple, consistent with 860-911's already-committed result,
UNCHANGED here), but this script ALSO extracts peaks[0] directly (no /2 -- per the user's
current hypothesis that with these newly-tuned fits, the first peak now cleanly resolves "q"
itself) as a second, separately-colored source, for wn>=900 only (860-890 have no first peak
to extract). This directly revisits the "open question" flagged in the 2026-07-16 log entry
(peaks[0] was unreliable at 900/911 with the OLD, un-retuned params) with the new tuning.

FFT-only this round (per user request) -- no CHT/Hankel/sqrtx/peak-spacing here. Uses the
'amp' channel (established default, FFT_CHANNEL_NB elsewhere in this project). E_F fixed at
0.69eV (not scanned) for a first look.
"""
import sys, os
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append('.')
sys.path.append('/Users/shizhe/envsetting')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import nanoftir_Shizhe as nanof
from snippet.nearfield import rp_schematic
from NearFieldOptics import Materials as M

plt.rcParams.update({
    'font.size': 12, 'font.family': 'Arial',
    'axes.linewidth': 1.5, 'xtick.major.width': 1.5, 'ytick.major.width': 1.5,
    'xtick.direction': 'in', 'ytick.direction': 'in', 'xtick.top': True,
    'ytick.right': True
})

SAVE_DIR = 'figures/rp_dispersion/graphene_4x1_manual_lp1_peak_spacing_crosscheck_v1/fft_all_wn'
os.makedirs(SAVE_DIR, exist_ok=True)

# FFT xr/q_guess/search_window/fit_window, read live from fitting_pipeline_graphene_4x1_lp1.ipynb
# on 2026-07-18 (user just finished retuning 920-1000cm-1 to match the established 900/911
# 2-peak [q, 2q] convention).
FFT_CHANNEL = 'amp'
FFT_PARAMS = {
    860: dict(xr=(0.31, 1.2), q_guess=[1.8], search_window=0.2, fit_window=0.5),
    870: dict(xr=(0.36, 1.5), q_guess=[2], search_window=0.3, fit_window=0.3),
    880: dict(xr=(0.34, 1.5), q_guess=[2], search_window=0.3, fit_window=0.3),
    890: dict(xr=(0.32, 1.5), q_guess=[2.5], search_window=0.3, fit_window=0.3),
    900: dict(xr=(0.33, 1.5), q_guess=[1.2, 2.2], search_window=0.3, fit_window=0.85),
    911: dict(xr=(0.32, 1.5), q_guess=[1.6, 3.0], search_window=0.3, fit_window=0.5),
    920: dict(xr=(0.18, 1.2), q_guess=[1.5, 3], search_window=0.3, fit_window=0.3),
    930: dict(xr=(0.3, 1.2), q_guess=[1.7, 3.2], search_window=0.3, fit_window=0.3),
    941: dict(xr=(0.1, 1.2), q_guess=[2, 3.8], search_window=0.3, fit_window=0.3),
    950: dict(xr=(0.12, 1.2), q_guess=[2, 4], search_window=0.3, fit_window=0.3),
    960: dict(xr=(0.1, 1.2), q_guess=[2, 4.00], search_window=0.3, fit_window=0.3),
    970: dict(xr=(0.14, 1.2), q_guess=[2.1, 4.2], search_window=0.3, fit_window=0.2),
    980: dict(xr=(0.1, 1.2), q_guess=[2.3, 4.6], search_window=0.3, fit_window=0.3),
    991: dict(xr=(0.18, 1.2), q_guess=[2.5, 5], search_window=0.3, fit_window=0.2),
    1000: dict(xr=(0.08, 1.2), q_guess=[4, 8.00], search_window=0.3, fit_window=0.3),
}

rows = []
for wn, fp in FFT_PARAMS.items():
    wn_str = f'{wn}cm-1'
    loaded = nanof.load_aligned_wn_signal(wn_str, 0, data_dir='data/graphene_4x1_manual', L_cutoff_fit=1.5, lp='lp1')
    df_target = loaded['df_target']
    amp_raw, phase_raw = df_target['O3A'], df_target['O3P']
    from scipy.signal import savgol_filter
    wlen = min(41, len(amp_raw) if len(amp_raw) % 2 != 0 else len(amp_raw) - 1)
    amp_osc = amp_raw - savgol_filter(amp_raw, window_length=wlen, polyorder=2)

    fft_out = nanof.plot_channel_fft(
        df_target['distance_um'], amp_osc, phase_raw, label='O3', wn=wn_str,
        xr=fp['xr'], q_range=(0, 10), window='hann', pad_factor=3.0, peak_method='lorentzian',
        q_guess=fp['q_guess'], search_window=fp['search_window'], fit_window=fp['fit_window'], plot=False)
    pk = fft_out[f'peaks_{FFT_CHANNEL}']['peaks']
    fw = fft_out[f'peaks_{FFT_CHANNEL}']['fwhm']

    # Established source, unchanged convention (2026-07-16): last/largest peak = "2q", /2-corrected.
    q_2q = pk[-1] / 10.0 / 2 if pk[-1] is not None else np.nan
    q_2q_fwhm = fw[-1] / 10.0 / 2 if fw[-1] is not None else np.nan

    # New: first peak, direct (no /2) -- only meaningful where 2 peaks were actually fit.
    if len(pk) >= 2 and pk[0] is not None:
        q_first = pk[0] / 10.0
        q_first_fwhm = fw[0] / 10.0
    else:
        q_first = q_first_fwhm = np.nan

    rows.append(dict(wn=wn, fft_q_2q=q_2q, fft_q_2q_fwhm=q_2q_fwhm,
                      fft_q_first=q_first, fft_q_first_fwhm=q_first_fwhm, n_peaks_fit=len(pk)))
    print(f'{wn_str}: 2q/2={q_2q:.3f}+-{q_2q_fwhm/2:.3f}  '
          + (f'first(direct)={q_first:.3f}+-{q_first_fwhm/2:.3f}' if not np.isnan(q_first) else 'first(direct)=n/a'))

df = pd.DataFrame(rows).sort_values('wn').reset_index(drop=True)
df.to_csv(f'{SAVE_DIR}/fft_full_wn_range.csv', index=False)
print(f'\nSaved {SAVE_DIR}/fft_full_wn_range.csv')

# --- rp overlay, FFT only, E_F=0.69eV fixed ---
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

fig, ax = plt.subplots(figsize=(8, 7))
rp_schematic(stack_4x1, materials_4x1_a, subs=SUBS, freqs=FREQS, qs=QS,
             Ef=EF_REF, gamma=GAMMA, figAx=(fig, ax), vmin=0, vmax=8, show_Ef=True, ytickspacing=20)

# established source (unchanged convention, same purple as 860-911's committed figures)
pts = df[['wn', 'fft_q_2q', 'fft_q_2q_fwhm']].dropna(subset=['wn', 'fft_q_2q'])
xerr = pts['fft_q_2q_fwhm'].values * 1e5 / 2
ax.errorbar(pts['fft_q_2q'].values * 1e5, pts['wn'].values, xerr=xerr, marker='v',
            color='#762a83', ms=8, mec='white', mew=1, ls='none', capsize=3, elinewidth=1.3,
            label='fft (amp, 2q/2)', zorder=6)

# new source: first peak, direct -- only wn>=900 have it
pts2 = df[['wn', 'fft_q_first', 'fft_q_first_fwhm']].dropna(subset=['wn', 'fft_q_first'])
xerr2 = pts2['fft_q_first_fwhm'].values * 1e5 / 2
ax.errorbar(pts2['fft_q_first'].values * 1e5, pts2['wn'].values, xerr=xerr2, marker='o',
            color='#d6604d', ms=8, mec='white', mew=1, ls='none', capsize=3, elinewidth=1.3,
            label='fft (amp, 1st peak, direct)', zorder=6)

ax.set_title(f'graphene_4x1_manual lp1 (edge), MoO3(a), E_F={EF_REF}eV (fixed ref), gamma={GAMMA}cm-1\n'
             f'FFT only, wn=860-1000cm-1', fontsize=11)
ax.legend(frameon=False, fontsize=9.5, loc='upper left')
fig.tight_layout()
out_path = f'{SAVE_DIR}/rp_fft_only_860to1000_Ef{EF_REF}eV.png'
fig.savefig(out_path, dpi=200, bbox_inches='tight')
plt.close(fig)
print(f'Saved {out_path}')
