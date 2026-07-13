"""
Peak-spacing cross-check for graphene_4x1_manual lp1, wn=860-911cm-1 -- the wn the user
re-tuned real-space xr_range/lam0_guess_um for in fitting_pipeline_graphene_4x1_lp1.ipynb
(CHT params x_start_cht/L_cutoff/k_fit_range_cm/k_linked_guess_cm confirmed UNCHANGED for
these wn via git diff, so CHT results below match the committed pkl exactly).

Motivation: user noticed a resolvable "q" (not just "2q") peak starting to appear in the
FFT spectrum around 911cm-1 (see rp_dispersion_plot.ipynb section 7 -- below 900cm-1 the
FFT spectrum is dominated by the 2q tip-launched channel). Wants a real-space cross-check:
find the oscillation peaks directly (skip the first one, which sits closest to the tip/edge
and may be contaminated by near-field/edge artifacts not representative of the clean
periodicity), then compare two interval definitions:
  - spacing(peak2, peak3) alone
  - mean(spacing(peak1, peak2), spacing(peak2, peak3))
against the already-established CHT/Hankel/1sqrtx q_p, to see which one the real-space
periodicity actually tracks.

VERSION_SUFFIX below -- bump before a new comparison round, do not overwrite past ones.
Prototype/investigation only -- does NOT touch the committed pkl or fit_comparison figures.
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
from scipy.signal import find_peaks
from scipy.optimize import minimize_scalar, curve_fit
import nanoftir_Shizhe as nanof
from snippet.nearfield import rp_schematic
from NearFieldOptics import Materials as M

plt.rcParams.update({
    'font.size': 12, 'font.family': 'Arial',
    'axes.linewidth': 1.5, 'xtick.major.width': 1.5, 'ytick.major.width': 1.5,
    'xtick.direction': 'in', 'ytick.direction': 'in', 'xtick.top': True,
    'ytick.right': True
})

VERSION_SUFFIX = 'v1'
DATASET = 'graphene_4x1_manual'
LP = 'lp1'  # "edge" linecut
DATA_DIR = f'data/{DATASET}'
SAVE_DIR = f'figures/rp_dispersion/{DATASET}_{LP}_peak_spacing_crosscheck_{VERSION_SUFFIX}'
os.makedirs(SAVE_DIR, exist_ok=True)

# Current (just-retuned) params, copied verbatim from fitting_pipeline_graphene_4x1_lp1.ipynb
# on 2026-07-13 -- CHT columns confirmed unchanged from the committed v1 pkl for these wn.
WN_PARAMS = {
    860: dict(x_start_cht=0.30, L_cutoff_cht=1.5, k_fit_range_cm=(0.1, 4.0), k_linked_guess_cm=1.0,
              xr_range_rs=(0.31, 1.5), lam0_guess_um=0.67),
    870: dict(x_start_cht=0.13, L_cutoff_cht=1.2, k_fit_range_cm=(0.2, 6.0), k_linked_guess_cm=1.0,
              xr_range_rs=(0.36, 1.5), lam0_guess_um=0.58),
    880: dict(x_start_cht=0.13, L_cutoff_cht=1.2, k_fit_range_cm=(0.2, 6.0), k_linked_guess_cm=1.6,
              xr_range_rs=(0.34, 1.5), lam0_guess_um=0.60),
    890: dict(x_start_cht=0.12, L_cutoff_cht=1.2, k_fit_range_cm=(0.2, 6.0), k_linked_guess_cm=1.8,
              xr_range_rs=(0.32, 1.5), lam0_guess_um=0.50),
    900: dict(x_start_cht=0.11, L_cutoff_cht=1.2, k_fit_range_cm=(0.2, 7.0), k_linked_guess_cm=1.8,
              xr_range_rs=(0.34, 1.5), lam0_guess_um=0.50),
    911: dict(x_start_cht=0.18, L_cutoff_cht=1.2, k_fit_range_cm=(0.2, 8.0), k_linked_guess_cm=1.8,
              xr_range_rs=(0.32, 1.5), lam0_guess_um=0.40),
}

rows = []
for wn, p in WN_PARAMS.items():
    wn_str = f'{wn}cm-1'
    loaded = nanof.load_aligned_wn_signal(wn_str, 0, data_dir=DATA_DIR, L_cutoff_fit=1.5, lp='lp1')
    x_f, sig_f = loaded['x_f'], loaded['sig_f']

    cht_results, fig_cht = nanof.fit_and_plot_cht(
        x_f, sig_f, wn_str, x_start_cht=p['x_start_cht'], L_cutoff_cht=p['L_cutoff_cht'],
        k_fit_range_cm=p['k_fit_range_cm'], k_linked_guess_cm=p['k_linked_guess_cm'])
    plt.close(fig_cht)
    cht_q = cht_results['q_p_1e5cm_1']

    amplp = pd.DataFrame({'distance_um': x_f, f'{wn_str}_O3A': sig_f})
    outs, fig_rs, _ = nanof.compare_cavity_models(
        amplp, f'{wn_str}_O3A', xr=p['xr_range_rs'], yc_um=1.9, fit_yc=False, edges='single',
        prefactors=('hankel', '1/sqrtx'), win=3, prom=0.01, lam0_guess=p['lam0_guess_um'],
        ylim=(sig_f.min() * 1.5, sig_f.max() * 0.5), figsize=(8, 6))
    plt.close(fig_rs)
    hankel_q = outs['hankel']['derived']['q_cm^-1'] / 1e5
    sqrtx_q = outs['1/sqrtx']['derived']['q_cm^-1'] / 1e5

    # Real-space peak spacing, skip the first peak (closest to tip/edge). Find peaks on
    # the FULL signal (not restricted to xr_range_rs) -- xr_range_rs's lower bound already
    # excludes the true first peak, so restricting to it before peak-finding made "peak1"
    # actually the 2nd real peak and "peak3" the 4th, which is deep in the damping tail
    # (weak, noisy). Using the full x_f/sig_f keeps peak1 as the genuine tip/edge peak.
    # Prominence threshold is scaled to the ptp of the DECAYING TAIL only (x>0.15um, past
    # the dominant first peak) -- using the full signal's ptp for the threshold would let
    # the huge first peak drown out the smaller peak2/3 (that's what happened on the first
    # attempt: threshold scaled to full-signal ptp filtered out real peak3 for 860/900/911).
    xw, sw = x_f, sig_f
    tail_mask = xw > 0.25  # past the first (dominant) peak for all 6 wn here (peak1 <=210nm)
    prom_ref = np.ptp(sw[tail_mask]) if tail_mask.any() else np.ptp(sw)
    peaks_idx, _ = find_peaks(sw, prominence=0.08 * prom_ref)

    def lorentzian(x, B, A, x0, gamma):
        return B + A / (1 + ((x - x0) / gamma) ** 2)

    def subpixel_peak_x(idx, half_win=3):
        """Fit a Lorentzian to a +-half_win window of points around a find_peaks index and
        return the fitted center -- removes the pixel-grid quantization (measured spacing
        was landing on exact integer multiples of the 15.38nm pixel pitch, e.g. 870/880cm-1
        both snapping to 18px) and is more robust to noise near the peak than just taking
        the raw grid maximum. Previously used a 3-point parabola; switched to Lorentzian per
        user request. First attempt (half_win=6, unbounded) let gamma blow up to 1000s of nm
        (a nearly-flat degenerate fit) for several wn -- visually confirmed via zoomed
        preview (user noticed raw data visibly poking out past the fit curve). Fixed by (a)
        narrowing the window to +-3 points (this peak shape isn't a sharp/narrow resonance,
        so a wide window just invites the fit to go broad/degenerate) and (b) bounding gamma
        to a physically sane range relative to the window span. Re-verified visually
        (scratchpad/peak2_zoom_narrowwin.png) before adopting -- gamma now 80-140nm, curves
        hug the actual peak shape. Falls back to the raw grid position if the fit fails."""
        lo, hi = max(0, idx - half_win), min(len(xw), idx + half_win + 1)
        xloc, yloc = xw[lo:hi], sw[lo:hi]
        B0 = np.min(yloc)
        A0 = yloc[idx - lo] - B0
        win_span = xloc[-1] - xloc[0]
        gamma0 = win_span / 4
        bounds = ([-np.inf, 0, xloc[0], (xloc[1] - xloc[0]) * 0.3], [np.inf, np.inf, xloc[-1], win_span * 1.5])
        try:
            popt, _ = curve_fit(lorentzian, xloc, yloc, p0=[B0, A0, xw[idx], gamma0], bounds=bounds, maxfev=5000)
            return popt[2]
        except Exception:
            return xw[idx]

    # Only peak2-peak3 (most reliable pair -- peak1 is skipped as tip/edge, peak4+ are
    # deep in the noisy damping tail per the 860/911cm-1 visual check), using sub-pixel
    # refined positions instead of the raw grid-point location.
    if len(peaks_idx) >= 3:
        x2_fit = subpixel_peak_x(peaks_idx[1])
        x3_fit = subpixel_peak_x(peaks_idx[2])
        interval_23 = x3_fit - x2_fit
        q_23 = (2 * np.pi / (interval_23 * 1e-4)) / 1e5       # spacing(peak2,peak3), fit position
    else:
        interval_23 = q_23 = np.nan
    q_avg12_23 = np.nan  # dropped per user request -- peak2-peak3 alone is the reliable pair

    rows.append(dict(wn=wn, cht_q=cht_q, hankel_q=hankel_q, sqrtx_q=sqrtx_q,
                      n_peaks=len(peaks_idx), interval_23_um_fit=interval_23, q_peak23_fit=q_23))
    print(f'{wn}cm-1: CHT={cht_q:.3f}  Hankel={hankel_q:.3f}  sqrtx={sqrtx_q:.3f}  '
          f'n_peaks={len(peaks_idx)}  q(peak2-3, subpixel fit)={q_23:.3f}')

df = pd.DataFrame(rows).sort_values('wn', ascending=False)

# Real-space peak-spacing interval was confirmed (previous run) to track the "2q"
# tip-launched channel, not "q" -- divide by 2 here so it's on the same q scale as
# CHT/Hankel/sqrtx for direct comparison.
df['q_peak23_fit_div2'] = df['q_peak23_fit'] / 2

df.to_csv(f'{SAVE_DIR}/peak_spacing_crosscheck.csv', index=False)
print(f'\nSaved {SAVE_DIR}/peak_spacing_crosscheck.csv')

# ---- Combined rp plot: 4 q sources x 6 wn on ONE Im(r_p) map, one panel per E_F ----
# CHT dropped from the comparison (redundant with Hankel/sqrtx here); peak-spacing sources
# use the /2-corrected ("q, not 2q") columns. Label reads "peak3-peak2" (later minus
# earlier, matching how the interval is actually computed: xp[2]-xp[1]), not "peak2-3".
GAMMA = 20
phonon_a = [(972, 4, 820, 4)]
phonon_b = [(1004, 2, 958, 2)]
MoO3_a = M.AnisotropicMaterial(eps_infinity=[4, 5.2, 2.4], phonon_params=[phonon_a, phonon_a, phonon_b])
materials_4x1 = {'SiO2': M.SiO2_300nm, 'MoO3': MoO3_a}
SUBS = M.Si
stack_4x1 = [('MoO3', 4.2), 'graphene', ('SiO2', 300)]

FREQS = np.linspace(850, 920, 300)
QS = np.linspace(1e4, 3.5e5, 400)

sources = [('hankel', 'hankel_q', '#1c7293', 's'), ('sqrtx', 'sqrtx_q', '#e08214', '^'),
           ('peak3-peak2 / 2 (subpixel fit)', 'q_peak23_fit_div2', '#4daf4a', 'D')]

# ---- Which doping level is actually correct: fit E_F per source (6 points each) ----
# Same "maximize sum of Im(r_p) at the (w,q) data points" recipe as rp_dispersion_plot.ipynb
# -- this is the direct answer to "which E_F makes these points fall ON the ridge peak",
# not eyeballing a handful of fixed-E_F panels.
def build_tm_layers(Ef, gamma):
    Efcm = Ef * 8065.5
    graphene = M.SingleLayerGraphene(chemical_potential=Efcm, gamma=gamma)
    return [(MoO3_a, 4.2e-7), graphene, (M.SiO2_300nm, 300e-7)]

def fit_Ef(w_data, q_data_1e5, gamma=GAMMA, ef_bounds=(0.1, 1.5)):
    w_data = np.asarray(w_data, dtype=float)
    q_data = np.asarray(q_data_1e5, dtype=float) * 1e5
    def neg_score(Ef):
        layersTM = M.LayeredMediaTM(*build_tm_layers(Ef, gamma), exit=SUBS)
        rp_grid = np.asarray(layersTM.reflection_p(w_data, q=q_data))
        return -float(np.sum(np.imag(np.diag(rp_grid))))
    res = minimize_scalar(neg_score, bounds=ef_bounds, method='bounded')
    return float(res.x)

fitted_ef = {}
print('\nFitted E_F per source (6 pts, wn=860-911cm-1):')
for label, col, color, marker in sources:
    pts_df = df[['wn', col]].dropna()
    ef_fit = fit_Ef(pts_df['wn'].values, pts_df[col].values)
    fitted_ef[label] = ef_fit
    print(f'  {label:>20s}: E_F = {ef_fit:.3f} eV')
print(f'  Raman-measured (edge): ~0.60-0.65 eV')

with open(f'{SAVE_DIR}/fitted_ef_per_source.csv', 'w') as f:
    f.write('source,fitted_Ef_eV\n')
    for label, ef_fit in fitted_ef.items():
        f.write(f'{label},{ef_fit:.4f}\n')

# 0.65eV = Raman-measured value for this dataset (graphene_4x1 lp1 = "edge");
# 0.60eV = lower end of the Raman-quoted range; 0.67eV = user-requested extra check point;
# 0.771eV = the previously-used full-15-wn CHT/MoO3(a) fit value -- one panel each, same
# 4 sources, for direct visual comparison. Fitted per-source E_F values (above) are the
# actual answer to "which doping level is correct"; these fixed-E_F panels are a
# complementary visual (same 4 points, ridge redrawn at each candidate E_F).
for ef_label, EF_REF in [('0.65eV_raman', 0.65), ('0.60eV_raman_lower', 0.60),
                          ('0.67eV', 0.67), ('0.771eV_cht_fit', 0.771)]:
    fig, ax = plt.subplots(figsize=(7, 6))
    rp_schematic(stack_4x1, materials_4x1, subs=SUBS, freqs=FREQS, qs=QS,
                 Ef=EF_REF, gamma=GAMMA, figAx=(fig, ax), vmin=0, vmax=8, show_Ef=True, ytickspacing=10)
    for label, col, color, marker in sources:
        pts_df = df[['wn', col]].dropna()
        if len(pts_df):
            # QS/the rp_schematic axes are in raw cm^-1 -- col values are in 1e5 cm^-1, convert
            ax.plot(pts_df[col].values * 1e5, pts_df['wn'].values, marker, color=color, ms=8,
                    mec='white', mew=1, ls='none', label=label, zorder=6)
    ax.set_title(f'{DATASET} {LP} (edge), MoO3(a), E_F={EF_REF}eV (fixed ref), gamma={GAMMA}cm-1\n'
                 f'wn=860-911cm-1 peak-spacing cross-check (peak-spacing sources /2-corrected)', fontsize=11)
    ax.legend(frameon=False, fontsize=9, loc='upper left')
    fig.tight_layout()
    out_path = f'{SAVE_DIR}/rp_peak_spacing_crosscheck_860to911_Ef{ef_label}.png'
    fig.savefig(out_path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f'Saved {out_path}')

# ---- Real-space Lorentzian-fit preview, all 6 wn, peak2/peak3 marked + 2q/q labeled ----
# Same fit recipe as subpixel_peak_x above (Lorentzian, +-3pt window, bounded gamma) --
# redone here (rather than reusing the per-wn closures above) purely for plotting.
def lorentzian_preview(x, B, A, x0, gamma):
    return B + A / (1 + ((x - x0) / gamma) ** 2)

def fit_peak_preview(xw, sw, idx, half_win=3):
    lo, hi = max(0, idx - half_win), min(len(xw), idx + half_win + 1)
    xloc, yloc = xw[lo:hi], sw[lo:hi]
    B0 = np.min(yloc)
    A0 = yloc[idx - lo] - B0
    win_span = xloc[-1] - xloc[0]
    gamma0 = win_span / 4
    bounds = ([-np.inf, 0, xloc[0], (xloc[1] - xloc[0]) * 0.3], [np.inf, np.inf, xloc[-1], win_span * 1.5])
    popt, _ = curve_fit(lorentzian_preview, xloc, yloc, p0=[B0, A0, xw[idx], gamma0], bounds=bounds, maxfev=5000)
    return popt, xloc, yloc

fig, axes = plt.subplots(2, 3, figsize=(16, 9))
for ax, wn in zip(axes.flat, WN_PARAMS.keys()):
    loaded = nanof.load_aligned_wn_signal(f'{wn}cm-1', 0, data_dir=DATA_DIR, L_cutoff_fit=1.5, lp='lp1')
    x_f, sig_f = loaded['x_f'], loaded['sig_f']
    tail_mask = x_f > 0.25
    prom_ref = np.ptp(sig_f[tail_mask])
    peaks_idx, _ = find_peaks(sig_f, prominence=0.08 * prom_ref)

    ax.plot(x_f * 1000, sig_f, '-', color='#555555', lw=1.2, alpha=0.7, zorder=1, label='data')
    ax.axhline(0, color='k', ls='--', lw=0.6, alpha=0.4)

    fit_x0 = {}
    colors = {'peak2': '#b2182b', 'peak3': '#1c7293'}
    for pk_i, label in zip([1, 2], ['peak2', 'peak3']):
        idx = peaks_idx[pk_i]
        popt, xloc, yloc = fit_peak_preview(x_f, sig_f, idx)
        fit_x0[label] = popt[2]
        xdense = np.linspace(xloc[0], xloc[-1], 200)
        ax.plot(xdense * 1000, lorentzian_preview(xdense, *popt), '-', color=colors[label], lw=2.2, zorder=4,
                label=f'{label} fit, x0={popt[2] * 1000:.1f}nm')
        ax.plot(popt[2] * 1000, lorentzian_preview(popt[2], *popt), '*', color=colors[label], ms=18,
                mec='black', mew=0.8, zorder=5)

    interval_23_preview = fit_x0['peak3'] - fit_x0['peak2']
    q_2q_preview = (2 * np.pi / (interval_23_preview * 1e-4)) / 1e5  # 1e5 cm^-1, pre /2-correction
    ax.text(0.97, 0.95, f'2q = {q_2q_preview:.2f}' + r'$\times10^5$cm$^{-1}$' +
            f'\nq = {q_2q_preview / 2:.2f}' + r'$\times10^5$cm$^{-1}$',
            transform=ax.transAxes, ha='right', va='top', fontsize=10.5, fontweight='bold',
            color='#4daf4a', bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='#4daf4a', alpha=0.9))

    ax.set_title(f'{wn}cm-1', fontsize=13, fontweight='bold')
    ax.set_xlabel('Distance from edge (nm)', fontweight='bold')
    ax.set_ylabel('Re signal (bg-sub, a.u.)', fontweight='bold')
    ax.set_xlim(0, 1200)
    ax.legend(fontsize=8.5, frameon=False, loc='upper left')
    ax.tick_params(direction='in', top=True, right=True)

fig.tight_layout()
out_path = f'{SAVE_DIR}/lorentzian_fit_realspace_preview_860to911.png'
fig.savefig(out_path, dpi=150, bbox_inches='tight')
plt.close(fig)
print(f'Saved {out_path}')
