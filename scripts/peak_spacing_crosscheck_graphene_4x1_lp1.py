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
from scipy.signal import find_peaks, savgol_filter
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

# FFT xr/q_guess/search_window/fit_window, re-read live from the notebook on 2026-07-16
# (user finished retuning all 6 wn with peak_method='lorentzian' + the new fit_window
# param, see envsetting/snippet/dispersion_fitting.py 2026-07-16 log entry).
#
# Physics-based labeling (per user correction 2026-07-16 -- label by which PEAK IS
# PHYSICALLY "2q", not by its index in q_guess): 860/870/880/890 have a single q_guess
# (targeting "2q" directly, since at these wn the FFT is 2q-dominated and there's no
# separately-resolvable "q" peak) -- fft_q = peaks[0]/2. 900/911 have TWO q_guess entries
# ([q, 2q] order), and per the 2026-07-16 visual/numeric check the SECOND (larger) one is
# still the reliable "2q" signal -- the first (smaller) one's fit was essentially
# unconstrained at 900cm-1 (FWHM=1.2) and undershot the ridge at 911cm-1, so it's dropped
# here rather than plotted as if it were a validated "q" measurement. fft_q = peaks[-1]/2
# for ALL wn now (last peak = 2q, works whether 1 or 2 peaks were fit).
#
# plot_channel_fft always Lorentzian-fits all 3 spectra (amp/phase/complex) internally in
# one call -- previously only the FFT_CHANNEL-selected one was pulled out and the other two
# were silently discarded (computed but not saved), so switching channels meant re-running
# this whole script (the slow part is CHT/Hankel/sqrtx/peak-spacing, not FFT itself, but
# there was no way to compare channels without a re-run). Now all 3 are saved to the CSV
# (fft_q_amp/fft_q_phase/fft_q_complex + _fwhm) so the notebook can switch interactively.
FFT_CHANNELS = ['amp', 'phase', 'complex']
FFT_PARAMS = {
    860: dict(xr=(0.31, 1.2), q_guess=[1.8], search_window=0.2, fit_window=0.5),
    870: dict(xr=(0.36, 1.5), q_guess=[2],   search_window=0.3, fit_window=0.3),
    880: dict(xr=(0.34, 1.5), q_guess=[2],   search_window=0.3, fit_window=0.3),
    890: dict(xr=(0.32, 1.5), q_guess=[2.5], search_window=0.3, fit_window=0.3),
    900: dict(xr=(0.33, 1.5), q_guess=[1.2, 2.2], search_window=0.3, fit_window=0.85),
    911: dict(xr=(0.32, 1.5), q_guess=[1.6, 3.0], search_window=0.3, fit_window=0.5),
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

    # FFT: Lorentzian peak fit (peak_method='lorentzian' + fit_window, see
    # envsetting/snippet/dispersion_fitting.py 2026-07-16 fix) on all 3 channels
    # (amp/phase/complex, all computed in this one call), with the user's per-wn
    # q_guess/search_window/fit_window. fft_q_{channel} is always peaks[-1]/2 (the
    # last/largest-q peak = "2q" physically, /2-corrected) -- see FFT_PARAMS note above
    # for why the 900/911 "first peak = q directly" candidate was dropped rather than
    # plotted.
    fp = FFT_PARAMS[wn]
    df_target = loaded['df_target']
    amp_raw, phase_raw = df_target['O3A'], df_target['O3P']
    wlen = min(41, len(amp_raw) if len(amp_raw) % 2 != 0 else len(amp_raw) - 1)
    amp_osc = amp_raw - savgol_filter(amp_raw, window_length=wlen, polyorder=2)
    fft_out = nanof.plot_channel_fft(
        df_target['distance_um'], amp_osc, phase_raw, label='O3', wn=wn_str,
        xr=fp['xr'], q_range=(0, 10), window='hann', pad_factor=3.0, peak_method='lorentzian',
        q_guess=fp['q_guess'], search_window=fp['search_window'], fit_window=fp['fit_window'], plot=False)

    fft_q_by_channel = {}
    fft_q_fwhm_by_channel = {}
    for ch in FFT_CHANNELS:
        pk, fw = fft_out[f'peaks_{ch}']['peaks'], fft_out[f'peaks_{ch}']['fwhm']
        fft_q_by_channel[ch] = pk[-1] / 10.0 / 2 if pk[-1] is not None else np.nan
        fft_q_fwhm_by_channel[ch] = fw[-1] / 10.0 / 2 if fw[-1] is not None else np.nan

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

    row = dict(wn=wn, cht_q=cht_q, hankel_q=hankel_q, sqrtx_q=sqrtx_q,
               n_peaks=len(peaks_idx), interval_23_um_fit=interval_23, q_peak23_fit=q_23)
    for ch in FFT_CHANNELS:
        row[f'fft_q_{ch}'] = fft_q_by_channel[ch]
        row[f'fft_q_{ch}_fwhm'] = fft_q_fwhm_by_channel[ch]
    rows.append(row)
    fft_str = '  '.join(f'FFT({ch})={fft_q_by_channel[ch]:.3f}+-{fft_q_fwhm_by_channel[ch]/2:.3f}' for ch in FFT_CHANNELS)
    print(f'{wn}cm-1: CHT={cht_q:.3f}  Hankel={hankel_q:.3f}  sqrtx={sqrtx_q:.3f}  '
          f'{fft_str}  q(peak2-3, subpixel fit)={q_23:.3f}')

df = pd.DataFrame(rows).sort_values('wn', ascending=False)

# Real-space peak-spacing interval was confirmed (previous run) to track the "2q"
# tip-launched channel, not "q" -- divide by 2 here so it's on the same q scale as
# CHT/Hankel/sqrtx for direct comparison.
df['q_peak23_fit_div2'] = df['q_peak23_fit'] / 2

df.to_csv(f'{SAVE_DIR}/peak_spacing_crosscheck.csv', index=False)
print(f'\nSaved {SAVE_DIR}/peak_spacing_crosscheck.csv')

# ---- rp plotting + E_F scan MOVED to rp_dispersion_plot.ipynb (section 8), 2026-07-16 ----
# per user request, so the source-selection/E_F-range/etc. can be tuned interactively
# without re-running the CHT/Hankel/sqrtx/peak-spacing/FFT fits above (which are the slow
# part). That notebook section loads peak_spacing_crosscheck.csv (saved above) and reuses
# MoO3_a/materials_4x1_a/stack_4x1/SUBS/GAMMA/build_tm_layers/fit_Ef already defined
# earlier in the notebook.

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
