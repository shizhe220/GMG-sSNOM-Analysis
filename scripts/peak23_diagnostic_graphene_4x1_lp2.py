"""
Real-space peak2-peak3 spacing diagnostic for graphene_4x1_manual lp2 (interior),
wn=860-1000cm-1 -- sibling of lp1's own diagnostic (2026-07-13/2026-07-22, never saved
as a permanent script -- only the output PNGs were kept). This one IS saved, so re-runs
and per-wn overrides don't require rebuilding the whole thing from scratch again.

Method (matching lp1 exactly): find_peaks on the full background-subtracted real-space
signal, prominence = PROM_REL * ptp(signal in the tail region x>0.25um) by default;
peak2/peak3 sub-pixel-refined via a bounded Lorentzian fit (B + A/(1+((x-x0)/gamma)^2))
over a +-3-point half-window around the raw find_peaks index. q from spacing: this is
the "2q" convention (peak2-peak3 spacing tracks the tip round-trip "2q" wave, same as
everywhere else in this project), so q = 2*pi/(x3_fit-x2_fit) is labeled "2q" and half
of that is "q".

This is a first pass with uniform default settings across all 15 wn -- expect some wn
(especially high-wn, many spurious/shallow peaks) to need per-wn prominence overrides or
a raw-grid fallback where the Lorentzian fit is degenerate, same as lp1's 941/1000 needed
(see PROM_ABS_OVERRIDE / RAW_GRID_OVERRIDE below, empty until the user reviews and flags
specific wn).
"""
import os
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys
sys.path.append('.')
sys.path.append('/Users/shizhe/envsetting')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
import nanoftir_Shizhe as nanof

plt.rcParams.update({'font.size': 11, 'font.family': 'Arial'})

SAVE_DIR = 'figures/rp_dispersion/graphene_4x1_manual_lp2_hankel_fft_compare_v1/peak23_diagnostic'
os.makedirs(SAVE_DIR, exist_ok=True)

WN_LOW = [860, 870, 880, 890, 900, 911]
WN_HIGH = [920, 930, 941, 950, 960, 970, 980, 991, 1000]

PROM_REL = 0.08          # default: relative to tail-region ptp
PROM_ABS_OVERRIDE = {}   # wn -> absolute prominence override, fill in after review
RAW_GRID_OVERRIDE = set()  # {(wn, 'peak2'|'peak3')}, fill in after review
HALF_WIN = 3

# Manual position hints from user visual review (2026-07-23), nm from edge. When set,
# the nearest actual local peak (from a low-threshold find_peaks pass, prominence=0.02x
# tail ptp) to this hint replaces the automatic peaks_idx[1]/[2] pick for that (wn,label).
# 950/960/970 had no exact position given by the user ("not sure, guess from the
# established two-wave/FFT q") -- those hints were back-calculated from the established
# q. 930/980/991/1000 revised again after checking raw y-values directly (930/1000: the
# user's original hint was on the rising edge of a lobe, not its actual peak; 980/991:
# the user caught that one point was sitting in a valley, not a peak) -- see log entry.
MANUAL_PEAK_XNM = {
    (870, 'peak3'): 780,
    (880, 'peak3'): 800,
    (890, 'peak2'): 424,   # not given directly; back-filled from peak3=700 minus expected spacing
    (890, 'peak3'): 700,
    (930, 'peak2'): 404,   # matches user's ~400 -- genuine local max of that lobe
    (930, 'peak3'): 727,   # revised from user's ~650 (rising edge, not the lobe's actual peak)
    (950, 'peak2'): 408,   # low confidence -- back-calculated from established q=1.91e5
    (950, 'peak3'): 593,
    (960, 'peak2'): 267,   # low confidence -- back-calculated from established q=2.06e5
    (960, 'peak3'): 421,
    (970, 'peak2'): 365,   # low confidence -- back-calculated from established q=2.65e5
    (970, 'peak3'): 504,
    (980, 'peak2'): 471,   # revised: 394 was a minor sub-peak, 471 is the lobe's true max
    (980, 'peak3'): 748,   # revised: 548 was sitting in a valley (user caught this)
    (991, 'peak2'): 375,   # revised: 283 was sitting in a valley (user caught this)
    (991, 'peak3'): 622,   # judgment call -- 452 is closer but implies an implausibly
                           # small spacing; 622 matches the periodicity trend of nearby wn
    (1000, 'peak2'): 364,  # matches user's ~400 well enough -- genuine local max
    (1000, 'peak3'): 534,  # revised from user's ~600 (already past the lobe's actual peak)
}


def lorentzian(x, B, A, x0, gamma):
    return B + A / (1 + ((x - x0) / gamma) ** 2)


def load_signal(wn):
    loaded = nanof.load_aligned_wn_signal(f'{wn}cm-1', 0, data_dir='data/graphene_4x1_manual',
                                           L_cutoff_fit=1.5, lp='lp2')
    return loaded['x_f'], loaded['sig_f']


def find_and_fit_peaks(wn, x_f, sig_f):
    x_nm = x_f * 1000
    tail_mask = x_f > 0.25
    prom_ref = np.ptp(sig_f[tail_mask]) if tail_mask.sum() > 3 else np.ptp(sig_f)
    prominence = PROM_ABS_OVERRIDE.get(wn, PROM_REL * prom_ref)

    peaks_idx, _ = find_peaks(sig_f, prominence=prominence)
    default_idx = dict(zip(['peak2', 'peak3'], peaks_idx[1:3] if len(peaks_idx) >= 3 else []))

    # low-threshold candidate list, used only to resolve MANUAL_PEAK_XNM hints to the
    # nearest actual local peak (rather than trusting the hint's exact nm value, which
    # came from an approximate visual read)
    fine_idx, _ = find_peaks(sig_f, prominence=0.02 * prom_ref)

    label_idx = {}
    for label in ('peak2', 'peak3'):
        if (wn, label) in MANUAL_PEAK_XNM and len(fine_idx) > 0:
            target = MANUAL_PEAK_XNM[(wn, label)]
            label_idx[label] = fine_idx[np.argmin(np.abs(x_nm[fine_idx] - target))]
        elif label in default_idx:
            label_idx[label] = default_idx[label]

    results = {}
    for label, idx in label_idx.items():
        if (wn, label) in RAW_GRID_OVERRIDE:
            results[label] = dict(x0=x_nm[idx], fit=None)
            continue
        lo = max(0, idx - HALF_WIN)
        hi = min(len(x_nm), idx + HALF_WIN + 1)
        xw, yw = x_nm[lo:hi], sig_f[lo:hi]
        try:
            p0 = [0, yw[HALF_WIN if idx - lo == HALF_WIN else len(yw) // 2] - np.median(yw),
                  x_nm[idx], (xw[-1] - xw[0]) / 4]
            popt, _ = curve_fit(lorentzian, xw, yw, p0=p0, maxfev=5000)
            x0_fit = popt[2]
            gamma_fit = abs(popt[3])
            amp_fit = abs(popt[1])
            span = xw[-1] - xw[0]
            local_ptp = np.ptp(yw)
            # degenerate: gamma blown past the fit window, x0 outside it, or amplitude
            # wildly larger than the local data range (runaway fit on near-flat/noisy data)
            if gamma_fit > span or not (xw[0] <= x0_fit <= xw[-1]) or amp_fit > 10 * max(local_ptp, 1e-6):
                results[label] = dict(x0=x_nm[idx], fit=None)  # degenerate, fall back
            else:
                results[label] = dict(x0=x0_fit, fit=(xw, popt))
        except Exception:
            results[label] = dict(x0=x_nm[idx], fit=None)
    return peaks_idx, results, prominence


def panel(ax, wn, x_f, sig_f, show_all_peaks=False):
    x_nm = x_f * 1000
    peaks_idx, results, prominence = find_and_fit_peaks(wn, x_f, sig_f)
    ax.plot(x_nm, sig_f, '-', color='gray', lw=1.3, label='data', zorder=1)
    ax.axhline(0, color='k', ls='--', lw=0.6, alpha=0.4)

    if show_all_peaks:
        ax.plot(x_nm[peaks_idx], sig_f[peaks_idx], 'o', color='lightgray', ms=6,
                mec='gray', zorder=2, label=f'all {len(peaks_idx)} found peaks')

    colors = {'peak2': '#b2182b', 'peak3': '#1c7293'}
    for label in ('peak2', 'peak3'):
        if label not in results:
            continue
        r = results[label]
        if r['fit'] is not None:
            xw, popt = r['fit']
            xx = np.linspace(xw[0], xw[-1], 60)
            ax.plot(xx, lorentzian(xx, *popt), '-', color=colors[label], lw=2.5,
                    label=f'{label} fit, x0={r["x0"]:.1f}nm', zorder=3)
        ax.plot(r['x0'], np.interp(r['x0'], x_nm, sig_f), '*', color=colors[label],
                ms=18, mec='k', mew=0.8, zorder=4)

    title = f'{wn}cm-1'
    if show_all_peaks:
        title += f' (n_peaks={len(peaks_idx)})'
    if wn in PROM_ABS_OVERRIDE:
        title += ' (custom prom)'
    ax.set_title(title, fontsize=12, fontweight='bold')
    ax.set_xlim(0, 1200)
    ax.set_xlabel('Distance from edge (nm)', fontweight='bold')
    ax.set_ylabel('Re signal (bg-sub, a.u.)', fontweight='bold')
    ax.legend(fontsize=8, loc='upper left')

    if 'peak2' in results and 'peak3' in results:
        x2, x3 = results['peak2']['x0'], results['peak3']['x0']
        q_2q = 2 * np.pi / ((x3 - x2) * 1e-7)  # nm -> cm, q in cm^-1
        ax.text(0.97, 0.95, f'2q = {q_2q/1e5:.2f}$\\times10^5$cm$^{{-1}}$\n'
                             f'q = {q_2q/2/1e5:.2f}$\\times10^5$cm$^{{-1}}$',
                transform=ax.transAxes, ha='right', va='top', fontsize=10, fontweight='bold',
                bbox=dict(boxstyle='round', fc='#eafaf1', ec='#4daf4a'))
        return q_2q / 2
    return np.nan


# --- panel 1: 860-911cm-1 (6 panels, 2x3) ---
fig1, axes1 = plt.subplots(2, 3, figsize=(18, 10))
q_low = {}
for ax, wn in zip(axes1.flat, WN_LOW):
    x_f, sig_f = load_signal(wn)
    q_low[wn] = panel(ax, wn, x_f, sig_f, show_all_peaks=False)
fig1.tight_layout()
out1 = f'{SAVE_DIR}/lorentzian_fit_realspace_preview_860to911_lp2.png'
fig1.savefig(out1, dpi=150, bbox_inches='tight')
plt.close(fig1)
print(f'Saved {out1}')

# --- panel 2: 920-1000cm-1 (9 panels, 3x3) ---
fig2, axes2 = plt.subplots(3, 3, figsize=(18, 14))
q_high = {}
for ax, wn in zip(axes2.flat, WN_HIGH):
    x_f, sig_f = load_signal(wn)
    q_high[wn] = panel(ax, wn, x_f, sig_f, show_all_peaks=True)
fig2.tight_layout()
out2 = f'{SAVE_DIR}/lorentzian_fit_realspace_preview_920to1000_lp2.png'
fig2.savefig(out2, dpi=150, bbox_inches='tight')
plt.close(fig2)
print(f'Saved {out2}')

print('\nq (peak2-peak3 / 2), 1e5 cm^-1:')
for wn, q in {**q_low, **q_high}.items():
    print(f'  {wn}: {q/1e5:.3f}' if not np.isnan(q) else f'  {wn}: n/a')
