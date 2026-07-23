"""
Shared crosscheck CSV for graphene_4x1_manual lp2 (interior), mirroring lp1's own
GMG#3_4x1um_GMGregion_leftedge_v2(07222026).csv exactly (same columns/order/convention,
see the sibling README for lp1). Values recomputed here (not hand-copied) from the
notebook's own live-tuned params for full precision and to keep this one script as the
single source of truth, matching the project's established practice for this file.

Column meanings -- see the sibling _README.txt for the full writeup. Key notes specific
to lp2, all flagged in the per-row `note` column (not a top comment -- see lp1's own
2026-07-22 lesson: a `#`-comment header breaks plain pd.read_csv() for anyone who
doesn't know to pass comment='#'):
- hankel_qp_cm-1 (single-wave): 5 of 15 wn (930,941,950,980,991) are bound-saturation
  artifacts of fit_cavity_prefactor_compare's default lam_bound_factor=(0.7,1.2) -- NOT
  genuine convergence (see log_and_conventions.md 2026-07-23 for the systematic check).
  2 more (911,960) were originally saturated too but have since been fixed by re-tuning
  lam0_guess (see hankel_2wave_qp_cm-1 below, and log entry). Do not trust the 5 flagged
  values as physical -- use hankel_2wave_qp_cm-1 instead for those wn (900-1000 already
  uses two-wave anyway).
- hankel_2wave_qp_cm-1 (900-1000 only): the two-wave model's "q" term, same convention
  as lp1's own column.
- fft_recommended_channel/_qp_cm-1/_fwhm_cm-1: amp channel for ALL 15 wn. (An earlier
  version of this script copied lp1's 991/1000-use-phase convention onto lp2 without
  actually checking lp2's own amp-vs-phase FWHM per wn -- caught by the user. Checked
  properly: lp2's best channel varies wn by wn with no clean pattern like lp1 had (e.g.
  1000cm-1's phase FWHM is actually *worse* than amp here, the opposite of lp1). Default
  to amp for all 15 until/unless a real per-wn channel study is done for lp2.)
- peak2to3_spacing_qp_cm-1: real-space peak2-peak3/2, all 15 wn, after several rounds
  of manual peak-position correction (see log entry) -- 930 and 1000 are genuine,
  unresolved outliers vs the other methods (flagged in their note), not a picking bug.
"""
import os
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys
sys.path.append('.')
sys.path.append('/Users/shizhe/envsetting')
import numpy as np
import pandas as pd
from scipy.signal import savgol_filter
import nanoftir_Shizhe as nanof

SAVE_DIR = 'figures/rp_dispersion/graphene_4x1_manual_lp2_hankel_fft_compare_v1'
os.makedirs(SAVE_DIR, exist_ok=True)
CSV_PATH = f'{SAVE_DIR}/GMG#3_4x1um_GMGregion_rightedge_v1(07232026).csv'
README_PATH = f'{SAVE_DIR}/GMG#3_4x1um_GMGregion_rightedge_v1(07232026)_README.txt'

WNS = [860, 870, 880, 890, 900, 911, 920, 930, 941, 950, 960, 970, 980, 991, 1000]

# --- single-wave real-space params, live from fitting_pipeline_graphene_4x1_lp2.ipynb
# (911/960 use the 2026-07-23 fix: lam0_guess retuned to escape bound-saturation; 960
# additionally needs a widened lam_bound_factor -- see log entry) ---
SINGLE_PARAMS = {
    860:  dict(xr=(0.22, 1.0), lam0=0.647),
    870:  dict(xr=(0.22, 1.2), lam0=0.64),
    880:  dict(xr=(0.22, 1.2), lam0=0.56),
    890:  dict(xr=(0.22, 1.2), lam0=0.52),
    900:  dict(xr=(0.22, 1.2), lam0=0.49),
    911:  dict(xr=(0.11, 1.2), lam0=0.405),
    920:  dict(xr=(0.14, 1.2), lam0=0.42),
    930:  dict(xr=(0.06, 1.2), lam0=0.4),
    941:  dict(xr=(0.04, 1.2), lam0=0.4),
    950:  dict(xr=(0.06, 1.2), lam0=0.28),
    960:  dict(xr=(0.25, 1.2), lam0=0.375, lam_bound_factor=(0.85, 1.15)),
    970:  dict(xr=(0.04, 1.2), lam0=0.4),
    980:  dict(xr=(0.02, 1.2), lam0=0.4),
    991:  dict(xr=(0.02, 1.5), lam0=0.3),
    1000: dict(xr=(0.1, 1.2), lam0=0.2),
}
BOUND_SATURATED = {930, 941, 950, 980, 991}  # confirmed no genuine minimum nearby (see log)

# --- two-wave params, live from the notebook ---
TWOWAVE_PARAMS = {
    900:  dict(xr=(0.24, 1.5), q_guess=[1.3, 2.6]),
    911:  dict(xr=(0.12, 1.5), q_guess=[1.5, 3]),
    920:  dict(xr=(0.02, 1.2), q_guess=[1.2, 3.0]),
    930:  dict(xr=(0.06, 1.2), q_guess=[2.1, 4.2]),
    941:  dict(xr=(0.04, 1.2), q_guess=[2, 4]),
    950:  dict(xr=(0.06, 1.2), q_guess=[2.2, 4.5]),
    960:  dict(xr=(0.02, 1.2), q_guess=[2.5, 5]),
    970:  dict(xr=(0.04, 1.2), q_guess=[2.4, 3.6]),
    980:  dict(xr=(0.02, 1.2), q_guess=[2, 4.5]),
    991:  dict(xr=(0.02, 1.5), q_guess=[2.5, 5]),
    1000: dict(xr=(0.12, 1.2), q_guess=[3.1, 6.3]),
}

# --- FFT params, live from the notebook ---
FFT_PARAMS = {
    860: dict(xr=(0.22, 1), q_guess=[2], search_window=0.3, fit_window=0.3),
    870: dict(xr=(0.22, 1.2), q_guess=[2], search_window=0.3, fit_window=0.5),
    880: dict(xr=(0.22, 1.2), q_guess=[2.5], search_window=0.3, fit_window=0.3),
    890: dict(xr=(0.22, 1.2), q_guess=[2.6], search_window=0.3, fit_window=0.3),
    900: dict(xr=(0.22, 1.2), q_guess=[2.8], search_window=0.3, fit_window=0.5),
    911: dict(xr=(0.11, 1.2), q_guess=[1.6, 3], search_window=0.3, fit_window=0.3),
    920: dict(xr=(0.14, 1.2), q_guess=[1.6, 3], search_window=0.3, fit_window=0.5),
    930: dict(xr=(0.06, 1.2), q_guess=[2, 4.70], search_window=0.3, fit_window=0.5),
    941: dict(xr=(0.04, 1.2), q_guess=[2, 4], search_window=0.3, fit_window=0.5),
    950: dict(xr=(0.06, 1.2), q_guess=[2, 4], search_window=0.3, fit_window=0.5),
    960: dict(xr=(0.25, 1.2), q_guess=[2.5, 5.00], search_window=0.3, fit_window=0.5),
    970: dict(xr=(0.04, 1.2), q_guess=[2.6, 3.6], search_window=0.3, fit_window=0.16),
    980: dict(xr=(0.02, 1.2), q_guess=[2, 4.00], search_window=0.3, fit_window=0.5),
    991: dict(xr=(0.02, 1.2), q_guess=[2.5, 5], search_window=0.3, fit_window=0.2),
    1000: dict(xr=(0.12, 1.2), q_guess=[3, 6.60], search_window=0.3, fit_window=0.2),
}
FFT_PHASE_CHANNEL = set()  # amp for all 15 -- see README/docstring for why this is now empty

# --- peak2-peak3 spacing, final values after manual correction (see log entry) ---
PEAK23_Q = {
    860: 113400, 870: 108900, 880: 103300, 890: 113400, 900: 165700, 911: 196700,
    920: 170400, 930: 104700, 941: 214200, 950: 191600, 960: 197700, 970: 223000,
    980: 111900, 991: 127100, 1000: 169500,
}
PEAK23_OUTLIER = {930, 1000}  # genuine, unresolved disagreement vs other methods

rows = []
for wn in WNS:
    wn_str = f'{wn}cm-1'
    loaded = nanof.load_aligned_wn_signal(wn_str, 0, data_dir='data/graphene_4x1_manual', L_cutoff_fit=1.5, lp='lp2')
    x_f, sig_f = loaded['x_f'], loaded['sig_f']

    # single-wave hankel + 1/sqrtx
    sp = SINGLE_PARAMS[wn]
    kw = dict(lam_bound_factor=sp['lam_bound_factor']) if 'lam_bound_factor' in sp else {}
    out_h = nanof.fit_cavity_prefactor_compare(x_f, sig_f, xr=sp['xr'], yc_um=1.9, fit_yc=False,
                                                edges='single', prefactor='hankel', win=3, prom=0.01,
                                                lam0_guess=sp['lam0'], **kw)
    out_s = nanof.fit_cavity_prefactor_compare(x_f, sig_f, xr=sp['xr'], yc_um=1.9, fit_yc=False,
                                                edges='single', prefactor='1/sqrtx', win=3, prom=0.01,
                                                lam0_guess=sp['lam0'], **kw)
    hankel_qp = out_h['derived']['q_cm^-1']
    sqrtx_qp = out_s['derived']['q_cm^-1']

    # two-wave hankel
    hankel_2wave_qp = np.nan
    if wn in TWOWAVE_PARAMS:
        tp = TWOWAVE_PARAMS[wn]
        out_2w = nanof.fit_hankel_two_wave(x_f, sig_f, xr=tp['xr'], q_guess=tp['q_guess'],
                                            share_q_imag=False, q_imag0=0.5, win=3, prom=0.01, verbose=False)
        hankel_2wave_qp = out_2w['derived']['q_cm^-1']

    # FFT, all 3 channels
    fp = FFT_PARAMS[wn]
    df_target = loaded['df_target']
    amp_raw, phase_raw = df_target['O3A'], df_target['O3P']
    wlen = min(41, len(amp_raw) if len(amp_raw) % 2 != 0 else len(amp_raw) - 1)
    amp_osc = amp_raw - savgol_filter(amp_raw, window_length=wlen, polyorder=2)
    fft_out = nanof.plot_channel_fft(
        df_target['distance_um'], amp_osc, phase_raw, label='O3', wn=wn_str,
        xr=fp['xr'], q_range=(0, 10), window='hann', pad_factor=3.0, peak_method='lorentzian',
        q_guess=fp['q_guess'], search_window=fp['search_window'], fit_window=fp['fit_window'], plot=False)
    fft_vals = {}
    for ch in ('amp', 'phase', 'complex'):
        pk = fft_out[f'peaks_{ch}']['peaks']
        fw = fft_out[f'peaks_{ch}']['fwhm']
        q = round(pk[-1] / 10.0 / 2 * 1e5, 1) if pk and pk[-1] is not None else np.nan
        fwhm = round(fw[-1] / 10.0 / 2 * 1e5, 1) if fw and fw[-1] is not None else np.nan
        fft_vals[ch] = (q, fwhm)

    rec_channel = 'phase' if wn in FFT_PHASE_CHANNEL else 'amp'
    rec_q, rec_fwhm = fft_vals[rec_channel]

    note_parts = []
    if wn in TWOWAVE_PARAMS:
        note_parts.append('uses two-wave hankel fit')
    if wn in BOUND_SATURATED:
        note_parts.append('hankel_qp_cm-1 is bound-saturated, NOT genuine -- use hankel_2wave_qp_cm-1')
    if wn in FFT_PHASE_CHANNEL:
        note_parts.append('fft uses phase ch. (amp too noisy)')
    if wn in PEAK23_OUTLIER:
        note_parts.append('peak2to3_spacing is a genuine outlier vs other methods (not a picking error)')

    rows.append(dict(**{
        'freq_cm-1': wn,
        'hankel_qp_cm-1': round(hankel_qp),
        'sqrtx_qp_cm-1': round(sqrtx_qp),
        'fft_amp_qp_cm-1': fft_vals['amp'][0],
        'fft_amp_fwhm_cm-1': fft_vals['amp'][1],
        'fft_phase_qp_cm-1': fft_vals['phase'][0],
        'fft_phase_fwhm_cm-1': fft_vals['phase'][1],
        'fft_complex_qp_cm-1': fft_vals['complex'][0],
        'fft_complex_fwhm_cm-1': fft_vals['complex'][1],
        'peak2to3_spacing_qp_cm-1': PEAK23_Q[wn],
        'hankel_2wave_qp_cm-1': round(hankel_2wave_qp) if not np.isnan(hankel_2wave_qp) else np.nan,
        'fft_recommended_channel': rec_channel,
        'fft_recommended_qp_cm-1': rec_q,
        'fft_recommended_fwhm_cm-1': rec_fwhm,
        'note': '; '.join(note_parts),
    }))
    print(f'{wn_str}: hankel={hankel_qp:.0f}  sqrtx={sqrtx_qp:.0f}  2wave={hankel_2wave_qp:.0f}  '
          f'fft_rec({rec_channel})={rec_q}  peak23={PEAK23_Q[wn]}')

df = pd.DataFrame(rows)
df.to_csv(CSV_PATH, index=False)
print(f'\nSaved {CSV_PATH}')

readme = """GMG#3_4x1um_GMGregion_rightedge_v1(07232026).csv -- column notes (2026-07-23, internal use)

graphene_4x1_manual, lp2 (right edge / interior), wn = 860-1000 cm-1. All q_p columns
are in raw cm^-1 (not this project's internal 1e5 cm^-1 convention). Sibling of lp1's
own GMG#3_4x1um_GMGregion_leftedge_v2(07222026).csv -- same columns/convention.

hankel_qp_cm-1, sqrtx_qp_cm-1
    Single-wave real-space fit (Hankel / 1-sqrt(x)). CAUTION: 5 of 15 wn
    (930,941,950,980,991) are bound-saturation artifacts of the fit's default
    lam_bound_factor=(0.7,1.2) -- NOT genuine convergence (verified via a wide
    lam_bound_factor search: no genuine local minimum exists near the other methods'
    consensus for these 5). Flagged in each row's `note`. Two more (911,960) were
    originally saturated too but were fixed by re-tuning lam0_guess (911) and lam0_guess
    + a widened lam_bound_factor just for that wn (960) -- both now genuine.

hankel_2wave_qp_cm-1
    Two-wave Hankel fit (edge-launched "q" term + tip-launched "2q" term). Filled for
    wn=900-1000 only -- use this instead of hankel_qp_cm-1 for that range (matches lp1's
    convention exactly).

fft_amp_qp_cm-1, fft_phase_qp_cm-1, fft_complex_qp_cm-1 (+ _fwhm_cm-1 each)
    FFT peak position (Lorentzian-fit) from each of the 3 channels. All 3 always
    computed and stored regardless of which is "recommended".

fft_recommended_channel, fft_recommended_qp_cm-1, fft_recommended_fwhm_cm-1
    "amp" for all 15 wn -- unlike lp1, no per-wn amp-vs-phase channel study has been
    done for lp2 yet, so this just defaults to the project's standard channel rather
    than guessing. (lp1's convention of using phase for 991/1000 does NOT carry over
    here -- checked, and lp2's per-channel FWHM doesn't follow the same pattern; e.g.
    1000cm-1's phase-channel FWHM is actually worse than amp's for lp2.)

peak2to3_spacing_qp_cm-1
    Real-space peak2-peak3 spacing / 2, independent of any fitting model. Filled for
    all 15 wn after several rounds of manual peak-position review (see
    log_and_conventions.md for the full correction history -- some wn needed the user's
    direct visual check to catch a mis-picked valley/shoulder). 930 and 1000cm-1 are
    genuine, unresolved outliers vs. hankel_2wave/fft (~40-50% lower) even after
    correction -- the real-space signal at those two wn has multiple close sub-peaks
    rather than one clean oscillation lobe, so this isn't a picking bug, just an open
    question about which sub-feature the spacing method should be reading. Flagged in
    their `note`.

note
    Short in-CSV flags for the caveats above -- kept brief since this file may be
    shared; see this README for the detailed version, for our own records.
"""
with open(README_PATH, 'w') as f:
    f.write(readme)
print(f'Saved {README_PATH}')
