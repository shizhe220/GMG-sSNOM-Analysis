"""Step 3: CHT / real-space (Hankel, 1/sqrtx) / FFT fitting per wavenumber.
Reuses nanoftir_Shizhe.py's existing fitting+plotting functions (matplotlib,
rendered via st.pyplot) so the math/visual style stays identical to the
notebook pipeline; the only new interactive bit is clicking on a |T(k)|
preview to set the CHT momentum guess.
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append('/Users/shizhe/envsetting')

import numpy as np
import pandas as pd
import streamlit as st
import plotly.graph_objects as go
from streamlit_plotly_events import plotly_events
from scipy.signal import savgol_filter

import nanoftir_Shizhe as nanoftir
from data_io import get_wn


def _get_xf_sigf(wn, l_cutoff_fit=1.5):
    entry = st.session_state['linecut_store'][wn]
    align = st.session_state['align_dict'].get(wn, 0.0)
    distance_um = (np.asarray(entry['distance_nm'], float) - align) / 1000.0
    y_mat = np.asarray(entry['amp'], float)
    window_len = min(41, len(y_mat) if len(y_mat) % 2 != 0 else len(y_mat) - 1)
    window_len = max(window_len, 5)
    y_bg = savgol_filter(y_mat, window_length=window_len, polyorder=2)
    y_osc = y_mat - y_bg
    mask = (distance_um >= 0) & (distance_um <= l_cutoff_fit)
    return distance_um[mask], y_osc[mask]


def render():
    st.header("3. Fit / FFT")

    store = st.session_state['linecut_store']
    aligned_wns = [wn for wn in store if wn in st.session_state['align_dict']]
    if not aligned_wns:
        st.info("Extract + align at least one wn (Steps 1-2) first.")
        return

    wns = sorted(aligned_wns, key=get_wn)
    wn = st.selectbox("Wavenumber", wns, key='s3_wn')
    st.session_state['current_wn'] = wn
    st.session_state['fit_params'].setdefault(wn, {})
    params = st.session_state['fit_params'][wn]

    x_f, sig_f = _get_xf_sigf(wn)

    tab_cht, tab_rs, tab_fft = st.tabs(["CHT", "Real-space (Hankel/1√x)", "FFT"])

    # ---------------- CHT ----------------
    with tab_cht:
        c1, c2 = st.columns([1, 2])
        with c1:
            x_start_cht = st.number_input("x_start_cht (um)", value=params.get('x_start_cht', 0.15), step=0.01, key='s3_xstart')
            l_cutoff_cht = st.number_input("L_cutoff_cht (um)", value=params.get('L_cutoff_cht', 1.2), step=0.05, key='s3_lcut')
            k_lo = st.number_input("k_fit_range low (1e5 cm-1)", value=params.get('k_fit_lo', 0.5), step=0.1, key='s3_klo')
            k_hi = st.number_input("k_fit_range high (1e5 cm-1)", value=params.get('k_fit_hi', 6.0), step=0.1, key='s3_khi')
            k_guess = st.number_input("k_linked_guess (1e5 cm-1) -- or click the chart", value=params.get('k_guess', 2.0), step=0.1, key='s3_kguess')
            st.caption("Click the |T(k)| preview to set the guess from the dominant peak.")

        with c2:
            mask_cht = x_f >= x_start_cht
            x_f_cht, sig_f_cht = x_f[mask_cht], sig_f[mask_cht]
            if len(x_f_cht) > 5:
                k_prev = np.linspace(0.1, 100, 400)
                _, T_prev = nanoftir.complex_hankel_transform(x_f_cht, sig_f_cht, L=l_cutoff_cht, k_array=k_prev)
                prev_fig = go.Figure(go.Scatter(x=k_prev / 10.0, y=np.abs(T_prev), mode='lines'))
                prev_fig.add_vrect(x0=k_lo, x1=k_hi, fillcolor='gray', opacity=0.2, line_width=0)
                prev_fig.add_vline(x=k_guess, line=dict(color='purple', dash='dash'))
                prev_fig.update_layout(height=280, margin=dict(l=10, r=10, t=20, b=10),
                                        xaxis_title='k (1e5 cm-1)', yaxis_title='|T(k)|')
                clicked = plotly_events(prev_fig, click_event=True, key=f's3_tk_{wn}')
                if clicked:
                    k_guess = float(clicked[0]['x'])
                    params['k_guess'] = k_guess
                    st.rerun()

        params.update(x_start_cht=x_start_cht, L_cutoff_cht=l_cutoff_cht,
                      k_fit_lo=k_lo, k_fit_hi=k_hi, k_guess=k_guess)

        if st.button("Run CHT fit", key='s3_run_cht'):
            results, fig = nanoftir.fit_and_plot_cht(
                x_f, sig_f, wn, x_start_cht=x_start_cht, L_cutoff_cht=l_cutoff_cht,
                k_fit_range_cm=(k_lo, k_hi), k_linked_guess_cm=k_guess)
            st.session_state['fit_results'].setdefault(wn, {})['cht'] = results
            st.pyplot(fig)
            st.json({k: (round(v, 4) if isinstance(v, float) else v) for k, v in results.items()
                     if k not in ('k_fit_range_cm',)})

    # ---------------- Real space ----------------
    with tab_rs:
        c1, c2 = st.columns([1, 2])
        with c1:
            xr_lo = st.number_input("xr low (um)", value=params.get('xr_lo', x_start_cht if 'x_start_cht' in params else 0.15), step=0.01, key='s3_xrlo')
            xr_hi = st.number_input("xr high (um)", value=params.get('xr_hi', 1.5), step=0.05, key='s3_xrhi')
            lam0_guess = st.number_input("lam0_guess (um)", value=params.get('lam0_guess', 0.3), step=0.01, format="%.3f", key='s3_lam0')
        params.update(xr_lo=xr_lo, xr_hi=xr_hi, lam0_guess=lam0_guess)

        with c2:
            st.caption("Run the fit to preview the real-space curves; q_guess isn't needed here -- "
                       "lam0_guess plus automatic peak-finding drive the initial guess instead.")

        if st.button("Run real-space fit", key='s3_run_rs'):
            amplp = pd.DataFrame({'distance_um': x_f, f'{wn}_O3A': sig_f})
            outs, fig, _ = nanoftir.compare_cavity_models(
                amplp, f'{wn}_O3A', xr=(xr_lo, xr_hi), yc_um=1.9, fit_yc=False, edges='single',
                prefactors=('hankel', '1/sqrtx'), win=3, prom=0.01, lam0_guess=lam0_guess,
                ylim=(sig_f.min() * 2, sig_f.max() * 1.5), figsize=(8, 6))
            results = st.session_state['fit_results'].setdefault(wn, {})
            for pf in ('hankel', '1/sqrtx'):
                p, d, met = outs[pf]['params'], outs[pf]['derived'], outs[pf]['metrics']
                damp_key = 'q_imag_um^-1' if pf == 'hankel' else 'alpha_env_um^-1'
                results[pf] = dict(
                    lambda_p_nm=p['lambda_p_um'] * 1000, q_p_1e5cm_1=d['q_cm^-1'] / 1e5,
                    damping=d['q_rad_per_um'] / p[damp_key] if p[damp_key] > 1e-9 else float('inf'),
                    rmse=met['rmse'], aic=met['aic'])
            st.pyplot(fig)

    # ---------------- FFT ----------------
    with tab_fft:
        c1, c2 = st.columns([1, 2])
        with c1:
            fft_xr_lo = st.number_input("FFT xr low (um)", value=params.get('fft_xr_lo', xr_lo if 'xr_lo' in params else 0.15), step=0.01, key='s3_fxrlo')
            fft_xr_hi = st.number_input("FFT xr high (um)", value=params.get('fft_xr_hi', 1.2), step=0.05, key='s3_fxrhi')
            q_guess_1 = st.number_input("q_guess (1e5 cm-1) -- or click the chart", value=params.get('fft_q1', 2.0), step=0.1, key='s3_fq1')
            q_guess_2 = st.number_input("2nd q_guess (optional, 0=skip)", value=params.get('fft_q2', 0.0), step=0.1, key='s3_fq2')

        entry = st.session_state['linecut_store'][wn]
        align = st.session_state['align_dict'].get(wn, 0.0)
        dist_full_um = (np.asarray(entry['distance_nm'], float) - align) / 1000.0
        amp_full = np.asarray(entry['amp'], float)
        phase_full = np.asarray(entry['phase'], float)
        window_len = min(41, len(amp_full) if len(amp_full) % 2 != 0 else len(amp_full) - 1)
        window_len = max(window_len, 5)
        amp_osc_full = amp_full - savgol_filter(amp_full, window_length=window_len, polyorder=2)

        with c2:
            mask_f = (dist_full_um >= fft_xr_lo) & (dist_full_um <= fft_xr_hi)
            if mask_f.sum() > 5:
                sig_seg = amp_osc_full[mask_f]
                n = len(sig_seg)
                dx = np.median(np.diff(dist_full_um[mask_f]))
                freqs = 2 * np.pi * np.fft.rfftfreq(n, d=dx) / 10.0  # -> 1e5 cm^-1
                spec = np.abs(np.fft.rfft(sig_seg - sig_seg.mean()))
                prev_fig = go.Figure(go.Scatter(x=freqs, y=spec, mode='lines'))
                prev_fig.add_vline(x=q_guess_1, line=dict(color='purple', dash='dash'))
                if q_guess_2 > 0:
                    prev_fig.add_vline(x=q_guess_2, line=dict(color='blue', dash='dash'))
                prev_fig.update_layout(height=280, margin=dict(l=10, r=10, t=20, b=10),
                                        xaxis_title='q (1e5 cm-1)', yaxis_title='FFT amp')
                clicked = plotly_events(prev_fig, click_event=True, key=f's3_fft_{wn}')
                if clicked:
                    params['fft_q1'] = float(clicked[0]['x'])
                    st.rerun()

        params.update(fft_xr_lo=fft_xr_lo, fft_xr_hi=fft_xr_hi, fft_q1=q_guess_1, fft_q2=q_guess_2)

        if st.button("Run FFT", key='s3_run_fft'):
            q_guess_list = [q_guess_1] + ([q_guess_2] if q_guess_2 > 0 else [])
            out = nanoftir.plot_channel_fft(
                pd.Series(dist_full_um), pd.Series(amp_osc_full), pd.Series(phase_full),
                label='O3', wn=wn, xr=(fft_xr_lo, fft_xr_hi), q_range=(0, 10),
                window='hann', pad_factor=3.0, q_guess=q_guess_list)
            fig = st.session_state.get('_last_mpl_fig')
            import matplotlib.pyplot as plt
            st.pyplot(plt.gcf())
            peaks = [p for p in out['peaks_complex']['peaks'] if p is not None]
            results = st.session_state['fit_results'].setdefault(wn, {})
            if peaks:
                q_fft = peaks[0] / 10.0
                results['fft'] = dict(q_p_1e5cm_1=q_fft, lambda_p_nm=(2 * np.pi / q_fft) * 100, damping=None)
            else:
                results['fft'] = dict(q_p_1e5cm_1=None, lambda_p_nm=None, damping=None)

    st.divider()
    if wn in st.session_state['fit_results']:
        st.subheader(f"Stored results for {wn}")
        st.json(st.session_state['fit_results'][wn])
