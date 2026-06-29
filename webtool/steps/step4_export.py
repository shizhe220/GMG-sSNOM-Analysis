"""Step 4: export linecuts (CSV) + fit results (pickle/CSV) and plot q vs wn."""
import os
import pickle
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt

from data_io import get_wn


def render():
    st.header("4. Export & Summary")

    store = st.session_state['linecut_store']
    if not store:
        st.info("Nothing extracted yet.")
        return

    out_dir = st.text_input("Output folder", value='data/graphene_3x1_manual')
    use_aligned = st.checkbox("Bake the Step-2 align offset into the exported distance_nm", value=True)

    if st.button("Export linecut CSVs", type='primary'):
        os.makedirs(out_dir, exist_ok=True)
        for wn, entry in store.items():
            shift = st.session_state['align_dict'].get(wn, 0.0) if use_aligned else 0.0
            df = pd.DataFrame({
                'distance_nm': np.asarray(entry['distance_nm']) - shift,
                'O3A': entry['amp'], 'O3P': entry['phase'],
                'Z_nm': entry['z'], 'Z_nm_corrected': np.asarray(entry['z']) - np.nanmin(entry['z']),
            })
            df.to_csv(os.path.join(out_dir, f'{wn}_AVG_lp1.csv'), index=False)
        st.success(f"Wrote {len(store)} CSVs to {out_dir}")

    results = st.session_state['fit_results']
    if not results:
        st.info("No fit results yet (Step 3) -- export the comparison table/pickle once you have some.")
        return

    rows = []
    for wn in sorted(results.keys(), key=get_wn):
        r = results[wn]
        row = {'wn': wn}
        for method in ('cht', 'hankel', '1/sqrtx', 'fft'):
            m = r.get(method, {})
            row[f'{method}_lambda_nm'] = m.get('lambda_p_nm')
            row[f'{method}_q_1e5cm-1'] = m.get('q_p_1e5cm-1', m.get('q_p_1e5cm_1'))
            row[f'{method}_damping'] = m.get('damping')
        rows.append(row)
    df_summary = pd.DataFrame(rows)
    st.subheader("Summary table")
    st.dataframe(df_summary)

    pkl_path = st.text_input("Pickle path", value='data/webtool_fit_results.pkl')
    csv_path = st.text_input("Summary CSV path", value='data/webtool_fit_summary.csv')
    if st.button("Save results (pickle + CSV)"):
        os.makedirs(os.path.dirname(pkl_path) or '.', exist_ok=True)
        with open(pkl_path, 'wb') as f:
            pickle.dump(results, f)
        df_summary.to_csv(csv_path, index=False)
        st.success(f"Saved {pkl_path} and {csv_path}")

    st.subheader("q vs wavenumber")
    wns_sorted = sorted(results.keys(), key=get_wn)
    wn_vals = [get_wn(w) for w in wns_sorted]
    methods = {'cht': ('CHT', '#b2182b', 'o'), 'hankel': ('Hankel', '#1c7293', 's'),
               '1/sqrtx': ('1/sqrt(x)', '#e08214', '^'), 'fft': ('FFT', '#542788', 'D')}

    fig, ax = plt.subplots(figsize=(6, 6))
    for key, (label, color, marker) in methods.items():
        qs = []
        for wn in wns_sorted:
            m = results.get(wn, {}).get(key, {})
            q = m.get('q_p_1e5cm-1', m.get('q_p_1e5cm_1'))
            qs.append(q)
        ax.plot(qs, wn_vals, linestyle='None', marker=marker, ms=8, mfc='none', mec=color, mew=1.8, label=label)
    ax.set_xlabel('Momentum q_p (1e5 cm-1)', fontweight='bold')
    ax.set_ylabel('Wavenumber (cm-1)', fontweight='bold')
    ax.tick_params(direction='in', top=True, right=True)
    ax.legend(frameon=False, fontsize=11)
    ax.set_title('Extracted q_p vs Wavenumber', fontweight='bold')
    st.pyplot(fig)
