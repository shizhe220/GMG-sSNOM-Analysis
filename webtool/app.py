"""GMG linecut/fit web tool -- entry point.

Run with:
    streamlit run webtool/app.py

Workflow: load a folder of *_AVG.npz files -> per-wn mapping+linecut
extraction (click on the map to place the cut) -> align (click on the
waterfall) -> fit/FFT (click to set q_guess) -> export + q-vs-wn summary.
Every step reads/writes the same st.session_state, so you can jump back to
any wavenumber at any step without losing progress.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append('/Users/shizhe/envsetting')

import streamlit as st

from state import init_state, wn_status, STATUS_EMOJI
from data_io import load_folder
from steps import step1_linecut, step2_align, step3_fit, step4_export

st.set_page_config(page_title="GMG Linecut & Fit Tool", layout="wide")
init_state()

with st.sidebar:
    st.title("GMG Linecut & Fit Tool")

    st.session_state['data_dir'] = st.text_input("npz folder", value=st.session_state['data_dir'])
    if st.button("Load folder"):
        try:
            scans, wn_list = load_folder(st.session_state['data_dir'])
            st.session_state['scans'] = scans
            st.session_state['wn_list'] = wn_list
            for wn in wn_list:
                scan = scans[wn]
                px = scan.get('x_pixelsize_nm', scan.get('pixelsize_nm'))
                st.session_state['pixelsize_um'][wn] = (px or 15.0) / 1000.0
            st.success(f"Loaded {len(wn_list)} wavenumbers.")
        except Exception as e:
            st.error(str(e))

    st.divider()
    step = st.radio("Step", [
        "1. Mapping & Linecut", "2. Align", "3. Fit / FFT", "4. Export & Summary",
    ])

    if st.session_state['wn_list']:
        st.divider()
        st.caption("Progress")
        for wn in st.session_state['wn_list']:
            status = wn_status(wn)
            st.write(f"{STATUS_EMOJI[status]} {wn} -- {status}")

if step.startswith("1"):
    step1_linecut.render()
elif step.startswith("2"):
    step2_align.render()
elif step.startswith("3"):
    step3_fit.render()
else:
    step4_export.render()
