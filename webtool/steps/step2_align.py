"""Step 2: align extracted linecuts by clicking on the waterfall (or typing
a number)."""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import streamlit as st
from streamlit_plotly_events import plotly_events

import plot_utils as pu
from data_io import get_wn


def render():
    st.header("2. Align")

    store = st.session_state['linecut_store']
    if not store:
        st.info("Extract at least one linecut in Step 1 first.")
        return

    wns = sorted(store.keys(), key=get_wn)
    wn = st.selectbox("Wavenumber to align", wns, key='s2_wn')
    st.session_state['current_wn'] = wn

    channel = st.radio("Channel", ['amp', 'phase', 'z'], horizontal=True, key='s2_channel')
    current_align = st.session_state['align_dict'].get(wn, 0.0)

    col1, col2 = st.columns([1, 3])
    with col1:
        st.metric("Current align (nm)", f"{current_align:.1f}")
        new_val = st.number_input("Set align (nm)", value=float(current_align), step=1.0, key='s2_numeric')
        if st.button("Apply numeric value"):
            st.session_state['align_dict'][wn] = new_val
            st.rerun()
        if st.button("Reset to 0"):
            st.session_state['align_dict'][wn] = 0.0
            st.rerun()
        st.caption("Or click on the highlighted curve below at the point that should become distance=0. "
                   "That shifts this wn's align value by the clicked x (added to whatever it already is).")

    with col2:
        fig = pu.make_waterfall_figure(store, st.session_state['align_dict'], channel=channel,
                                        highlight_wn=wn, height=600)
        clicked = plotly_events(fig, click_event=True, key=f's2_wf_{wn}_{channel}')
        if clicked:
            x_click = clicked[0]['x']
            st.session_state['align_dict'][wn] = current_align + x_click
            st.rerun()
