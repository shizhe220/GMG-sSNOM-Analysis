"""Step 1: pick a wavenumber, view its mapping, place a linecut by clicking
on the map (or typing numbers), optionally flatten Z, extract, add to the
running waterfall."""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append('/Users/shizhe/envsetting')

import numpy as np
import streamlit as st
from streamlit_plotly_events import plotly_events

from snippet.extract_linecut import radial_rect_profile
from snippet.linecut_extraction import suggest_flat_points, fit_plane_and_subtract
import plot_utils as pu
from data_io import get_wn


def _default_geom(scan, px):
    lp = scan['info'].get('lineprofile', {}).get('profiles', [])
    if lp:
        x0px, y0px = lp[0]['native_start_px']
        x1px, y1px = lp[0]['native_end_px']
        center = (x0px * px, y0px * px)
        angle = float(np.degrees(np.arctan2((y1px - y0px) * px, (x1px - x0px) * px)))
        radius = float(np.hypot((x1px - x0px) * px, (y1px - y0px) * px))
        width = lp[0].get('width_px', 5) * px
        return center, angle, radius, width
    ny, nx = scan['O3A'].shape
    return (0.1 * nx * px, 0.5 * ny * px), 0.0, 0.6 * nx * px, 5 * px


def render():
    st.header("1. Mapping & Linecut Extraction")

    if not st.session_state['scans']:
        st.info("Load a data folder from the sidebar first.")
        return

    wn = st.selectbox("Wavenumber", st.session_state['wn_list'], key='s1_wn')
    st.session_state['current_wn'] = wn
    scan = st.session_state['scans'][wn]
    px = st.session_state['pixelsize_um'][wn]

    geom = st.session_state['linecut_geom'].get(wn)
    if geom is None:
        center, angle, radius, width = _default_geom(scan, px)
        geom = dict(center_um=center, angle_deg=angle, radius_um=radius, rect_width_um=width)
        st.session_state['linecut_geom'][wn] = geom

    col_map, col_ctrl = st.columns([2, 1])

    with col_ctrl:
        st.subheader("Map channel")
        map_channel = st.radio("Show", ['amp', 'phase', 'z'], horizontal=True, key='s1_channel')

        st.subheader("Z-plane correction")
        correct_z = st.checkbox("Flatten Z before extracting", key='s1_correct_z')
        z_data = scan['Z']
        if correct_z:
            if st.button("Auto-suggest 3 flat points"):
                pts_px = suggest_flat_points(scan['Z'])
                st.session_state['plane_points_um'][wn] = [(x * px, y * px) for x, y in pts_px]
            pts_um = st.session_state['plane_points_um'].get(wn, [])
            st.caption("Click 3+ points on the Z map below to mark flat spots (or use auto-suggest).")
            if pts_um:
                st.write(", ".join(f"({x:.3f}, {y:.3f})" for x, y in pts_um))
                if st.button("Clear plane points"):
                    st.session_state['plane_points_um'][wn] = []
            if len(pts_um) >= 3:
                pts_px = [(x / px, y / px) for x, y in pts_um]
                z_data, _, _ = fit_plane_and_subtract(scan['Z'], pts_px)

        st.subheader("Cut geometry")
        cx = st.number_input("center x (um)", value=float(geom['center_um'][0]), step=0.01, format="%.4f")
        cy = st.number_input("center y (um)", value=float(geom['center_um'][1]), step=0.01, format="%.4f")
        angle = st.number_input("angle (deg)", value=float(geom['angle_deg']), step=1.0)
        radius = st.number_input("radius (um)", value=float(geom['radius_um']), step=0.05, format="%.4f")
        width = st.number_input("width (um)", value=float(geom['rect_width_um']), step=0.005, format="%.4f")
        geom = dict(center_um=(cx, cy), angle_deg=angle, radius_um=radius, rect_width_um=width)
        st.session_state['linecut_geom'][wn] = geom

        st.caption("Click twice on the map (start, then end of the cut) to set the geometry from the image instead.")
        if st.button("Clear click points"):
            st.session_state['click_points'] = []

        if st.button("Extract & add to waterfall", type='primary'):
            amp = scan['O3A']; phase = scan['O3P']
            kw = dict(center=geom['center_um'], radius=geom['radius_um'],
                      angle_deg=geom['angle_deg'], rect_width=geom['rect_width_um'], pixelsize=px)
            d_um, la = radial_rect_profile(amp, **kw)
            _, lp_ = radial_rect_profile(phase, **kw)
            _, lz = radial_rect_profile(z_data, **kw)
            st.session_state['linecut_store'][wn] = dict(
                distance_nm=d_um * 1000.0, amp=la, phase=lp_, z=lz)
            st.session_state['align_dict'].setdefault(wn, 0.0)
            st.success(f"Stored linecut for {wn} ({len(d_um)} points).")

    with col_map:
        display_data = z_data if (map_channel == 'z' and correct_z) else scan[{'amp': 'O3A', 'phase': 'O3P', 'z': 'Z'}[map_channel]]
        fig = pu.make_map_figure(display_data, px, map_channel, title=f"{wn} {map_channel}")
        pu.add_line_overlay(fig, geom['center_um'], geom['angle_deg'], geom['radius_um'], geom['rect_width_um'])
        if correct_z and map_channel == 'z':
            pu.add_points_overlay(fig, st.session_state['plane_points_um'].get(wn, []))

        clicked = plotly_events(fig, click_event=True, key=f's1_map_{wn}_{map_channel}')
        if clicked:
            x, y = clicked[0]['x'], clicked[0]['y']
            if correct_z and map_channel == 'z':
                pts = st.session_state['plane_points_um'].setdefault(wn, [])
                pts.append((x, y))
            else:
                pts = st.session_state['click_points']
                pts.append((x, y))
                if len(pts) >= 2:
                    (x0, y0), (x1, y1) = pts[-2], pts[-1]
                    new_angle = float(np.degrees(np.arctan2(y1 - y0, x1 - x0)))
                    new_radius = float(np.hypot(x1 - x0, y1 - y0))
                    st.session_state['linecut_geom'][wn] = dict(
                        center_um=(x0, y0), angle_deg=new_angle, radius_um=new_radius,
                        rect_width_um=geom['rect_width_um'])
                    st.session_state['click_points'] = []
            st.rerun()

        if wn in st.session_state['linecut_store']:
            st.caption("Current 1D extraction preview")
            entry = st.session_state['linecut_store'][wn]
            import plotly.graph_objects as go
            prof_fig = go.Figure()
            prof_fig.add_trace(go.Scatter(x=entry['distance_nm'], y=entry['amp'], name='amp'))
            prof_fig.update_layout(height=200, margin=dict(l=10, r=10, t=10, b=10),
                                    xaxis_title='Distance (nm)', yaxis_title='Amp (a.u.)')
            st.plotly_chart(prof_fig, use_container_width=True)

    if st.session_state['linecut_store']:
        st.subheader("Running waterfall (all extracted so far)")
        wf_channel = st.radio("Waterfall channel", ['amp', 'phase', 'z'], horizontal=True, key='s1_wf_channel')
        wf_fig = pu.make_waterfall_figure(st.session_state['linecut_store'], st.session_state['align_dict'],
                                           channel=wf_channel, highlight_wn=wn)
        st.plotly_chart(wf_fig, use_container_width=True)
