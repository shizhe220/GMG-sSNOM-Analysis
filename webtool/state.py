"""Session-state schema for the linecut/fit web tool. Everything lives in
st.session_state so any step can read/write any wn's data and you can jump
back and forth without losing progress."""
import streamlit as st


def init_state():
    defaults = {
        'data_dir': 'data/processed_3x1um',
        'scans': {},            # wn -> Scan (from loadnpz)
        'wn_list': [],          # sorted descending
        'pixelsize_um': {},     # wn -> float

        'plane_corrected_z': {},  # wn -> 2D array (only present once corrected)
        'plane_points_um': {},    # wn -> list of (x_um, y_um)

        'linecut_geom': {},     # wn -> dict(center_um, angle_deg, radius_um, rect_width_um)
        'linecut_store': {},    # wn -> dict(distance_nm, amp, phase, z)
        'align_dict': {},       # wn -> float (nm)

        'fit_params': {},       # wn -> dict of fit input params (per method)
        'fit_results': {},      # wn -> dict(cht=..., hankel=..., sqrtx=..., fft=...)

        'click_points': [],     # scratch: clicked (x,y) for the line-placement UI
        'plane_click_points': [],  # scratch: clicked (x,y) for plane-correction points
        'q_guess_clicks': [],   # scratch: clicked q for fit q_guess

        'current_wn': None,
    }
    for k, v in defaults.items():
        if k not in st.session_state:
            st.session_state[k] = v


def wn_status(wn):
    """Return a short status string for the sidebar/progress table."""
    has_cut = wn in st.session_state['linecut_store']
    has_align = wn in st.session_state['align_dict']
    has_fit = wn in st.session_state['fit_results'] and bool(st.session_state['fit_results'][wn])
    if has_fit:
        return 'fitted'
    if has_align:
        return 'aligned'
    if has_cut:
        return 'cut'
    return 'pending'


STATUS_EMOJI = {'pending': '⬜', 'cut': '🟨', 'aligned': '🟦', 'fitted': '🟩'}
