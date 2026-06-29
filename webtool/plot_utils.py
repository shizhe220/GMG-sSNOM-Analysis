"""Plotly figure builders shared across webtool steps. Mirrors the
established style (Sky for amp, bwr_r for phase, gray for z; bold labels)
from snippet.linecut_extraction / nanoftir_Shizhe, just re-targeted at Plotly
so maps/spectra can be clicked on.
"""
import sys
sys.path.append('/Users/shizhe/envsetting')

import numpy as np
import plotly.graph_objects as go
import matplotlib.cm as mcm
from snippet import Sky


def mpl_cmap_to_plotly(cmap, n=256):
    if isinstance(cmap, str):
        cmap = mcm.get_cmap(cmap)
    colors = cmap(np.linspace(0, 1, n))
    return [[i / (n - 1), f'rgb({int(r*255)},{int(g*255)},{int(b*255)})']
            for i, (r, g, b, a) in enumerate(colors)]


_AMP_SCALE = mpl_cmap_to_plotly(Sky)
_PHASE_SCALE = mpl_cmap_to_plotly('bwr_r')
_Z_SCALE = mpl_cmap_to_plotly('gray')
CHANNEL_SCALE = {'amp': _AMP_SCALE, 'phase': _PHASE_SCALE, 'z': _Z_SCALE}
CHANNEL_LABEL = {'amp': 'Amplitude (a.u.)', 'phase': 'Phase (a.u.)', 'z': 'Height (nm)'}


def make_map_figure(data2d, pixelsize_um, channel, title='', vmin=None, vmax=None,
                     height=380):
    """2D map as a Plotly heatmap in physical (um) coordinates -- click events
    from streamlit_plotly_events come back in these same um units."""
    ny, nx = data2d.shape
    x = (np.arange(nx) + 0.5) * pixelsize_um
    y = (np.arange(ny) + 0.5) * pixelsize_um

    fig = go.Figure(go.Heatmap(
        z=data2d, x=x, y=y,
        colorscale=CHANNEL_SCALE.get(channel, 'gray'),
        zmin=vmin, zmax=vmax,
        colorbar=dict(title=CHANNEL_LABEL.get(channel, '')),
    ))
    fig.update_yaxes(autorange='reversed', scaleanchor='x', scaleratio=1,
                      title='Y (um)', title_font=dict(size=13))
    fig.update_xaxes(title='X (um)', title_font=dict(size=13))
    fig.update_layout(title=dict(text=title, font=dict(size=14)),
                       height=height, margin=dict(l=10, r=10, t=40, b=10))
    return fig


def add_line_overlay(fig, center_um, angle_deg, radius_um, width_um=None, color='cyan'):
    """Draw the proposed/placed linecut on a map figure (in-place)."""
    ang = np.radians(angle_deg)
    x0, y0 = center_um
    x1 = x0 + radius_um * np.cos(ang)
    y1 = y0 + radius_um * np.sin(ang)
    fig.add_trace(go.Scatter(x=[x0, x1], y=[y0, y1], mode='lines+markers',
                              line=dict(color=color, width=4), marker=dict(size=8, color=color),
                              name='cut', showlegend=False))
    if width_um:
        # show the strip width as two parallel guide lines
        dx, dy = -np.sin(ang) * width_um / 2, np.cos(ang) * width_um / 2
        for sgn in (1, -1):
            fig.add_trace(go.Scatter(
                x=[x0 + sgn * dx, x1 + sgn * dx], y=[y0 + sgn * dy, y1 + sgn * dy],
                mode='lines', line=dict(color=color, width=1, dash='dot'),
                showlegend=False, hoverinfo='skip'))
    return fig


def add_points_overlay(fig, points_um, color='red', label_prefix='P'):
    if not points_um:
        return fig
    xs, ys = zip(*points_um)
    fig.add_trace(go.Scatter(
        x=xs, y=ys, mode='markers+text', marker=dict(size=10, color=color, symbol='x'),
        text=[f'{label_prefix}{i+1}' for i in range(len(xs))], textposition='top center',
        showlegend=False))
    return fig


def make_waterfall_figure(linecut_store, align_dict, channel='amp', stack_offset=1.0,
                           normalize='mean', cmap_name='coolwarm_r', highlight_wn=None,
                           height=550):
    """Single-channel stacked waterfall, click-to-align friendly (one trace per wn,
    so plotly_events can report which curve+x was clicked)."""
    from data_io import get_wn

    wns = sorted(linecut_store.keys(), key=get_wn)
    fig = go.Figure()
    if not wns:
        fig.update_layout(height=height, title='No linecuts extracted yet')
        return fig

    n = len(wns)
    cmap = mcm.get_cmap(cmap_name)
    colors = [f'rgb({int(r*255)},{int(g*255)},{int(b*255)})'
              for r, g, b, a in cmap(np.linspace(0, 1, max(n, 2)))][:n]

    def _norm(arr):
        arr = np.asarray(arr, float)
        if normalize == 'mean':
            return arr - np.nanmean(arr)
        if normalize == 'minmax':
            lo, hi = np.nanmin(arr), np.nanmax(arr)
            return (arr - lo) / (hi - lo + 1e-12)
        return arr

    for i, wn in enumerate(wns):
        entry = linecut_store[wn]
        d = np.asarray(entry['distance_nm']) - align_dict.get(wn, 0.0)
        y = _norm(entry[channel]) + i * stack_offset
        is_hl = (wn == highlight_wn)
        fig.add_trace(go.Scatter(
            x=d, y=y, mode='lines', name=wn,
            line=dict(color=colors[i], width=4 if is_hl else 2),
            opacity=1.0 if (highlight_wn is None or is_hl) else 0.35,
        ))
    fig.add_vline(x=0, line=dict(color='gray', dash='dash'))
    fig.update_layout(
        height=height, showlegend=True,
        xaxis_title='Distance (nm)', yaxis_title=CHANNEL_LABEL.get(channel, ''),
        font=dict(size=13),
        margin=dict(l=10, r=10, t=30, b=10),
    )
    return fig
