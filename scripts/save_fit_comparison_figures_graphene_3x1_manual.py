"""
Generate the CHT-vs-Hankel-vs-1-sqrtx fit-comparison figures (combined single-
panel overlay + a/b/c/d panel breakdown) for all 15 wavenumbers of the
graphene_3x1_manual dataset, using the parameters currently tuned in
fitting_pipeline_graphene_3x1_manual.ipynb. Saves to
figures/graphene_3x1_manual/fit_comparison/{wn}_{combined,panels}.png.

Sibling of scripts/save_fit_results_graphene_3x1_manual.py -- reuses the same
notebook-anchor parameter-parsing logic, but calls fit_and_plot_cht/
compare_cavity_models directly (not run_wn_comparison) so the raw per-method
result dicts are available for plot_fit_comparison_combined/panels. Panel (a)
(map + linecut ROI) is read from data/graphene_3x1_manual/linecut_geometry.csv;
falls back to no panel (a) if a wn's row is missing.
"""
import json
import re
import os
import sys
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append('.')
sys.path.append('/Users/shizhe/envsetting')
import matplotlib
matplotlib.use('Agg')
import pandas as pd
import nanoftir_Shizhe as nanof
from loadnpz import loadnpz

NOTEBOOK_PATH = 'fitting_pipeline_graphene_3x1_manual.ipynb'
DATA_DIR = 'data/graphene_3x1_manual'
SAVE_DIR = 'figures/graphene_3x1_manual/fit_comparison'
MAP_NPZ_DIR = 'data/processed_3x1um'
GEOM_CSV = f'{DATA_DIR}/linecut_geometry.csv'

os.makedirs(SAVE_DIR, exist_ok=True)

# graphene_3x1_manual's CSVs were already aligned during manual extraction --
# no further per-wn shift needed here (matches save_fit_results_*.py).
align_dict = {f'{wn}cm-1': 0 for wn in
              (860, 870, 880, 890, 900, 911, 920, 930, 941, 950, 960, 970, 980, 991, 1000)}

# ---- 1. Parse the current (user-tuned) per-wn parameters straight out of the notebook ----
with open(NOTEBOOK_PATH) as f:
    nb = json.load(f)

select_idx = {}
for i, c in enumerate(nb['cells']):
    m = re.search(r"target_wn\s*=\s*'(\d+cm-1)'", ''.join(c['source']))
    if m:
        select_idx[int(m.group(1).replace('cm-1', ''))] = i
cht_idx = {i + 2: wn for wn, i in select_idx.items()}
rs_idx = {i + 4: wn for wn, i in select_idx.items()}

def parse_cht(idx):
    src = ''.join(nb['cells'][idx]['source'])
    return dict(
        x_start_cht=float(re.search(r'x_start_cht\s*=\s*([\d.]+)', src).group(1)),
        L_cutoff_cht=float(re.search(r'L_cutoff\s*=\s*([\d.]+)', src).group(1)),
        k_fit_range_cm=tuple(float(v) for v in re.search(r'k_fit_range_cm\s*=\s*\(([^)]+)\)', src).group(1).split(',')),
        k_linked_guess_cm=float(re.search(r'k_linked_guess_cm\s*=\s*([\d.]+)', src).group(1)),
    )

def parse_rs(idx):
    src = ''.join(nb['cells'][idx]['source'])
    return dict(
        xr_range_rs=tuple(float(v) for v in re.search(r'xr_range\s*=\s*\(([^)]+)\)', src).group(1).split(',')),
        lam0_guess_um=float(re.search(r'lam0_guess_um\s*=\s*([\d.]+)', src).group(1)),
    )

wn_params = {}
for wn_val in sorted(cht_idx.values(), reverse=True):
    wn = f'{wn_val}cm-1'
    cht_cell = [k for k, v in cht_idx.items() if v == wn_val][0]
    rs_cell = [k for k, v in rs_idx.items() if v == wn_val][0]
    wn_params[wn] = dict(**parse_cht(cht_cell), **parse_rs(rs_cell))

# ---- 2. Load linecut geometry (for panel (a)) ----
geom_df = pd.read_csv(GEOM_CSV) if os.path.exists(GEOM_CSV) else None

# ---- 3. Per-wn: fit + build comparison figures ----
for wn, p in wn_params.items():
    print(f'Running {wn} ...')
    wn_val = int(wn.replace('cm-1', ''))

    loaded = nanof.load_aligned_wn_signal(wn, align_dict[wn], data_dir=DATA_DIR, L_cutoff_fit=1.5)
    x_f, sig_f = loaded['x_f'], loaded['sig_f']

    cht_results, fig_cht = nanof.fit_and_plot_cht(
        x_f, sig_f, wn, x_start_cht=p['x_start_cht'], L_cutoff_cht=p['L_cutoff_cht'],
        k_fit_range_cm=p['k_fit_range_cm'], k_linked_guess_cm=p['k_linked_guess_cm'])

    amplp = pd.DataFrame({'distance_um': x_f, f'{wn}_O3A': sig_f})
    outs, fig_rs, _ = nanof.compare_cavity_models(
        amplp, f'{wn}_O3A', xr=p['xr_range_rs'], yc_um=1.9, fit_yc=False, edges='single',
        prefactors=('hankel', '1/sqrtx'), win=3, prom=0.01, lam0_guess=p['lam0_guess_um'],
        ylim=(sig_f.min() * 1.5, sig_f.max() * 0.5), figsize=(8, 6))

    map_panel = None
    if geom_df is not None:
        row = geom_df[geom_df['wn'] == wn]
        if len(row):
            row = row.iloc[0]
            scan_map = loadnpz(f'{MAP_NPZ_DIR}/{wn}_AVG.npz')
            map_panel = dict(
                amp_data=scan_map['O3A'], pixelsize_um=scan_map['x_pixelsize_nm'] / 1000.0,
                center_um=(row['center_x_um'], row['center_y_um']),
                angle_deg=row['angle_deg'], radius_um=row['radius_um'],
                rect_width_um=row['rect_width_um'], title='O3A + linecut ROI',
            )

    fig1, curves = nanof.plot_fit_comparison_combined(
        x_f, sig_f, cht_results, outs['hankel'], outs['1/sqrtx'], wn,
        save_path=f'{SAVE_DIR}/{wn}_combined.png')
    fig2, _ = nanof.plot_fit_comparison_panels(
        x_f, sig_f, cht_results, outs['hankel'], outs['1/sqrtx'], wn,
        map_panel=map_panel, save_path=f'{SAVE_DIR}/{wn}_panels.png')

    import matplotlib.pyplot as plt
    plt.close(fig_cht); plt.close(fig_rs); plt.close(fig1); plt.close(fig2)

    print(f"  {wn}: CHT={curves['cht']['lambda_p_nm']:.1f}nm  "
          f"Hankel={curves['hankel']['lambda_p_nm']:.1f}nm  "
          f"1/sqrtx={curves['1/sqrtx']['lambda_p_nm']:.1f}nm  "
          f"map_panel={'yes' if map_panel else 'no'}")

print(f'\nDone. Figures saved to {SAVE_DIR}/')
