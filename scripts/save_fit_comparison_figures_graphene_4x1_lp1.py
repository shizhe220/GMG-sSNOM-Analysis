"""
Generate the CHT-vs-Hankel-vs-1-sqrtx fit-comparison figures (combined single-
panel overlay + a/b/c/d panel breakdown) for all 15 wavenumbers of the
graphene_4x1_manual lp1 (left edge) dataset, using the parameters currently
tuned in fitting_pipeline_graphene_4x1_lp1.ipynb. Saves to
figures/graphene_4x1_manual/lp1/fit_comparison/{wn}_{combined,panels}.png.

Sibling of scripts/save_fit_comparison_figures_graphene_3x1_manual.py --
same structure, but reads the combined multi-wn npz (loadnpz_combined) for
panel (a)'s map, and both load_aligned_wn_signal/linecut_geometry.csv lookups
are scoped to lp='lp1' (data/graphene_4x1_manual/ holds both lp1 and lp2 CSVs
for every wn, so this distinction matters -- see log_and_conventions.md
2026-07-08 for the bug this caused when it was missing).
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
from loadnpz_combined import loadnpz

LP = 'lp1'
NOTEBOOK_PATH = f'fitting_pipeline_graphene_4x1_{LP}.ipynb'
DATA_DIR = 'data/graphene_4x1_manual'
SAVE_DIR = f'figures/graphene_4x1_manual/{LP}/fit_comparison'
MAP_NPZ_PATH = 'data/processed_4x1um/860-1000_AVG.npz'
GEOM_CSV = f'{DATA_DIR}/linecut_geometry.csv'

os.makedirs(SAVE_DIR, exist_ok=True)

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

# ---- 2. Load linecut geometry (for panel (a)), scoped to this lp ----
geom_df = pd.read_csv(GEOM_CSV) if os.path.exists(GEOM_CSV) else None
if geom_df is not None:
    geom_df = geom_df[geom_df['lp'] == LP]

scan = loadnpz(MAP_NPZ_PATH)

# ---- 3. Per-wn: fit + build comparison figures ----
for wn, p in wn_params.items():
    print(f'Running {wn} ...')
    wn_val = int(wn.replace('cm-1', ''))

    loaded = nanof.load_aligned_wn_signal(wn, align_dict[wn], data_dir=DATA_DIR, L_cutoff_fit=1.5, lp=LP)
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
            map_panel = dict(
                amp_data=scan['maps'][wn_val]['O3A'], pixelsize_um=scan['pixelsize_nm'] / 1000.0,
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
