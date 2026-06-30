"""Export a CSV summary from data/fit_results_graphene_4x1_lp1.pkl: per-wn
wavelength from CHT, Hankel, and 1/sqrtx, the CHT momentum fit window and
real-space distance fit range used, and CHT's deviation (%) from each
real-space method.

Sibling of scripts/export_comparison_csv_graphene_3x1_manual.py, pointed at
the graphene_4x1_manual lp1 (left edge) dataset."""
import json
import re
import pickle
import os
import pandas as pd

os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

NOTEBOOK_PATH = 'fitting_pipeline_graphene_4x1_lp1.ipynb'
PICKLE_PATH = 'data/fit_results_graphene_4x1_lp1.pkl'
OUT_CSV = 'data/cht_vs_realspace_wavelength_comparison_graphene_4x1_lp1.csv'

with open(PICKLE_PATH, 'rb') as f:
    results = pickle.load(f)

with open(NOTEBOOK_PATH) as f:
    nb = json.load(f)

select_idx = {}
for i, c in enumerate(nb['cells']):
    m = re.search(r"target_wn\s*=\s*'(\d+cm-1)'", ''.join(c['source']))
    if m:
        select_idx[int(m.group(1).replace('cm-1', ''))] = i
cht_idx = {i + 2: wn for wn, i in select_idx.items()}

def parse_x_range(idx):
    src = ''.join(nb['cells'][idx]['source'])
    x_start = float(re.search(r'x_start_cht\s*=\s*([\d.]+)', src).group(1))
    l_cutoff = float(re.search(r'L_cutoff\s*=\s*([\d.]+)', src).group(1))
    return (x_start, l_cutoff)

rows = []
for idx, wn_val in sorted(cht_idx.items(), key=lambda kv: -kv[1]):
    wn = f'{wn_val}cm-1'
    r = results[wn]
    lam_cht, lam_hankel, lam_sqrtx = r['cht']['lambda_p_nm'], r['hankel']['lambda_p_nm'], r['sqrtx']['lambda_p_nm']
    x_range_um = parse_x_range(idx)

    rows.append(dict(
        wn_cm1=wn_val,
        lambda_cht_nm=round(lam_cht, 1),
        lambda_hankel_nm=round(lam_hankel, 1),
        lambda_sqrtx_nm=round(lam_sqrtx, 1),
        k_fit_range_1e5cm1=str(r['cht']['k_fit_range_cm']),
        x_fit_range_um=str(x_range_um),
        cht_vs_hankel_pct=round(100 * abs(lam_cht - lam_hankel) / lam_hankel, 1),
        cht_vs_sqrtx_pct=round(100 * abs(lam_cht - lam_sqrtx) / lam_sqrtx, 1),
    ))

df = pd.DataFrame(rows)
df.to_csv(OUT_CSV, index=False)
print(df.to_string(index=False))
print(f'\nSaved {OUT_CSV}')
