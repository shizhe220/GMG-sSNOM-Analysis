"""Export a CSV summary from data/fit_results.pkl: per-wn wavelength from CHT,
Hankel, and 1/sqrtx, the CHT momentum fit window and real-space distance fit
range used, and CHT's deviation (%) from each real-space method."""
import json
import re
import pickle
import pandas as pd

with open('data/fit_results.pkl', 'rb') as f:
    results = pickle.load(f)

with open('fitting_pipeline.ipynb') as f:
    nb = json.load(f)

cht_idx = {7: 860, 15: 870, 23: 880, 31: 890, 39: 900, 47: 911, 55: 920, 63: 930,
           71: 941, 79: 950, 87: 960, 95: 970, 103: 980, 111: 991, 119: 1000}

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
df.to_csv('data/cht_vs_realspace_wavelength_comparison.csv', index=False)
print(df.to_string(index=False))
print('\nSaved data/cht_vs_realspace_wavelength_comparison.csv')
