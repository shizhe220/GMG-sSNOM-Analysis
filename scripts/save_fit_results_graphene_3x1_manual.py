"""
Regenerate the cht/realspace/fft figures for all 15 wavenumbers using the
parameters currently tuned in fitting_pipeline_graphene_3x1_manual.ipynb, and
save a results dict (lambda_p, momentum q, damping per method) to
data/fit_results_graphene_3x1_manual.pkl.

Sibling of scripts/save_fit_results.py, pointed at the graphene_3x1_manual
dataset -- no "old window" comparison here since there's no tuning history
for this dataset yet (that comparison in the original script was specifically
about the 960-1000cm-1 k_fit_range regression discovered for graphene_3x1).
"""
import json
import re
import pickle
import os
import sys
# Run from anywhere -- chdir to the GMG repo root (this script's parent dir) so
# the relative paths below ('data/...', 'figures/...', notebook path) resolve
# the same way they did back when this script lived at the repo root.
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append('.')
sys.path.append('/Users/shizhe/envsetting')
import matplotlib
matplotlib.use('Agg')
import nanoftir_Shizhe as nanof

NOTEBOOK_PATH = 'fitting_pipeline_graphene_3x1_manual.ipynb'
DATA_DIR = 'data/graphene_3x1_manual'
SAVE_DIR = 'figures/graphene_3x1_manual'
PICKLE_PATH = 'data/fit_results_graphene_3x1_manual.pkl'

# graphene_3x1_manual's CSVs were already aligned during manual extraction
# (extract_and_plot's align_bg_loc_nm/realign_nm baked the shift into the
# exported distance_nm) -- no further per-wn shift needed here by default.
align_dict = {f'{wn}cm-1': 0 for wn in
              (860, 870, 880, 890, 900, 911, 920, 930, 941, 950, 960, 970, 980, 991, 1000)}

# ---- 1. Parse the current (user-tuned) per-wn parameters straight out of the notebook ----
with open(NOTEBOOK_PATH) as f:
    nb = json.load(f)

# Don't hardcode cell indices -- they shift whenever a cell gets inserted/removed
# anywhere earlier in the notebook. Instead, locate each wn's own
# "target_wn = '...'" select cell and use it as an anchor: the per-wn template
# always lays out select -> CHT (+2) -> real-space (+4) -> FFT (+6).
select_idx = {}
for i, c in enumerate(nb['cells']):
    m = re.search(r"target_wn\s*=\s*'(\d+cm-1)'", ''.join(c['source']))
    if m:
        select_idx[int(m.group(1).replace('cm-1', ''))] = i
cht_idx = {i + 2: wn for wn, i in select_idx.items()}
rs_idx = {i + 4: wn for wn, i in select_idx.items()}
fft_idx = {i + 6: wn for wn, i in select_idx.items()}

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

def parse_fft(idx):
    src = ''.join(nb['cells'][idx]['source'])
    q_guess = [float(v) for v in re.search(r'q_guess\s*=\s*\[([^\]]+)\]', src).group(1).split(',')]
    fft_xr = tuple(float(v) for v in re.search(r'xr=\(([^)]+)\)', src).group(1).split(','))
    return dict(fft_q_guess=q_guess, fft_xr=fft_xr)

wn_params = {}
for wn_val in sorted(cht_idx.values(), reverse=True):
    wn = f'{wn_val}cm-1'
    cht_cell = [k for k, v in cht_idx.items() if v == wn_val][0]
    rs_cell = [k for k, v in rs_idx.items() if v == wn_val][0]
    fft_cell = [k for k, v in fft_idx.items() if v == wn_val][0]
    wn_params[wn] = dict(**parse_cht(cht_cell), **parse_rs(rs_cell), **parse_fft(fft_cell))

# ---- 2. Run the current notebook parameters, saving all 3 figures per wn ----
all_results = {}
for wn, p in wn_params.items():
    print(f'Running {wn} ...')
    results, _ = nanof.run_wn_comparison(
        wn, align_shift_nm=align_dict[wn], k_linked_guess_cm=p['k_linked_guess_cm'],
        k_fit_range_cm=p['k_fit_range_cm'], x_start_cht=p['x_start_cht'], L_cutoff_cht=p['L_cutoff_cht'],
        lam0_guess_um=p['lam0_guess_um'], xr_range_rs=p['xr_range_rs'],
        fft_xr=p['fft_xr'], fft_q_guess=p['fft_q_guess'],
        data_dir=DATA_DIR, save_dir=SAVE_DIR, show=False)

    all_results[wn] = {
        'cht': {
            'lambda_p_nm': results['lambda_p_nm'],
            'q_p_1e5cm-1': results['q_re'] / 10.0,
            'damping': results['damping'],
            'k_fit_range_cm': results['k_fit_range_cm'],
            'k_linked_guess_cm': results['k_linked_guess_cm'],
        },
        'hankel': {
            'lambda_p_nm': results['hankel_lambda_p_nm'],
            'q_p_1e5cm-1': results['hankel_q_p_1e5cm-1'],
            'damping': results['hankel_damping'],
        },
        'sqrtx': {
            'lambda_p_nm': results['1/sqrtx_lambda_p_nm'],
            'q_p_1e5cm-1': results['1/sqrtx_q_p_1e5cm-1'],
            'damping': results['1/sqrtx_damping'],
        },
        'fft': {
            'lambda_p_nm': results['fft_lambda_p_nm'],
            'q_p_1e5cm-1': results['fft_q_p_1e5cm-1'],
            'damping': results['fft_damping'],
        },
    }

# ---- 3. Save ----
with open(PICKLE_PATH, 'wb') as f:
    pickle.dump(all_results, f)

print(f'\nSaved {PICKLE_PATH} and {len(all_results)*3} figures to {SAVE_DIR}/')
print(f"\n{'wn':>8} | {'CHT nm':>8} | {'Hankel nm':>9} | {'1/sqrtx nm':>10}")
for wn, r in all_results.items():
    print(f"{wn:>8} | {r['cht']['lambda_p_nm']:>8.1f} | {r['hankel']['lambda_p_nm']:>9.1f} | {r['sqrtx']['lambda_p_nm']:>10.1f}")
