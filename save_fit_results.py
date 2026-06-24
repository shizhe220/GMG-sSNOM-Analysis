"""
One-off run: regenerate the cht/realspace/fft figures for all 15 wavenumbers
using the parameters currently tuned in fitting_pipeline.ipynb, and assemble a
results dict (lambda_p, momentum q, damping per method) plus a CHT
old-window-vs-new-window deviation table vs the real-space Hankel benchmark.
Saved to data/fit_results.pkl (pickle) for later inspection / CSV export.
"""
import json
import re
import pickle
import sys
sys.path.append('/Users/shizhe/envsetting')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import nanoftir_Shizhe as nanof

align_dict = {
    '860cm-1': 60, '870cm-1': 75, '880cm-1': 75, '890cm-1': 75, '900cm-1': 75,
    '911cm-1': 65, '920cm-1': 65, '930cm-1': 75, '941cm-1': 63, '950cm-1': 75,
    '960cm-1': 65, '970cm-1': 65, '980cm-1': 45, '991cm-1': 61, '1000cm-1': 55,
}

# CHT parameters as originally auto-tuned (before this session's manual k_fit_range
# adjustments) -- used only as the "old window" reference point in the deviation table.
old_cht_params = {
    '860cm-1': dict(k_fit_range_cm=(0.5, 5.0), k_linked_guess_cm=1.5),
    '870cm-1': dict(k_fit_range_cm=(0.5, 6.0), k_linked_guess_cm=1.60),
    '880cm-1': dict(k_fit_range_cm=(0.5, 6.0), k_linked_guess_cm=1.75),
    '890cm-1': dict(k_fit_range_cm=(0.5, 6.0), k_linked_guess_cm=1.90),
    '900cm-1': dict(k_fit_range_cm=(0.5, 6.0), k_linked_guess_cm=2.10),
    '911cm-1': dict(k_fit_range_cm=(0.5, 6.0), k_linked_guess_cm=2.20),
    '920cm-1': dict(k_fit_range_cm=(0.5, 6.0), k_linked_guess_cm=2.30),
    '930cm-1': dict(k_fit_range_cm=(0.5, 6.0), k_linked_guess_cm=2.35),
    '941cm-1': dict(k_fit_range_cm=(0.5, 6.0), k_linked_guess_cm=2.60),
    '950cm-1': dict(k_fit_range_cm=(0.5, 6.0), k_linked_guess_cm=2.70),
    '960cm-1': dict(k_fit_range_cm=(0.5, 6.0), k_linked_guess_cm=3.00),
    '970cm-1': dict(k_fit_range_cm=(0.5, 6.0), k_linked_guess_cm=3.20),
    '980cm-1': dict(k_fit_range_cm=(0.5, 6.0), k_linked_guess_cm=3.50),
    '991cm-1': dict(k_fit_range_cm=(0.5, 6.0), k_linked_guess_cm=3.80),
    '1000cm-1': dict(k_fit_range_cm=(0.5, 6.0), k_linked_guess_cm=4.00),
}
old_x_start_cht = 0.22
old_L_cutoff_cht = 1.2

# ---- 1. Parse the current (user-tuned) per-wn parameters straight out of the notebook ----
with open('fitting_pipeline.ipynb') as f:
    nb = json.load(f)

cht_idx = {7: 860, 15: 870, 23: 880, 31: 890, 39: 900, 47: 911, 55: 920, 63: 930,
           71: 941, 79: 950, 87: 960, 95: 970, 103: 980, 111: 991, 119: 1000}
rs_idx = {9: 860, 17: 870, 25: 880, 33: 890, 41: 900, 49: 911, 57: 920, 65: 930,
          73: 941, 81: 950, 89: 960, 97: 970, 105: 980, 113: 991, 121: 1000}
fft_idx = {11: 860, 19: 870, 27: 880, 35: 890, 43: 900, 51: 911, 59: 920, 67: 930,
           75: 941, 83: 950, 91: 960, 99: 970, 107: 980, 115: 991, 123: 1000}

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

# ---- 2. Run the NEW (current notebook) parameters, saving all 3 figures per wn ----
all_results = {}
for wn, p in wn_params.items():
    print(f'Running {wn} (new/current params) ...')
    new_results, _ = nanof.run_wn_comparison(
        wn, align_shift_nm=align_dict[wn], k_linked_guess_cm=p['k_linked_guess_cm'],
        k_fit_range_cm=p['k_fit_range_cm'], x_start_cht=p['x_start_cht'], L_cutoff_cht=p['L_cutoff_cht'],
        lam0_guess_um=p['lam0_guess_um'], xr_range_rs=p['xr_range_rs'],
        fft_xr=p['fft_xr'], fft_q_guess=p['fft_q_guess'],
        save_dir='figures', show=False)

    print(f'Running {wn} (old baseline CHT window, no save) ...')
    old_p = old_cht_params[wn]
    old_cht_results, old_fig = nanof.fit_and_plot_cht(
        nanof.load_aligned_wn_signal(wn, align_dict[wn])['x_f'],
        nanof.load_aligned_wn_signal(wn, align_dict[wn])['sig_f'],
        wn, x_start_cht=old_x_start_cht, L_cutoff_cht=old_L_cutoff_cht,
        k_fit_range_cm=old_p['k_fit_range_cm'], k_linked_guess_cm=old_p['k_linked_guess_cm'])
    plt.close(old_fig)

    hankel_lambda = new_results['hankel_lambda_p_nm']
    old_dev_pct = 100 * abs(old_cht_results['lambda_p_nm'] - hankel_lambda) / hankel_lambda
    new_dev_pct = 100 * abs(new_results['lambda_p_nm'] - hankel_lambda) / hankel_lambda

    all_results[wn] = {
        'cht': {
            'lambda_p_nm': new_results['lambda_p_nm'],
            'q_p_1e5cm-1': new_results['q_re'] / 10.0,
            'damping': new_results['damping'],
            'k_fit_range_cm': new_results['k_fit_range_cm'],
            'k_linked_guess_cm': new_results['k_linked_guess_cm'],
        },
        'hankel': {
            'lambda_p_nm': new_results['hankel_lambda_p_nm'],
            'q_p_1e5cm-1': new_results['hankel_q_p_1e5cm-1'],
            'damping': new_results['hankel_damping'],
        },
        'sqrtx': {
            'lambda_p_nm': new_results['1/sqrtx_lambda_p_nm'],
            'q_p_1e5cm-1': new_results['1/sqrtx_q_p_1e5cm-1'],
            'damping': new_results['1/sqrtx_damping'],
        },
        'fft': {
            'lambda_p_nm': new_results['fft_lambda_p_nm'],
            'q_p_1e5cm-1': new_results['fft_q_p_1e5cm-1'],
            'damping': new_results['fft_damping'],
        },
        'cht_old_window': {
            'lambda_p_nm': old_cht_results['lambda_p_nm'],
            'damping': old_cht_results['damping'],
            'k_fit_range_cm': old_p['k_fit_range_cm'],
            'k_linked_guess_cm': old_p['k_linked_guess_cm'],
        },
        'cht_deviation_vs_hankel_pct': {
            'old_window': round(old_dev_pct, 1),
            'new_window': round(new_dev_pct, 1),
        },
    }

# ---- 3. Save ----
with open('data/fit_results.pkl', 'wb') as f:
    pickle.dump(all_results, f)

print('\nSaved data/fit_results.pkl')
print('\nSummary (CHT lambda deviation from Hankel benchmark, old window vs new window):')
print(f"{'wn':>8} | {'old_dev%':>9} | {'new_dev%':>9} | new k_fit_range_cm")
for wn, r in all_results.items():
    dev = r['cht_deviation_vs_hankel_pct']
    print(f"{wn:>8} | {dev['old_window']:>9} | {dev['new_window']:>9} | {r['cht']['k_fit_range_cm']}")
