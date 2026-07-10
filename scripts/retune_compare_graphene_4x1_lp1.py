"""
Regenerate fit_comparison + rp_dispersion figures for graphene_4x1_manual lp1
("edge") from the CURRENT tuned parameters in fitting_pipeline_graphene_4x1_lp1.ipynb,
into a separate versioned folder (does not touch/overwrite any existing committed
output) -- for side-by-side before/after comparison after a re-tuning pass.

Set VERSION_SUFFIX below before each new comparison round.
"""
import json, re, os, sys
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append('.')
sys.path.append('/Users/shizhe/envsetting')
import matplotlib
matplotlib.use('Agg')
import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar
import nanoftir_Shizhe as nanof
from loadnpz_combined import loadnpz
from snippet.nearfield import rp_schematic
from NearFieldOptics import Materials as M

VERSION_SUFFIX = 'v2'  # bump this each retuning round; never reuse a suffix already on disk

NOTEBOOK_PATH = 'fitting_pipeline_graphene_4x1_lp1.ipynb'
DATA_DIR = 'data/graphene_4x1_manual'
LP = 'lp1'
GEOM_CSV = f'{DATA_DIR}/linecut_geometry.csv'
MAP_NPZ_PATH = 'data/processed_4x1um/860-1000_AVG.npz'

PICKLE_PATH = f'data/fit_results_graphene_4x1_lp1_{VERSION_SUFFIX}.pkl'
FIT_COMPARISON_DIR = f'figures/graphene_4x1_manual/lp1_{VERSION_SUFFIX}/fit_comparison'
RP_DIR = f'figures/rp_dispersion_{VERSION_SUFFIX}'
os.makedirs(FIT_COMPARISON_DIR, exist_ok=True)
os.makedirs(RP_DIR, exist_ok=True)

# guard against accidentally clobbering a previous comparison round
for p in (PICKLE_PATH,):
    if os.path.exists(p):
        raise FileExistsError(f"{p} already exists -- bump VERSION_SUFFIX before rerunning")

align_dict = {f'{wn}cm-1': 0 for wn in
              (860, 870, 880, 890, 900, 911, 920, 930, 941, 950, 960, 970, 980, 991, 1000)}

# ---- 1. Parse the CURRENT (just re-tuned) per-wn parameters from the notebook ----
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

# ---- 2. Fit + build fit_comparison figures, collect pkl-shaped results ----
geom_df = pd.read_csv(GEOM_CSV)
geom_df = geom_df[geom_df['lp'] == LP]
scan = loadnpz(MAP_NPZ_PATH)

all_results = {}
cht_results_by_wn, hankel_out_by_wn, sqrtx_out_by_wn = {}, {}, {}

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
    plt.close(fig_cht); plt.close(fig_rs)

    cht_results_by_wn[wn] = cht_results
    hankel_out_by_wn[wn] = outs['hankel']
    sqrtx_out_by_wn[wn] = outs['1/sqrtx']

    all_results[wn] = {
        'cht': {'lambda_p_nm': cht_results['lambda_p_nm'], 'q_p_1e5cm-1': cht_results['q_p_1e5cm_1'],
                'damping': cht_results['damping'], 'k_fit_range_cm': cht_results['k_fit_range_cm'],
                'k_linked_guess_cm': cht_results['k_linked_guess_cm']},
    }
    for pf, key in (('hankel', 'hankel'), ('1/sqrtx', 'sqrtx')):
        o = outs[pf]
        pp, d, met = o['params'], o['derived'], o['metrics']
        damp_key = 'q_imag_um^-1' if pf == 'hankel' else 'alpha_env_um^-1'
        all_results[wn][key] = {
            'lambda_p_nm': pp['lambda_p_um'] * 1000, 'q_p_1e5cm-1': d['q_cm^-1'] / 1e5,
            'damping': d['q_rad_per_um'] / pp[damp_key] if pp[damp_key] > 1e-9 else np.inf,
        }

    map_panel = None
    row = geom_df[geom_df['wn'] == wn]
    if len(row):
        row = row.iloc[0]
        map_panel = dict(
            amp_data=scan['maps'][wn_val]['O3A'], pixelsize_um=scan['pixelsize_nm'] / 1000.0,
            center_um=(row['center_x_um'], row['center_y_um']), angle_deg=row['angle_deg'],
            radius_um=row['radius_um'], rect_width_um=row['rect_width_um'], title='O3A + linecut ROI',
        )

    fig1, curves = nanof.plot_fit_comparison_combined(
        x_f, sig_f, cht_results, outs['hankel'], outs['1/sqrtx'], wn,
        save_path=f'{FIT_COMPARISON_DIR}/{wn}_combined.png')
    fig2, _ = nanof.plot_fit_comparison_panels(
        x_f, sig_f, cht_results, outs['hankel'], outs['1/sqrtx'], wn,
        map_panel=map_panel, save_path=f'{FIT_COMPARISON_DIR}/{wn}_panels.png')
    plt.close(fig1); plt.close(fig2)
    print(f"  {wn}: CHT={curves['cht']['lambda_p_nm']:.1f}nm  Hankel={curves['hankel']['lambda_p_nm']:.1f}nm  "
          f"1/sqrtx={curves['1/sqrtx']['lambda_p_nm']:.1f}nm")

with open(PICKLE_PATH, 'wb') as f:
    pickle.dump(all_results, f)
print(f'\nSaved {PICKLE_PATH} and {len(all_results)*2} fit_comparison figures to {FIT_COMPARISON_DIR}/')

# ---- 3. rp_dispersion, same MoO3/stack/Ef-fit recipe as rp_dispersion_plot.ipynb ----
modes_xy = [(972, 4, 820, 4)]
modes_z = [(1004, 2, 958, 2)]
MoO3 = M.AnisotropicMaterial(eps_infinity=[4, 5.2, 2.4], phonon_params=[modes_xy, modes_xy, modes_z])
materials_4x1 = {'SiO2': M.SiO2_300nm, 'MoO3': MoO3}
SUBS = M.Si
stack_4x1 = [('MoO3', 4.2), 'graphene', ('SiO2', 300)]
GAMMA = 20
METHODS = ['cht', 'hankel', 'sqrtx', 'fft']  # fft not in all_results here (fit_comparison doesn't compute fft) -- skip it below

def build_pts(pkl_dict):
    out = {}
    for m in ('cht', 'hankel', 'sqrtx'):
        rows = [{'w(cm-1)': int(wn.replace('cm-1', '')), 'q(cm-1)': r[m]['q_p_1e5cm-1'] * 1e5}
                for wn, r in pkl_dict.items()]
        out[f'{m}_qp'] = pd.DataFrame(rows).sort_values('w(cm-1)').reset_index(drop=True)
    return out

pts = build_pts(all_results)

def build_tm_layers(stack, materials, Ef, gamma):
    Efcm = Ef * 8065.5
    graphene = M.SingleLayerGraphene(chemical_potential=Efcm, gamma=gamma)
    tm_layers = []
    for item in stack:
        if isinstance(item, str) and item.lower() == 'graphene':
            tm_layers.append(graphene)
        else:
            name, thickness_nm = item
            tm_layers.append((materials[name], thickness_nm * 1e-7))
    return tm_layers

def fit_Ef(stack, materials, subs, w_data, q_data, gamma=GAMMA, ef_bounds=(0.1, 1.5)):
    w_data = np.asarray(w_data, dtype=float)
    q_data = np.asarray(q_data, dtype=float)
    def neg_score(Ef):
        layersTM = M.LayeredMediaTM(*build_tm_layers(stack, materials, Ef, gamma), exit=subs)
        rp_grid = np.asarray(layersTM.reflection_p(w_data, q=q_data))
        return -float(np.sum(np.imag(np.diag(rp_grid))))
    res = minimize_scalar(neg_score, bounds=ef_bounds, method='bounded')
    return float(res.x)

ef_by_method = {}
for m in ('cht', 'hankel', 'sqrtx'):
    df = pts[f'{m}_qp']
    ef_by_method[m] = fit_Ef(stack_4x1, materials_4x1, SUBS, df['w(cm-1)'].values, df['q(cm-1)'].values)
print('E_F (eV), lp1 edge, retuned:', {k: round(v, 3) for k, v in ef_by_method.items()})

FREQS = np.linspace(850, 1010, 300)
QS = np.linspace(1e4, 5.5e5, 400)
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
for ax, method in zip(axes, ('cht', 'hankel', 'sqrtx')):
    rp_schematic(stack_4x1, materials_4x1, subs=SUBS, freqs=FREQS, qs=QS,
                 Ef=ef_by_method[method], gamma=GAMMA, figAx=(fig, ax), vmin=0, vmax=8,
                 pts=pts, pts_chn=method, show_Ef=True, ytickspacing=20)
    ax.set_title(f'MoO3 (4.2nm)/G/SiO2/Si -- edge ({VERSION_SUFFIX})\nmethod={method}, E_F={ef_by_method[method]:.2f}eV (fit), gamma={GAMMA}cm-1', fontsize=11)
fig.tight_layout()
rp_path = f'{RP_DIR}/graphene_4x1_lp1_edge_{VERSION_SUFFIX}_rp_dispersion.png'
fig.savefig(rp_path, dpi=200, bbox_inches='tight')
plt.close(fig)
print(f'Saved {rp_path}')
