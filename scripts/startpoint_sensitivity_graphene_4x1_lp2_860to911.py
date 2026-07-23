"""
Starting-point sensitivity check for graphene_4x1_manual lp2, wn=860-911cm-1. Sibling of
startpoint_sensitivity_graphene_4x1_lp1_860to911.py (2026-07-17), same methodology, run on
lp2's data instead. lam0_guess_um values below are read live from
fitting_pipeline_graphene_4x1_lp2.ipynb's own per-wn real-space cells on 2026-07-23 --
lp2's data was re-extracted that same day with the perpendicular averaging width doubled
(0.077->0.154um) to reduce noise, and these are the user's freshly-retuned lam0_guess
values against the new linecuts (still a rough/first pass, not final).

GOOD_REGION_NM is filled in AFTER inspecting this script's own output (unlike lp1's script,
where it was already known) -- see log_and_conventions.md for the value chosen and why.
"""
import sys, os
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append('.')
sys.path.append('/Users/shizhe/envsetting')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import nanoftir_Shizhe as nanof

plt.rcParams.update({
    'font.size': 12, 'font.family': 'Arial',
    'axes.linewidth': 1.5, 'xtick.major.width': 1.5, 'ytick.major.width': 1.5,
    'xtick.direction': 'in', 'ytick.direction': 'in', 'xtick.top': True,
    'ytick.right': True
})

SAVE_DIR = 'figures/graphene_4x1_manual/lp2/startingpoint_sensitivity'
os.makedirs(SAVE_DIR, exist_ok=True)

XR1 = 1.5
YC_UM = 1.9
WIN, PROM = 3, 0.01

# same widened bound as the lp1 scan -- lets the scan explore its actual sensitivity
# instead of silently clamping to +-20-30% of lam0_guess (see lp1 script for the full
# rationale/verification).
LAM_BOUND_FACTOR = (0.4, 2.0)

LAM0_GUESS = {
    860: 0.647, 870: 0.64, 880: 0.56, 890: 0.52, 900: 0.49, 911: 0.6,
}

start_points = np.round(np.arange(0.0, 0.50 + 1e-9, 0.02), 2)

# Chosen after inspecting the first pass of this script's own output (2026-07-23): unlike
# lp1, hankel is NOT stable for 900/911cm-1 anywhere in 0-500nm here (900cm-1 never
# plateaus, climbs the whole way; 911cm-1 hops between local minima, e.g. crashing to
# ~0.52e5 right at lp1's old 260-340nm window). 1/sqrtx is comparably much better-behaved
# for lp2 across all 6 wn from roughly 260-380nm (860/870/890 near-flat, 900 finally
# plateaus ~1.13-1.33, 911 flat except one single-point dip right at 420nm, just past this
# window). See log entry for the full hankel-vs-sqrtx comparison and why hankel's own
# instability here is treated as more evidence for the two-wave model at 900cm-1+, not a
# start-point problem to fix by moving this window.
# Narrowed to 200-300nm on user's request (2026-07-23): lp2's linecut is only ~1000nm
# total (vs lp1's longer cuts), so a window extending to 380nm eats too much of the
# already-short usable range -- 200-300nm still sits past the initial trough for all 6 wn
# (see the realspace preview) while leaving more of the cut available for the fit itself.
GOOD_REGION_NM = (200, 300)

dfs_by_wn = {}
wn_list = list(LAM0_GUESS.keys())

for wn in wn_list:
    wn_str = f'{wn}cm-1'
    loaded = nanof.load_aligned_wn_signal(wn_str, 0, data_dir='data/graphene_4x1_manual', L_cutoff_fit=1.5, lp='lp2')
    x_f, sig_f = loaded['x_f'], loaded['sig_f']

    rows = []
    for xr0 in start_points:
        row = {'start_um': xr0}
        for pf in ('hankel', '1/sqrtx'):
            out = nanof.fit_cavity_prefactor_compare(
                x_f, sig_f, xr=(xr0, XR1), yc_um=YC_UM, fit_yc=False, edges='single',
                prefactor=pf, win=WIN, prom=PROM, lam0_guess=LAM0_GUESS[wn],
                lam_bound_factor=LAM_BOUND_FACTOR)
            row[f'{pf}_q_1e5cm-1'] = out['derived']['q_cm^-1'] / 1e5
            row[f'{pf}_rmse'] = out['metrics']['rmse']
        rows.append(row)

    df = pd.DataFrame(rows).sort_values('start_um').reset_index(drop=True)
    df.to_csv(f'{SAVE_DIR}/{wn_str}_q_vs_start_point.csv', index=False)
    print(f'{wn_str}: saved {SAVE_DIR}/{wn_str}_q_vs_start_point.csv  '
          f'(hankel q range {df["hankel_q_1e5cm-1"].min():.3f}-{df["hankel_q_1e5cm-1"].max():.3f})')
    dfs_by_wn[wn] = df

    fig, (ax, axr) = plt.subplots(2, 1, figsize=(7, 6.5), sharex=True,
                                    gridspec_kw={'height_ratios': [2, 1], 'hspace': 0.06})
    colors = {'hankel': '#1c7293', '1/sqrtx': '#e08214'}
    labels = {'hankel': 'Hankel', '1/sqrtx': r'1/$\sqrt{x}$'}
    for pf in ('hankel', '1/sqrtx'):
        col = f'{pf}_q_1e5cm-1'
        ax.plot(df['start_um'] * 1000, df[col], 'o-', color=colors[pf], ms=6, lw=1.6, label=labels[pf])
        axr.plot(df['start_um'] * 1000, df[f'{pf}_rmse'], 'o-', color=colors[pf], ms=5, lw=1.4)
    ax.set_ylabel(r'$q_p$ ($10^5$ cm$^{-1}$)', fontweight='bold')
    ax.set_title(f'{wn_str} (lp2): fit q_p vs starting point (right bound fixed at {XR1}um)', fontsize=11.5)
    ax.legend(frameon=False, fontsize=10, loc='best')
    ax.tick_params(direction='in', top=True, right=True)
    axr.set_xlabel('Fit starting point, distance from edge (nm)', fontweight='bold')
    axr.set_ylabel('RMSE (a.u.)', fontweight='bold')
    axr.tick_params(direction='in', top=True, right=True)
    fig.tight_layout()
    out_path = f'{SAVE_DIR}/{wn_str}_q_vs_start_point.png'
    fig.savefig(out_path, dpi=200, bbox_inches='tight')
    plt.close(fig)

cmap = plt.get_cmap('viridis')
for pf, pf_tag in (('hankel', 'hankel'), ('1/sqrtx', 'sqrtx')):
    fig_all, (ax_all, axr_all) = plt.subplots(
        2, 1, figsize=(8, 7.5), sharex=True, gridspec_kw={'height_ratios': [2, 1], 'hspace': 0.06})
    col = f'{pf}_q_1e5cm-1'
    for i, wn in enumerate(wn_list):
        df = dfs_by_wn[wn]
        c = cmap(i / max(len(wn_list) - 1, 1))
        ax_all.plot(df['start_um'] * 1000, df[col], 'o-', ms=4, lw=1.4, color=c, label=f'{wn}cm-1')
        axr_all.plot(df['start_um'] * 1000, df[f'{pf}_rmse'], 'o-', ms=3.5, lw=1.2, color=c)

    ax_all.axvspan(*GOOD_REGION_NM, color='#4daf4a', alpha=0.12, zorder=0)
    axr_all.axvspan(*GOOD_REGION_NM, color='#4daf4a', alpha=0.12, zorder=0)
    ax_all.set_ylabel(r'$q_p$ ($10^5$ cm$^{-1}$)', fontweight='bold')
    ax_all.set_title(f'860-911cm-1 (lp2): fit q_p vs starting point ({labels[pf]})', fontsize=11.5)
    ax_all.legend(frameon=False, fontsize=9, loc='best', ncol=2)
    ax_all.tick_params(direction='in', top=True, right=True)
    axr_all.set_xlabel('Fit starting point, distance from edge (nm)', fontweight='bold')
    axr_all.set_ylabel('RMSE (a.u.)', fontweight='bold')
    axr_all.tick_params(direction='in', top=True, right=True)
    fig_all.tight_layout()
    out_path_all = f'{SAVE_DIR}/860to911cm-1_q_vs_start_point_overlay_{pf_tag}.png'
    fig_all.savefig(out_path_all, dpi=200, bbox_inches='tight')
    plt.close(fig_all)
    print(f'Saved {out_path_all}')

fig_rs, axes_rs = plt.subplots(2, 3, figsize=(15, 8), sharex=True)
for ax, wn in zip(axes_rs.flat, wn_list):
    wn_str = f'{wn}cm-1'
    loaded = nanof.load_aligned_wn_signal(wn_str, 0, data_dir='data/graphene_4x1_manual', L_cutoff_fit=1.5, lp='lp2')
    x_f, sig_f = loaded['x_f'], loaded['sig_f']
    ax.plot(x_f * 1000, sig_f, '-', color='#555555', lw=1.2, alpha=0.8)
    ax.axhline(0, color='k', ls='--', lw=0.6, alpha=0.4)
    ax.axvspan(*GOOD_REGION_NM, color='#4daf4a', alpha=0.15, zorder=0)
    ax.set_title(f'{wn_str} (lp2)', fontsize=12, fontweight='bold')
    ax.set_xlim(0, 1000)
    ax.tick_params(direction='in', top=True, right=True)
for ax in axes_rs[-1]:
    ax.set_xlabel('Distance from edge (nm)', fontweight='bold')
for ax in axes_rs[:, 0]:
    ax.set_ylabel('Amp (bg-sub, a.u.)', fontweight='bold')
fig_rs.suptitle(f'lp2 real-space signal, {GOOD_REGION_NM[0]}-{GOOD_REGION_NM[1]}nm region shaded', fontsize=13)
fig_rs.tight_layout()
out_path_rs = f'{SAVE_DIR}/860to911cm-1_realspace_with_good_region.png'
fig_rs.savefig(out_path_rs, dpi=200, bbox_inches='tight')
plt.close(fig_rs)
print(f'Saved {out_path_rs}')
