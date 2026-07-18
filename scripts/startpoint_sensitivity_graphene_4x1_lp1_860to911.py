"""
Starting-point sensitivity check for graphene_4x1_manual lp1, wn=860-911cm-1 (the same 6 wn
as peak_spacing_crosscheck_graphene_4x1_lp1.py). Generalizes the 860cm-1-only prototype
(startpoint_sensitivity_graphene_4x1_lp1_860cm-1.py, now superseded/removed) once its result
was reviewed and approved.

Motivation: the real-space Hankel/1-sqrtx fit window's *lower* bound (xr[0], "how close to
the edge does the fit start") is a hand-tuned parameter that directly sets q_p, but there was
no visual/quantitative check of how sensitive q_p actually is to this choice. For 860cm-1,
scanning xr[0] from 0 to 0.50um in 0.02um steps (xr[1] fixed at 1.5um) revealed a smooth,
non-monotonic oscillation (period ~160-200nm, roughly a fraction of lambda_p at this wn) with
two full periods over 0-400nm, then hankel/1-sqrtx start to diverge from each other past
~450nm even though RMSE is lowest there (see log entry -- likely too few oscillation periods
left in the fit window for a robust fit, not "genuinely better"). This script re-runs that
same scan for the other 5 wn (870/880/890/900/911) to see whether the oscillation
period/amplitude/divergence-onset scale consistently with wn (hence with lambda_p).

Uses fit_cavity_prefactor_compare directly (not the compare_cavity_models wrapper) since we
only need q_p per start point, not a full comparison figure at every scan point. Does NOT
touch fitting_pipeline_graphene_4x1_lp1.ipynb's own tuned xr_range -- side investigation only.
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

SAVE_DIR = 'figures/graphene_4x1_manual/lp1/startingpoint_sensitivity'
os.makedirs(SAVE_DIR, exist_ok=True)

XR1 = 1.5   # right bound, fixed -- matches the currently-tuned value for all 6 wn
YC_UM = 1.9
WIN, PROM = 3, 0.01

# fit_cavity_prefactor_compare's default lam_bound_factor=(0.7,1.2) clamps lambda to
# +-20-30% of the single fixed lam0_guess -- fine for production fits (lam0_guess is
# tuned per wn's already-chosen xr[0]), but too tight here: scanning xr[0] far from its
# tuned value can genuinely want a lambda outside that narrow window, and the optimizer
# clamps to the bound edge instead of converging to a real local minimum (silent
# artifact -- confirmed via exact-duplicate q values across neighboring start points:
# 911cm-1 had 11/26 points clamped, 900cm-1 had 2/26, 860-890 had 0/26 with the default).
# Widened here (scan only, NOT touched in fitting_pipeline_graphene_4x1_lp1.ipynb) so the
# scan can explore its actual sensitivity rather than hitting an artificial wall.
LAM_BOUND_FACTOR = (0.4, 2.0)

# lam0_guess_um per wn, read live from fitting_pipeline_graphene_4x1_lp1.ipynb's own
# per-wn real-space cells on 2026-07-17 (win_size/prominence/yc_um confirmed identical
# -- 3/0.01/1.9 -- across all 6, only lam0_guess_um differs per wn).
LAM0_GUESS = {
    860: 0.67, 870: 0.58, 880: 0.60, 890: 0.50, 900: 0.50, 911: 0.40,
}

start_points = np.round(np.arange(0.0, 0.50 + 1e-9, 0.02), 2)

# Region flagged by inspection as where all 6 wn's normalized q_p trends climb together
# fairly consistently (past the shared 40-80/200-240nm dips, before the >350nm divergence
# between hankel/sqrtx) -- highlighted on both overlay plots below for reference only,
# not used anywhere in the fits themselves.
GOOD_REGION_NM = (260, 340)

dfs_by_wn = {}
wn_list = list(LAM0_GUESS.keys())

for wn in wn_list:
    wn_str = f'{wn}cm-1'
    loaded = nanof.load_aligned_wn_signal(wn_str, 0, data_dir='data/graphene_4x1_manual', L_cutoff_fit=1.5, lp='lp1')
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

    # Per-wn figure: q_p vs start point (top) + RMSE vs start point (bottom, same x-axis)
    # -- same "data panel + quality panel" layout as compare_cavity_models's (a)/(b), so a
    # q_p dip can be checked directly against whether RMSE spikes at the same start point
    # (real bad fit) or not (could just be method disagreement/noise). No GOOD_REGION_NM
    # shading here -- that's a cross-wn observation, only meaningful on the combined
    # overlay plots below, not on a single wn's own curve.
    fig, (ax, axr) = plt.subplots(2, 1, figsize=(7, 6.5), sharex=True,
                                    gridspec_kw={'height_ratios': [2, 1], 'hspace': 0.06})
    colors = {'hankel': '#1c7293', '1/sqrtx': '#e08214'}
    labels = {'hankel': 'Hankel', '1/sqrtx': r'1/$\sqrt{x}$'}
    for pf in ('hankel', '1/sqrtx'):
        col = f'{pf}_q_1e5cm-1'
        ax.plot(df['start_um'] * 1000, df[col], 'o-', color=colors[pf], ms=6, lw=1.6, label=labels[pf])
        axr.plot(df['start_um'] * 1000, df[f'{pf}_rmse'], 'o-', color=colors[pf], ms=5, lw=1.4)
    ax.set_ylabel(r'$q_p$ ($10^5$ cm$^{-1}$)', fontweight='bold')
    ax.set_title(f'{wn_str}: fit q_p vs starting point (right bound fixed at {XR1}um)', fontsize=11.5)
    ax.legend(frameon=False, fontsize=10, loc='best')
    ax.tick_params(direction='in', top=True, right=True)
    axr.set_xlabel('Fit starting point, distance from edge (nm)', fontweight='bold')
    axr.set_ylabel('RMSE (a.u.)', fontweight='bold')
    axr.tick_params(direction='in', top=True, right=True)
    fig.tight_layout()
    out_path = f'{SAVE_DIR}/{wn_str}_q_vs_start_point.png'
    fig.savefig(out_path, dpi=200, bbox_inches='tight')
    plt.close(fig)

# Combined 6-wn overlay, one per prefactor -- ACTUAL q_p (not normalized), so the 6 wn sit
# on their own absolute q_p bands rather than collapsed onto a shared relative scale.
# GOOD_REGION_NM shaded for reference (no legend label -- purely a visual marker, not a
# claimed quantitative result). RMSE subpanel below, same per-wn color coding.
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
    ax_all.set_title(f'860-911cm-1: fit q_p vs starting point ({labels[pf]})', fontsize=11.5)
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

# Real-space data preview, all 6 wn, GOOD_REGION_NM shaded directly on the actual signal --
# so the abstract "good region" from the q_p scan can be read against which part of each
# wn's own oscillation it actually falls on.
fig_rs, axes_rs = plt.subplots(2, 3, figsize=(15, 8), sharex=True)
for ax, wn in zip(axes_rs.flat, wn_list):
    wn_str = f'{wn}cm-1'
    loaded = nanof.load_aligned_wn_signal(wn_str, 0, data_dir='data/graphene_4x1_manual', L_cutoff_fit=1.5, lp='lp1')
    x_f, sig_f = loaded['x_f'], loaded['sig_f']
    ax.plot(x_f * 1000, sig_f, '-', color='#555555', lw=1.2, alpha=0.8)
    ax.axhline(0, color='k', ls='--', lw=0.6, alpha=0.4)
    ax.axvspan(*GOOD_REGION_NM, color='#4daf4a', alpha=0.15, zorder=0)
    ax.set_title(wn_str, fontsize=12, fontweight='bold')
    ax.set_xlim(0, 1000)
    ax.tick_params(direction='in', top=True, right=True)
for ax in axes_rs[-1]:
    ax.set_xlabel('Distance from edge (nm)', fontweight='bold')
for ax in axes_rs[:, 0]:
    ax.set_ylabel('Amp (bg-sub, a.u.)', fontweight='bold')
fig_rs.suptitle(f'Real-space signal, {GOOD_REGION_NM[0]}-{GOOD_REGION_NM[1]}nm region shaded', fontsize=13)
fig_rs.tight_layout()
out_path_rs = f'{SAVE_DIR}/860to911cm-1_realspace_with_good_region.png'
fig_rs.savefig(out_path_rs, dpi=200, bbox_inches='tight')
plt.close(fig_rs)
print(f'Saved {out_path_rs}')
