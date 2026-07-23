"""
Cross-check lp2's raw linecut O3A against a labmate's (Weicheng) independently-extracted
linecut of the same graphene_4x1_manual lp2 (interior) region -- same underlying scan,
different extraction pipeline/starting point (2026-07-23). See log_and_conventions.md
for the full writeup; this script was originally run ad hoc (not saved) and is being
saved as a permanent script now so it can be re-run/tweaked without rebuilding from
scratch.

Method: per wn, cross-correlate our O3A against Weicheng's (over a scan of candidate
shifts) to find the best alignment -- the two extractions' starting points genuinely
differ, confirmed to be a consistent ~370nm offset across all 13 shared wn. Both curves
z-scored before comparison since absolute amplitude scales/normalization differ between
the two extractions; only shape is comparable.
"""
import os
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams.update({
    'font.size': 11, 'font.family': 'Arial',
    'axes.linewidth': 1.3, 'xtick.direction': 'in', 'ytick.direction': 'in',
    'xtick.top': True, 'ytick.right': True
})

OUT_PATH = 'figures/graphene_4x1_manual/lp2/compare_with_weicheng_linecut.png'

wc = pd.read_csv("data/Weicheng's linecut/4x1um_lp2/MoO3_interior_S3amp_linecuts_waterfall.csv")
wc.columns = [c.replace('_cm?????', '') for c in wc.columns]
wns = ['860', '870', '880', '890', '900', '911', '920', '930', '941', '950', '960', '970', '980']

shifts = np.arange(-0.1, 0.6, 0.005)
fig, axes = plt.subplots(4, 4, figsize=(17, 12), sharex=True)

for ax, wn in zip(axes.flat, wns):
    ours = pd.read_csv(f'data/graphene_4x1_manual/{wn}cm-1_AVG_lp2.csv')
    x_ours, y_ours = ours['distance_nm'].values / 1000.0, ours['O3A'].values
    x_wc, y_wc = wc['pixel_index'].values, wc[wn].values

    x_common = np.arange(0.0, 1.4, 0.005)
    y_ours_i = np.interp(x_common, x_ours, y_ours, left=np.nan, right=np.nan)

    best_shift, best_corr = None, -2
    for s in shifts:
        y_wc_i = np.interp(x_common, x_wc - s, y_wc, left=np.nan, right=np.nan)
        mask = ~np.isnan(y_ours_i) & ~np.isnan(y_wc_i)
        if mask.sum() < 50:
            continue
        c = np.corrcoef(y_ours_i[mask], y_wc_i[mask])[0, 1]
        if c > best_corr:
            best_corr, best_shift = c, s

    y_wc_aligned = np.interp(x_common, x_wc - best_shift, y_wc, left=np.nan, right=np.nan)
    z_ours = (y_ours_i - np.nanmean(y_ours_i)) / np.nanstd(y_ours_i)
    z_wc = (y_wc_aligned - np.nanmean(y_wc_aligned)) / np.nanstd(y_wc_aligned)

    ax.plot(x_common * 1000, z_ours, color='#1c7293', lw=1.6, label='ours')
    ax.plot(x_common * 1000, z_wc, color='#d6604d', lw=1.6, alpha=0.85, label='Weicheng')
    ax.set_title(f'{wn}cm-1  (shift={best_shift*1000:.0f}nm, r={best_corr:.3f})', fontsize=10.5)
    ax.set_xlim(0, 1200)
    ax.set_xticks(np.arange(0, 1201, 200))
    ax.tick_params(labelsize=9, labelbottom=True)  # sharex hides these by default -- force on
    ax.set_xlabel('Distance (nm, our x=0 edge convention)', fontsize=9)

for ax in axes.flat[len(wns):]:
    ax.axis('off')
axes[0, 0].legend(fontsize=9, frameon=False)
for ax in axes[:, 0]:
    ax.set_ylabel('Amp (z-scored)', fontsize=9)
fig.suptitle("lp2: our O3A vs Weicheng's O3A (each wn independently shift-aligned, z-scored for shape comparison)",
             fontsize=13, y=1.0)
fig.tight_layout()
fig.savefig(OUT_PATH, dpi=170, bbox_inches='tight')
print(f'Saved {OUT_PATH}')
