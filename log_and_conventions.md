# Log & Conventions

## Folder Structure (as of 2026-06-29)

```
GMG/
├── fitting_pipeline.ipynb        CHT/Hankel/1-sqrtx/FFT fitting, per wn (15 cells each)
├── manual_linecut_pipeline.ipynb Manual per-wn linecut extraction from processed_3x1um/*.npz
├── nanoftir_Shizhe.py             Core fitting module (CHT, real-space, FFT) -- imported by
│                                  both notebooks + webtool/. Kept at root: moving it would
│                                  require editing every notebook's import cell for no benefit.
├── loadnpz.py                     Reader for the FwdBwd .npz summary format (2D channels +
│                                  ROI/line-profile metadata). Also kept at root for the same reason.
├── log_and_conventions.md         This file.
│
├── data/
│   ├── processed_3x1um/           Raw .npz per wn (full 2D O2/O3/O4 amp+phase + Z maps)
│   ├── graphene_3x1/               Old CSV-based line profiles (pre-npz pipeline)
│   ├── graphene_3x1_manual/        CSV export from manual_linecut_pipeline.ipynb (same columns
│   │                                as graphene_3x1/, drop-in replacement -- point fitting_pipeline's
│   │                                data_dir here to fit the manually-extracted+aligned linecuts)
│   ├── fit_results.pkl             Per-wn lambda/q/damping for CHT/Hankel/1-sqrtx/FFT
│   └── cht_vs_realspace_wavelength_comparison.csv
│
├── figures/
│   ├── cht/ realspace/ fft/        Per-wn fit figures (15 each)
│   ├── cht_diagnostics/            |T(k)| peak-location diagnostics (960-1000cm-1 k_fit_range issue)
│   ├── manual_linecut/             Z-plane-correction checks + waterfall from manual_linecut_pipeline
│   ├── overview/                   Amp/Phase/Z overview + waterfall snapshots
│   └── q_vs_wn_*.png               Combined + per-method momentum-vs-wavenumber summary plots
│
├── scripts/                        Standalone regeneration scripts (not imported by either notebook).
│   │                                Each chdir's to the repo root on startup (computed from __file__),
│   │                                so they can be run from anywhere, not just when cwd is the repo root.
│   ├── save_fit_results.py         Re-parses fitting_pipeline.ipynb's current per-wn CHT/RS/FFT params,
│   │                                regenerates all 45 figures + fit_results.pkl. Locates each wn's
│   │                                cells dynamically (via the "target_wn = '...'" anchor cell), not by
│   │                                hardcoded index -- those go stale if cells get inserted upstream.
│   ├── export_comparison_csv.py    fit_results.pkl -> cht_vs_realspace_wavelength_comparison.csv
│   └── make_overview_pptx.py       Builds slides/GMG_CHT_overview.pptx from the figures/ + comparison CSV
│
├── docs/                           Reference papers (Woessner et al. Nat. Mater. 2015 + SI)
├── slides/                         Presentation decks (.pptx)
├── webtool/                        Streamlit v1 interactive tool (click-to-place-linecut/align/q_guess).
│                                    Has known rough edges -- not yet the primary workflow.
└── dev_scripts_archive/            Pre-Claude-Code scratch scripts (gitignored, not part of the repo)
```

**Shared code that lives outside this repo**: `/Users/shizhe/envsetting/snippet/` -- in particular
`extract_linecut.py` (radial_rect_profile/radial_sector_profile, the 2D-map-to-1D-linecut math) and
`linecut_extraction.py` (Z-plane correction, extract_and_plot, plot_waterfall_3channel,
plot_mapping_waterfall). Project-agnostic on purpose (works for G-CIPS-style data too), so it's a
separate git checkout rather than living inside GMG/.

## Progress Log
*   **[2026-06-15 16:48]** Analyzed *Woessner et al. Nature Materials 2015* paper and SI.
*   **[2026-06-15 17:05]** Reviewed user's `nanoftir_Shizhe.py` to compare existing real-space fitting and FFT methods with the paper's Complex Hankel Transform (CHT) method.
*   **[2026-06-15 17:21]** Created implementation plan to add CHT to `nanoftir_Shizhe.py` and build an interactive `.ipynb` pipeline.
*   **[2026-06-15 17:34]** Updated CSV headers in `/data` to explicitly state `# Pixel size: 15 nm`.
*   **[2026-06-15 17:44]** Updated data processing to auto-detect boundary (substrate cutoff) at ~0.33 um.
*   **[2026-06-15 17:50]** Implemented paper-accurate APS plotting style with damping envelopes and created waterfall plot.
*   **[2026-06-16 02:08]** Compared CHT method with old `hankel` and `1/sqrtx` methods for $860\text{ cm}^{-1}$, creating a structured CSV `data/fitting_results_comparison.csv`.
*   **[2026-06-16 02:22]** User noted peak-to-peak distance is ~200nm, implying $\lambda_p = 400\text{ nm}$. Refined boundary to $x = 0.35\text{ \mu m}$ and adjusted initial guess $lam0\_guess$ to 400 nm.
*   **[2026-06-16 09:47]** Reconstructed `.ipynb` pipeline using `compare_cavity_models`. Simulated single edge by setting $xr[1] = 100.0\text{ \mu m}$ and $fit\_yc=False$. Adjusted CHT plot annotation positions.
*   **[2026-06-16 09:56]** Added native single-edge mode to `nanoftir_Shizhe.py`. Added `edges='single'` argument to `compare_cavity_models` and `fit_cavity_prefactor_compare` to eliminate infinity-boundary artifacts. Regenerated `.ipynb` to reflect this new elegant parameter.
*   **[2026-06-16 10:07]** Fixed `compare_cavity_models` `TypeError` (removed invalid `robust` arg). Fixed critical logical bug where `xr=(0.35, 1.2)` was shifting an already zero-referenced distance array; updated to `xr=(0.0, 1.2)`. Added explicit parameter blocks in `.ipynb` for both CHT and Real Space fits, and added in-line matplotlib plotting for CHT. Enhanced `compare_cavity_models` to output $\lambda$, $q$, and $\gamma^{-1}$ alongside AIC in the plot text box.
*   **[2026-06-16 10:35]** Explained CHT mathematical derivation (Complex Hankel Transform integral and Apodization Window). Upgraded the CHT plot in the Notebook to a 3-panel figure: Panel (c) explicitly visualizes the momentum space ($k$-space) transform $|T(k)|$, showing the resonance peak and the shaded fitting range `k_fit_range`.
*   **[2026-06-16 10:45]** Replaced `plt.show(fig)` with just `fig` in the Real Space plotting cell to fix rendering bugs in Jupyter widget backend. Modified `k_fit_range` for CHT to `(20, 35)` to avoid the massive low-$k$ background from incomplete background subtraction pulling the fit. Added Cell 4 to `fitting_pipeline.ipynb` that invokes `plot_channel_fft` to allow direct comparison between the CHT transformation and Standard 1D FFT.
*   **[2026-06-16 10:55]** Resolved Jupyter module caching bug (`TypeError: unexpected keyword argument 'k_plot_range'`) by injecting `importlib.reload(nanoftir)` into Cell 1. Exposed `q_guess=[5, 30]` in the FFT comparison cell for user tuning.
*   **[2026-06-17 10:25]** Implemented a complete batch data alignment pipeline for the `graphene_3x1` dataset across 15 wavenumbers. Replaced single-file boundary truncation with a multi-file interactive alignment cell using an `align_dict`.
*   **[2026-06-17 10:38]** Fixed a critical bug where the new pipeline fed raw `O3A` amplitude to the CHT algorithm by restoring the `savgol_filter` background subtraction step before fitting.
*   **[2026-06-17 10:49]** Upgraded waterfall visualization aesthetics: enforced `Z_nm_custom = Z_nm_corrected - min`, implemented separate `normalize_A` and `normalize_P` flags, utilized `coolwarm` colormap, removed Y-ticks, and exposed vertical spacing variables `offset_A`, `offset_P`, and `offset_Z`.
*   **[2026-06-17 13:26]** Reversed the sorting of wavenumbers in the waterfall plot generation (`reverse=True`) to adhere to the convention of plotting higher wavenumbers at the top and lower wavenumbers at the bottom.

## Conventions & Rules

### 1. Extensive Logging Policy
*   **Requirement**: Every modification to code, parameters, or logic must be recorded in the Progress Log with an exact timestamp `[YYYY-MM-DD HH:MM]`.
*   **Detail Level**: Logs must specify *what* parameter was changed (e.g., $xr$, $lam0\_guess$, boundary index), *why* it was changed, and *where* it was modified.

### 1. Code modification policy
*   **Non-destructive**: All new core mathematical functions (e.g., `complex_hankel_transform`) will be appended to existing modules (`nanoftir_Shizhe.py`) without altering the signature or behavior of existing functions (like `joint_fit` or `plot_nf_fft`).

### 2. Plotting Style (APS Standard)
All generated plots MUST strictly adhere to the following `matplotlib` configuration based on the user's provided example:
```python
plt.rcParams.update({
    'font.size': 12, 
    'font.family': 'Arial',
    'axes.linewidth': 1.5, 
    'xtick.major.width': 1.5, 
    'ytick.major.width': 1.5,
    'xtick.direction': 'in', 
    'ytick.direction': 'in', 
    'xtick.top': True,
    'ytick.right': True
})
```
*   **Labels**: Use `fontweight='bold'` for `xlabel` and `ylabel`.
*   **Legends**: `frameon=False`, `fontsize=12`.
*   **Math**: Use LaTeX formatting in strings (e.g., `r'Frequency (cm$^{-1}$)'`).
*   **Panel Labels**: Placed at top-left inside the panel `(a), (b), (c)` using `ax.text` with `fontweight='bold'`, `fontsize=14`.


## 2026-06-17 CHT Enhancements
- Added multi-peak support (num_peaks parameter) to the Complex Hankel Transform (CHT) fitting algorithm in nanoftir_Shizhe.py.
- Updated CHT plots in fitting_pipeline.ipynb to use 10^5 cm^-1 for momentum k to match FFT and real-space fits.
- Changed CHT plot annotations to output Peak 1, Peak 2, etc., instead of assuming Peak is approx 2q_p to avoid confusion.


## 2026-06-18 Batch Processing
- Automated insertion of CHT/Real-space/FFT template cells into fitting_pipeline.ipynb for all wavenumbers.
- Batch-generated and categorized all figures for 15 wavenumbers into the figures/ directory.

## 2026-06-23 Refactor + CHT k_fit_range diagnosis (Claude Code)
- **Refactor (no logic change)**: extracted the CHT fit+plot code that used to be copy-pasted in every
  one of the 15 per-wavenumber notebook cells into `fit_and_plot_cht()` in `nanoftir_Shizhe.py`, plus
  `load_aligned_wn_signal()` for data load/align/background-subtraction and `run_wn_comparison()` to
  orchestrate CHT + real-space (hankel/1-sqrtx) + FFT for one wavenumber, with an optional `save_dir` to
  write all three figures. Each wavenumber keeps its own notebook cell/markdown header (same outline,
  individually re-runnable); only the CHT cell shrank, from ~8000 chars to ~12 lines of parameters + one
  function call. Verified by re-executing the full notebook headlessly: all 15 wn reproduced the exact
  same lambda_p/damping/RMSE as before the refactor.
- **Root cause found for the high-wn CHT fit quality issue**: the CHT `k_fit_range_cm` had been left at a
  fixed (0.5, 6.0) for every wavenumber, but as `k_linked_guess_cm` grows with wn, the tip-launched (2*q_p)
  peak moves past k=6 and falls outside the fit window -- worst at 980/991/1000cm-1, where 991cm-1's CHT
  damping diverged to ~3641 (degenerate fit, no real signal in-window). Diagnostic |T(k)| plots saved to
  `figures/cht_diagnostics/`.
- User manually re-tuned `x_start_cht` / `k_fit_range_cm` / `k_linked_guess_cm` per wavenumber in the
  notebook (now generally widening/shifting the upper k bound with the guess, and narrowing the lower
  bound to dodge a spurious low-k peak). Deviation of CHT lambda_p vs the real-space Hankel benchmark
  dropped substantially for the worst cases (960cm-1: 17.3%->3.6%, 970cm-1: 10.7%->0.8%, 980cm-1:
  14.3%->4.6%); 991cm-1's damping is now a physical 5.2 instead of 3641. 1000cm-1, 890cm-1, 950cm-1 got
  slightly worse and may need another look (see CSV below).
- Re-ran `save_fit_results.py` with the user's tuned parameters to regenerate all 45 figures
  (`figures/cht/`, `figures/realspace/`, `figures/fft/`) and saved per-wn lambda/momentum/damping for
  CHT, Hankel, 1/sqrtx, and FFT (None where not applicable) to `data/fit_results.pkl`.
- Exported `data/cht_vs_realspace_wavelength_comparison.csv`: per-wn lambda_p (CHT/Hankel/1-sqrtx),
  the CHT k_fit_range_cm and x-distance fit range used, and CHT's % deviation from each real-space method.
- Cleanup: deleted stale one-off scripts superseded by the above (`append_cells.py`, `batch_run.py`,
  `dump_cells.txt`). `data/fitting_results_comparison.csv` (860cm-1-only, pre-refactor) and
  `figures/fit_result_860_*.png` / `figures/waterfall_gmg3.png` (pre-15-wn-pipeline) are now superseded
  by the files above but were left in place pending user confirmation before deleting tracked files.

## 2026-06-29 Manual linecut pipeline, shared snippet module, folder reorganization
- Built `manual_linecut_pipeline.ipynb`: per-wn manual linecut extraction from the richer
  `data/processed_3x1um/*.npz` maps (full O2/O3/O4 amp+phase + Z, not just one pre-extracted profile),
  with optional inline Z-plane correction (3-point fit, points drawn on the Topography panel), two-stage
  alignment, and CSV export to `data/graphene_3x1_manual/` in the same column format as the old data.
  Chosen per-wn (not one shared cut geometry) because horizontal align alone can't correct for drift
  *perpendicular* to the cut between separate scans -- confirmed real ~16px edge-position drift across
  the 15 npz scans, but only ~3px tilt within a single scan.
- Moved the reusable extraction/plotting code (Z-plane correction, `extract_and_plot`,
  `plot_waterfall_3channel`, `plot_mapping_waterfall`, `lineprofile_to_store`) out of any one project's
  notebook into `/Users/shizhe/envsetting/snippet/linecut_extraction.py` (+ `extract_linecut.py` for the
  underlying radial profile math) -- shared across GMG and G-CIPS-style projects instead of a
  slightly-different copy duplicated in each.
- Fixed: `plot_waterfall_3channel` no longer draws per-point markers by default (was rendering as a
  scatter cloud); `extract_and_plot`'s map aspect ratio bug (`aspect='auto'` stretched 3x1um maps into
  squares); Z-plane-correction point units (`plane_points_px` -> `plane_points_um`, matching the map's
  own axes instead of needing manual pixel conversion); `cmap_phase` default `'RdBu'` -> `'bwr_r'`.
- Built a v1 Streamlit web tool (`webtool/`) as an interactive alternative to typing every parameter --
  click on a map to place a linecut, click on the waterfall to align, click a spectrum to set q_guess.
  Has known rough edges (no headless browser available to test actual click interactions); not yet the
  primary workflow.
- Reorganized the folder: utility scripts -> `scripts/`, reference PDFs -> `docs/`, presentation decks
  -> `slides/`. `nanoftir_Shizhe.py`/`loadnpz.py`/both notebooks stay at the repo root (imported
  everywhere; moving them would mean editing notebook import cells for no benefit). Each moved script
  now `chdir`s to the repo root on startup so it runs correctly from any cwd.
- Found + fixed a real bug surfaced by the reorg: `save_fit_results.py` and `export_comparison_csv.py`
  hardcoded notebook cell indices to find each wn's CHT/real-space/FFT cells. Those went stale when the
  "Amp/Phase/Z overview" section was inserted into `fitting_pipeline.ipynb` earlier, silently
  mis-mapping cells for every wn except 860 (sits before the insertion point). Replaced with a dynamic
  lookup keyed off each wn's own `target_wn = '...'` cell. Re-verified against known-good lambda/q/damping.

## 2026-07-08 Multi-method fit comparison plots (CHT vs Hankel vs 1/√x)

- **Motivation**: with 3 real-space fitting methods (CHT, Hankel, 1/√x) all reporting a lambda_p per wn,
  it wasn't visually obvious which model actually tracks the data best -- only summary numbers were being
  compared (CSV/pptx). Added two new comparison figures per wn, styled after Woessner et al. 2015 SI
  Figure S4 (single-panel all-fits-overlaid + separate a/b/c/d panels with one map-location panel).
- **Baseline alignment (key design decision)**: CHT fits the background-subtracted oscillation (`sig_f`,
  mean ~0); `compare_cavity_models`'s Hankel/1-sqrtx fit whatever `y` column they're given with a free
  baseline `B`. In `fitting_pipeline_graphene_3x1_manual.ipynb`'s existing Cell 9, that `y` column is
  already `sig_f` (not raw O3A) -- confirmed by reading the cell -- so `B` comes out near 0 already; the
  three fits differ from the data by (at most) that small constant `B`. All three dense fit curves are
  plotted as `y_dense - B` (`y_dense - 0` for CHT, which has no free baseline) so everything sits on the
  same background-subtracted axis as the plotted data points, no additional background reconstruction
  needed. Self-checked by interpolating each dense fit curve onto the actual data x-positions and
  confirming the residual mean is small relative to the data's own std (CHT +0.020 vs data std 0.163,
  Hankel +0.016 vs 0.145, 1/sqrtx +0.020 vs 0.145 for 860cm-1) -- no systematic offset survives.
- **Q factor**: unified "damping"/"γp⁻¹" (already computed identically everywhere in the codebase as
  Re(q_p)/Im(q_p), using alpha_env as an Im(q_p) proxy for non-Hankel prefactors) under a single label
  "Q" in the new plots' legends, alongside λp and qp (in 1e5 cm⁻¹, matching the existing CSV/print
  convention). No new physics/formula -- purely a display-label unification per user request.
- **New code (non-destructive -- no existing function signature/behavior changed)**:
  - `nanoftir_Shizhe.py`: `fit_and_plot_cht`'s returned `results` dict gained 4 new keys
    (`x_dense_um`, `fit_dense`, `q_p_1e5cm_1`, `Q`) so its already-computed dense fit curve can be reused
    for overlay plots without re-fitting. `plot_fit_comparison_combined()` (single-panel, all 3 fits +
    data) and `plot_fit_comparison_panels()` (a/b/c/d layout, panel (a) optional) added at end of file.
  - `/Users/shizhe/envsetting/snippet/linecut_extraction.py`: `plot_linecut_location(ax, ...)` added --
    a single-Axes version of `extract_and_plot`'s amp-map + ROI-rectangle overlay, for embedding as one
    panel of a larger figure instead of its own standalone 3-panel figure.
- **Panel (a) (map + linecut location) data source**: the ROI geometry (`center_um`/`angle_deg`/
  `radius_um`/`rect_width_um`) used during manual extraction is not saved into the exported CSVs --
  it only exists in `manual_linecut_pipeline.ipynb`'s per-wn `extract_and_plot(...)` call. For the
  860cm-1 prototype, read directly from that notebook's cell: `center_um=(0.2, 0.8)`,
  `angle_deg=DEFAULT_ANGLE_DEG=0.0`, `radius_um=DEFAULT_RADIUS_UM=2.085`,
  `rect_width_um=DEFAULT_WIDTH_UM=0.075` (DEFAULT_* values come from that notebook's Cell 6 printed
  output). `map_panel=None` skips panel (a) entirely (falls back to a 1x3 layout) for datasets/wn where
  this geometry isn't available (e.g. built-in npz line profiles never manually placed).
- **Prototyped for one wn only so far**: new cells added to `fitting_pipeline_graphene_3x1_manual.ipynb`
  right after the existing 860cm-1 section (2 new code cells + 1 markdown header, nothing existing
  edited), producing `figures/graphene_3x1_manual/fit_comparison/860cm-1_{combined,panels}.png`.
  Verified by executing the actual (unmodified) notebook cells 0-16 headlessly end-to-end -- reproduces
  the exact same lambda_p/q_p/Q as the already-saved `data/fit_results_graphene_3x1_manual.pkl` (CHT
  451.6nm, Hankel 456.0nm, 1/sqrtx 403.1nm for 860cm-1). Not yet rolled out to the other 14 wn or to the
  4x1 lp1/lp2 datasets -- pending user review of the 860cm-1 prototype.

## 2026-07-08 Export linecut ROI geometry to CSV (replaces hardcoded/notebook-parsing lookup)

- **Problem**: the 860cm-1 prototype above got panel (a)'s geometry (`center_um`/`angle_deg`/
  `radius_um`/`rect_width_um`) by manually reading it out of `manual_linecut_pipeline.ipynb`'s source --
  doesn't scale to the other 14 wn or the 4x1 lp1/lp2 datasets, and silently goes stale if extraction
  parameters are ever re-tuned.
- **Considered writing the manual geometry back into the raw npz** (it already has an `lp1`/`lp2` slot
  and an `info.lineprofile.profiles[]` metadata block that could hold it). Rejected: that slot holds the
  *instrument/processing-software's own auto-extracted* profile (confirmed via `info.lineprofile.
  profiles[0].native_start_px/end_px` in `data/processed_3x1um/860cm-1_AVG.npz`) -- the whole reason
  manual extraction exists is that this auto placement often isn't good enough, so overwriting it would
  destroy the auto-vs-manual comparison and there'd be no way to recover the original. `np.savez` also
  has no in-place partial update (always rewrites the full archive), so every re-tune of any one wn would
  mean rewriting that wn's entire raw 2D-map file -- higher blast radius than corrupting a small sidecar,
  and opaque/undiffable in git (binary blob) vs a small text file.
- **Fix**: added a geometry-export cell (non-destructive append, after the existing CSV-export cell) to
  both `manual_linecut_pipeline.ipynb` and `manual_linecut_pipeline_graphene_4x1.ipynb`. Writes one
  combined `linecut_geometry.csv` per dataset (one row per wn, or per wn+edge for 4x1's lp1/lp2) into the
  same folder as the linecut CSVs (`data/graphene_3x1_manual/linecut_geometry.csv`,
  `data/graphene_4x1_manual/linecut_geometry.csv`) -- one file per dataset rather than one per wn, both
  to match the existing "one combined CSV, one row per wn" convention already used by the comparison CSVs
  and so the full history is visible in a single file's git log instead of 15+ separate file histories.
  Columns: `wn, lp, center_x_um, center_y_um, angle_deg, radius_um, rect_width_um, align_bg_loc_nm,
  realign_nm` (the last two document the alignment shift already baked into each CSV's `distance_nm`,
  free to include alongside the ROI geometry).
- Generated both files now by executing the two extraction notebooks headlessly end-to-end (Agg backend,
  stripped `%matplotlib`/`get_ipython()` magic calls that don't work outside a real kernel). Verified via
  `git status`/`git diff --stat` that this re-run reproduced the existing committed linecut CSVs
  byte-for-byte (fully deterministic given unchanged cell parameters) -- only the new
  `linecut_geometry.csv` appeared as untracked, nothing existing was altered.
- Updated `fitting_pipeline_graphene_3x1_manual.ipynb`'s 860cm-1 panel-(a) cell (added this session, not
  original pipeline code) to read `{data_dir}/linecut_geometry.csv` and look up the row matching
  `target_wn`, instead of the hardcoded tuple from before -- falls back to `map_panel = None` if the file
  or that wn's row is missing. Re-verified headlessly: reproduces the identical `center_um=(0.2, 0.8),
  angle_deg=0.0` panel (a) as the hardcoded version.

## 2026-07-08 Cross-project index + dispersion-fitting engine migrated to shared snippet/

- Set up cross-project infrastructure after discovering `nanoftir_Shizhe.py` (the CHT/Hankel/
  1-sqrtx/FFT dispersion-fitting engine, originally authored by Yinming Shao) had silently
  forked into 3 diverged copies across GMG (most current), CrSBr_Analysis (634 lines behind),
  and an older ancestor `nanoftir.py` in G-CIPS. Added `~/.claude/CLAUDE.md` (global, always
  loaded) pointing at `/Users/shizhe/envsetting/PROJECTS.md` (index of all of Shizhe's
  experimental-data-analysis projects and what each depends on in the shared `snippet/`
  library) and `dispersion_fitting_migration_plan.md` (the migration plan, written so a future
  session in any directory can pick it up without relying on this conversation's memory).
- **Phase 1 executed**: extracted `fit_cavity_prefactor_compare`, `compare_cavity_models`,
  `plot_channel_fft`, `complex_hankel_transform`, `fit_cht_peaks`, `fit_and_plot_cht`,
  `plot_fit_comparison_combined`, `plot_fit_comparison_panels` out of `nanoftir_Shizhe.py` into
  `/Users/shizhe/envsetting/snippet/dispersion_fitting.py`. `load_aligned_wn_signal` and
  `run_wn_comparison` stayed here (GMG-specific `data_dir` defaults). `nanoftir_Shizhe.py`
  shrank 4028 -> 2818 lines and now re-exports the moved names
  (`from snippet.dispersion_fitting import (...)`) -- every existing notebook cell
  (`nanoftir.fit_and_plot_cht(...)` etc.) kept working with zero edits.
- **Verification**: re-ran all 3 `save_fit_results_graphene_*.py` scripts post-migration and
  diffed every value in the resulting `.pkl` files against the pre-migration versions -- zero
  numeric difference across every wn/method. Re-executed the 860cm-1 fit-comparison prototype
  cells headlessly, identical result. Caught and fixed a real bug during the move itself (not
  in the shipped result): the line-deletion script initially processed two overlapping ranges
  out of order, which would have deleted the wrong 393-line block (dropping `plot_stacked_fft`,
  which should have stayed, while leaving part of `plot_channel_fft`, which should have moved) --
  caught by the syntax + remaining-function-list check before any commit; reverted via
  `git checkout` and redone with ranges sorted strictly descending plus an assertion.
- **Not done**: Phase 2 (CrSBr_Analysis / G-CIPS adopting the shared module) -- their copies
  have genuinely diverged, not just fallen behind, so this needs a function-by-function diff
  audit first. See `dispersion_fitting_migration_plan.md` Phase 2 for the plan.

## 2026-07-08 Critical bug: lp2 pipeline was fitting lp1 data

- Found while building batch fit-comparison figures for graphene_4x1: `load_aligned_wn_signal`'s
  glob pattern was hardcoded to `*_AVG_lp1.csv`. Since `data/graphene_4x1_manual/` holds both
  lp1 and lp2 CSVs for every wn, `save_fit_results_graphene_4x1_lp2.py` had been silently
  loading and fitting **lp1 (left edge) data** with lp2-tuned CHT/real-space parameters the
  entire time -- the "lp2" pkl/CSV/PPTX delivered earlier this session did not describe the
  right edge at all.
- This explains the previously-flagged oddity: 1000cm-1's CHT-vs-Hankel deviation was 60.1%
  under the bug, now 6.5% after the fix -- consistent with the mismatch being wrong-data, not a
  real fit-quality problem.
- Fix: added `lp='lp1'` parameter to `load_aligned_wn_signal` and `run_wn_comparison` (default
  preserves old behavior for single-linecut datasets), glob pattern now `*_AVG_{lp}.csv`.
  `save_fit_results_graphene_4x1_lp{1,2}.py` now pass `lp` explicitly.
- Verified: lp1 rerun is byte-identical (0 values changed -- it was already reading the right
  file, just not by design). lp2 rerun: 157 of ~225 values changed, largest single lambda_p
  delta ~181nm. Regenerated lp2's CSV/PPTX from the corrected pkl. Some wn (980/970/960/911/900)
  still show large CHT-vs-Hankel deviation (42-48%) even with correct data -- worth revisiting
  those k_fit_range/xr tunings separately, not something this fix addresses.
