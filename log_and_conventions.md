# Progress Log

Stable folder structure + conventions/rules moved to `CLAUDE.md` (2026-07-10) -- this file is
now the dated history only (why things are the way they are), not auto-loaded every session.
Read it deliberately when historical context is needed.

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

## 2026-07-08 Im(r_p) dispersion maps: extracted q_p vs simulated multilayer reflection

- New `rp_dispersion_plot.ipynb`: for each of the 3 linecut datasets (graphene_3x1_manual =
  "Gr/SiO2(300nm)/Si"; graphene_4x1_manual lp1 = "MoO3(4.2nm)/Gr/SiO2/Si -- edge"; lp2 = "...--
  interior"), computes Im(r_p)(q, omega) via `snippet.nearfield.rp_schematic`
  (`NearFieldOptics.Materials.AnisotropicMaterial` + `LayeredMediaTM`, same engine
  `response_Shizhe.ipynb` uses for doped graphene/CIPS and `CrSBr_Rp_calculation` uses for
  CrSBr) and overlays the CHT/Hankel/1-sqrtx/FFT q_p points from `data/fit_results_graphene_*.pkl`.
- MoO3 modeled as `M.AnisotropicMaterial(eps_infinity=[4,5.2,2.4], phonon_params=[...])` per
  https://onlinelibrary.wiley.com/doi/abs/10.1002/adma.201908176 (in-plane x=y: wLO=972 gLO=4
  wTO=820 gTO=4; out-of-plane z: wLO=1004 gLO=2 wTO=958 gTO=2) -- confirmed the library's
  `get_phonons()` expects `(w_LO, g_LO, w_TO, g_TO)` per mode, which is already the order the
  reference paper's numbers were given in (no reordering needed, unlike CrSBr_Rp.py's own
  `convert_modes()` which had to reorder from a different raw tuple convention).
- **E_F fit per method per dataset (12 fits)**, not one E_F per dataset: initial attempt used a
  single E_F per dataset (from the user's reference table: edge=0.65eV, interior=0.62eV,
  graphene-only=0.48eV) applied to all 4 methods' panels -- only the method that table value
  actually corresponded to (apparently CHT) hugged the ridge; Hankel/1-sqrtx/FFT visibly missed
  it, since each method's own q_p extraction doesn't agree exactly at a given wn. Switched to
  fitting E_F independently per method: for a trial E_F, evaluate Im(r_p) at the exact (q_i, w_i)
  data points (diagonal of the freq x q grid `reflection_p` returns, converted to plain
  `np.ndarray` first -- passing the library's native AWA-wrapped array straight into
  `scipy.optimize.minimize_scalar` throws `IndexError` from `common/baseclasses.py`'s
  `inheriting_operator`), sum, maximize via bounded 1D `minimize_scalar`. gamma=20cm-1 fixed
  throughout (per the reference figures). Results: 3x1 {cht:0.51, hankel:0.50, sqrtx:0.50,
  fft:0.65}, lp1/edge {cht:0.77, hankel:0.77, sqrtx:0.77, fft:0.70}, lp2/interior {cht:0.72,
  hankel:0.92, sqrtx:0.96, fft:0.51} -- close to but not identical to the user's reference table
  (different fitting method/objective, not expected to match exactly).
- **Self-check finding (not a bug)**: 3x1_manual and lp1(edge) show all 4 methods' points
  tracking the ridge cleanly after the per-method E_F fit. lp2(interior) does not -- points at
  950-990cm-1 sit systematically off the ridge across all 4 methods, while 860-940 and 1000cm-1
  track fine. This lines up exactly with the 2026-07-08 entry above flagging
  980/970/960/911/900cm-1 as still having large CHT-vs-Hankel disagreement even post-bugfix --
  independent confirmation via a completely different (rp-simulation-based) method that those
  specific wn's fits are the weak point, not new noise introduced by this notebook.
- `ytickspacing=10` passed to `rp_schematic` (default 200cm-1 -- with our ~160cm-1-wide freq
  window that rendered only one labeled y-tick). Later widened to `ytickspacing=20` per user
  feedback (10 was too dense).

## 2026-07-09/10 Restyled compare_cavity_models; new fit-comparison overview PPTX

- **`compare_cavity_models` restyle** (`snippet/dispersion_fitting.py`): matched the visual
  style of `plot_fit_comparison_combined`/`plot_fit_comparison_panels` -- gray data dots (was
  black), per-curve lambda_p/q_p/Q folded into the legend (was a separate stats textbox), bold
  panel letters (a)/(b), bold axis labels/title. `ylim` default changed from a stale hardcoded
  `(0.5, 0.9)` to `None` (auto-scale) -- grepped all 65 call sites across every GMG
  notebook/script and confirmed every one already passes an explicit, per-wn-computed `ylim`
  (e.g. `(sig_f.min()*1.5, sig_f.max()*0.5)`), so xlim/ylim was *already* being judged
  automatically from each wn's own data, never a fixed standard or manually eyeballed value --
  same answer applies to `plot_fit_comparison_combined`/`panels`, which also default to
  `xlim=ylim=None` (matplotlib auto-scale) and were never called with an explicit override.
  Purely cosmetic -- `fit_cavity_prefactor_compare` (the fitting math) untouched, return
  signature `(outs, fig, (ax, axr))` unchanged. Verified by re-running all 3
  `save_fit_results_graphene_*.py` scripts and diffing every pkl value pre/post -- zero
  numeric difference. Regenerated the 45 affected `realspace/{wn}_realspace.png` figures.
- **New `scripts/make_fit_comparison_pptx.py`** -> `slides/GMG_fit_comparison_overview.pptx`
  (49 slides, font Aptos): title + one section-divider per dataset + one slide per wn
  (`{wn}cm-1_combined.png` on top, `{wn}cm-1_panels.png` below). Caught a real layout bug in
  self-check before shipping: a fixed-guess `top_h=3.55in` for the top picture, combined with
  the bottom picture's own aspect ratio at full content width, pushed the bottom image 0.13in
  past the slide edge on every wn slide -- no PPTX renderer available in this environment, so
  verified via python-pptx shape-bounds introspection (`shape.left/top/width/height` vs
  `prs.slide_width/height`) across all 195 shapes in the deck instead of visual inspection.
  Fixed by deriving `top_h` from slide height minus the aspect-locked bottom picture height
  rather than guessing a fixed number.

## 2026-07-10 v2 retune comparison for graphene_4x1 lp1 (edge) [logged retroactively 2026-07-13]

- User re-tuned `xr_range` (Hankel/1-sqrtx fit window) for most wn, plus `x_start_cht`/
  `L_cutoff_cht` for 860cm-1, directly in `fitting_pipeline_graphene_4x1_lp1.ipynb` --
  **lp1 (edge) only**, not lp2 or 3x1_manual.
- New `scripts/retune_compare_graphene_4x1_lp1.py`: parses the notebook's *current* tuned
  params (same `target_wn = '...'` anchor-cell regex used by the other `save_fit_results_*`
  scripts), reruns the fit, and writes into separate `_v2`-suffixed paths (`data/
  fit_results_graphene_4x1_lp1_v2.pkl`, `figures/graphene_4x1_manual/lp1_v2/`, `figures/
  rp_dispersion_v2/`) -- does **not** touch/overwrite the original non-v2 outputs, so old
  and new sit side by side. `VERSION_SUFFIX` constant at the top, bump before each new
  retuning round; script raises `FileExistsError` if the target pkl already exists (refuses
  to silently clobber a previous comparison round).
- **Self-check** (diffed old vs new lambda_p + inter-method agreement, all 15 wn): net
  improvement for most wn, e.g. 941cm-1 agreement 25.5%->2.3%, 860cm-1 4.1%->0.8%. **1000cm-1
  got worse** (1.4%->32.0%) -- narrowing `xr_range` there made Hankel/1-sqrtx lock onto a
  shorter spurious oscillation near the noisy edge region instead of the dominant one.
  Flagged to user at the time, left unresolved -- 1000cm-1's `xr_range` likely needs
  separate attention.
- **v2 was never merged into the main (non-v2) pkl/figures** -- it's a side-by-side
  comparison snapshot only, not adopted. The live `data/fit_results_graphene_4x1_lp1.pkl`
  and its CSV export (`data/cht_vs_realspace_wavelength_comparison_graphene_4x1_lp1.csv`,
  generated by `scripts/export_comparison_csv_graphene_4x1_lp1.py`) are both **still the
  original, pre-retune (v1) numbers**, untouched since 2026-06-30 (commit `285eced`) --
  confirmed via `git log` showing the `fb772bb` retune commit only *added* new v2-suffixed
  files, never modified these. The pre-retune parameter *values* themselves (the old
  `xr_range` etc., as opposed to the results they produced) are no longer visible in the
  live notebook (its cells were edited in place to the new values) but remain fully
  recoverable from git history: `git show 02f9a2b:fitting_pipeline_graphene_4x1_lp1.ipynb`
  is the version immediately before the retune commit `fb772bb`.
- Commit message (`fb772bb`) also notes `rp_dispersion_plot.ipynb`'s saved execution
  outputs were swept in from a normal interactive run in the same session -- unrelated to
  the retune itself, no code changes, no errors.
- **Process gap noted**: this change was documented in the commit message but not in this
  Progress Log at the time, despite CLAUDE.md's Extensive Logging Policy -- backfilled now.

## 2026-07-13 MoO3 a-axis vs c-axis in-plane parameter comparison in `rp_dispersion_plot.ipynb`

- Question raised: our MoO3 `AnisotropicMaterial` (`eps_infinity=[4,5.2,2.4]`,
  `phonon_params=[modes_xy, modes_xy, modes_z]`, modes_xy=(972,4,820,4)) previously cited
  `10.1002/adma.201908176` -- checked user-supplied papers in `docs/` instead (Zheng et al.
  2019, *Sci. Adv.* 5, eaav8690, "A mid-infrared biaxial hyperbolic van der Waals crystal";
  its own SI `docs/aav8690_sm.pdf` has Table S1, a "this work" column with the same x/z
  numbers plus the previously-missing y column). Confirmed alpha-MoO3 crystallographic
  convention: x=[100]=a-axis, y=[001]=c-axis, z=[010]=b-axis (out-of-plane, vdW stacking
  direction). Table S1 "this work": a-axis eps_inf=4.0/LO=972/TO=820/Gamma=4, c-axis
  eps_inf=5.2/LO=851/TO=545/Gamma=4, b-axis eps_inf=2.4/LO=1004/TO=958/Gamma=2.
- **Confirmed our (and an independently-written colleague copy of the same) `phonon_params`
  reuses the a-axis's LO/TO for both the x and y slots, while `eps_infinity` keeps each
  axis's own value (4.0 vs 5.2)** -- user confirmed this is intentional, not a typo.
  Traced `NearFieldOptics.Materials.material_types.BaseAnisotropicMaterial.epsilon()` ->
  `ordinary_epsilon()` -> `ordinary_component()` = `numpy.mean(diag[:2])`: `reflection_p`'s
  underlying uniaxial TM solver has no azimuthal-angle parameter and only ever uses the
  x/y **average** of the principal epsilon tensor, so a true direction-dependent biaxial
  response isn't representable here regardless of how many distinct axis values are fed in
  -- only a uniaxial (x=y) choice is physically meaningful. With x=y sharing one phonon
  lineshape but different eps_infinity, the effective in-plane response is mathematically
  a single uniaxial material with eps_infinity=mean(4.0,5.2)=4.6 and whichever axis's
  LO/TO/Gamma was chosen for both slots -- not literally "a-axis" or "c-axis" alone.
- Since the true crystalline orientation of this flake relative to the linecut direction is
  not currently known, extended `rp_dispersion_plot.ipynb` to build both `MoO3_a`
  (`phonon_params=[phonon_a, phonon_a, phonon_b]`) and `MoO3_c` (`phonon_params=[phonon_c,
  phonon_c, phonon_b]`), fit E_F and render the 4-panel dispersion figure for both, for lp1
  and lp2 (graphene_3x1_manual unaffected, no MoO3 in that stack). Titles now read
  `MoO3(a)...`/`MoO3(c)...` so the two are unambiguous at a glance. Old unsuffixed
  `graphene_4x1_lp{1,2}_*_rp_dispersion.png` (the previous, unlabeled a-axis-only output)
  removed, superseded by `..._MoO3a.png`/`..._MoO3c.png`.
- **Self-check**: MoO3(a) E_F values reproduce the previously-committed ones exactly (lp1:
  cht/hankel/sqrtx/fft = 0.771/0.773/0.768/0.704, matching pre-change to 3 decimals) --
  confirms this was a pure extension, not a change to the existing pathway.
- **Result**: switching a->c does NOT bring E_F closer to the Raman-measured ~0.6-0.65eV --
  it moves further away in every case. lp1/edge: cht 0.771->0.839, hankel 0.773->0.829,
  sqrtx 0.768->0.837, fft 0.704->0.735. lp2/interior: cht 0.723->0.738, hankel
  0.915->0.960, sqrtx 0.960->0.959, fft 0.510->0.795. So the a/c ambiguity is not the
  explanation for the E_F-vs-Raman gap either way; c-axis is the worse of the two options
  by this metric. Axis choice (which one reflects the actual flake orientation) still
  unresolved -- both figures are kept side by side so this can be revisited once/if the
  crystal orientation is independently determined (e.g. polarized Raman per Zheng et al.
  2019 Note S4).
- Added a 7th section to `rp_dispersion_plot.ipynb`: a single-panel Im(r_p) vs q lineout at
  fixed omega=900cm-1 and a single fixed E_F=0.65eV (top of Raman range, common reference
  point -- deliberately NOT any method's own fit, since per-method E_F fits partly
  compensate for the a/c difference and would obscure how much the axis choice shifts the
  dispersion by itself), MoO3(a) solid black vs MoO3(c) dashed blue. First version also
  overlaid the extracted q_p at 900cm-1 (cht/hankel/sqrtx/fft) as vertical lines, but
  cht/hankel sit almost exactly on top of each other there (q=1.20 both) making inline
  labels illegible, and per user feedback the data overlay was dropped entirely to keep
  the figure to just the axis-choice comparison -- simpler reads better than a fixed
  corner text box here. Result: MoO3(c)'s peak sits at visibly larger q than MoO3(a)'s
  and is slightly broader, at the same fixed E_F -- the direct geometric reason a fixed
  q_p needs a higher fitted E_F under the c-axis assumption. Peak position marked with a
  dot and its lambda_p = 2*pi/q_peak (nm) put directly in the legend label per curve:
  MoO3(a) lambda_p=439nm vs MoO3(c) lambda_p=409nm at this wn/E_F. Saved to `figures/
  rp_dispersion/imrp_vs_q_900cm-1_MoO3_axis_compare.png`.
- Extended the same cell to also overlay each wn's actual FFT spectrum on a `twinx()` right
  axis, and to loop over 860/900/980cm-1 (one saved PNG per wn) instead of just 900. FFT
  spectrum is the **complex-FFT** (`plot_channel_fft(..., plot=False)`'s `q_complex`/
  `spec_complex` return values) -- confirmed this is the same one `run_wn_comparison` uses
  for the `fft` method's own `q_p_1e5cm-1` in the pkl (`nanoftir_Shizhe.py:2816`,
  `fft_peaks = fft_out['peaks_complex']['peaks']`), not the amp-FFT or phase-FFT panels
  `plot_channel_fft` also computes. Per-wn `xr`/`q_guess` FFT params copied verbatim from
  `fitting_pipeline_graphene_4x1_lp1.ipynb`'s per-wn FFT cells (anchor-cell-relative offset
  +6, same pattern used elsewhere for CHT/real-space params): 860cm-1 xr=(0.15,1.2)
  q_guess=[1,2]; 900cm-1 xr=(0.11,1.2) q_guess=[1.2,2.2]; 980cm-1 xr=(0.14,1.2)
  q_guess=[3.5,7.0]. Raw data read directly from `data/graphene_4x1_manual/{wn}cm-1_AVG_lp1.csv`
  (`distance_um = distance_nm/1000`, matches `aligned_data_dict` construction in the fitting
  pipeline notebook since align_dict shift=0 for all lp1 wn). Saved to `figures/
  rp_dispersion/imrp_vs_q_{860,900,980}cm-1_MoO3_axis_compare.png`.
- **Follow-up finding, prompted by user noticing the raw FFT peak looked ~2x too high
  compared to the Im(r_p) ridge**: verified numerically (not guessed) that at all three wn,
  the FFT spectrum's dominant/global peak is consistently ~2x the independently-established
  q_p from CHT/Hankel/1-sqrtx (which agree tightly with each other at each wn) -- 860cm-1:
  global peak=1.80 vs 2x(0.93-0.96)=1.86-1.92; 900cm-1: peak=2.49 vs 2x(1.17-1.20)=2.34-2.40;
  980cm-1: peak=4.74 vs 2x(2.36-2.38)=4.72-4.76 (near-exact). Smaller secondary peaks further
  out (860cm-1: 4.20; 900cm-1: 4.99=2x2.49) are consistent with windowing/padding sidelobes
  or higher harmonics (4q), not a resolvable independent "q" feature.
  **Conclusion: the FFT spectrum's strongest feature is the "2q" (tip-launched round-trip)
  channel, not "q" (edge-launched)** -- direct spectral confirmation of the q/2q
  coexistence hypothesis raised earlier this session (see the two-term Hankel prototype
  entries above), this time via the existing FFT machinery rather than a new model.
  This also explains a chunk of `fft`'s known unreliability as a method throughout this
  project: checked what `fft`'s own pkl-adopted q_p (nearest `q_guess[0]`, "q not 2q"
  convention) actually corresponds to at these 3 wn -- at 860/900cm-1 it isn't a genuine
  local peak at all (the spectrum is monotonically rising through the whole `q_guess[0]`
  search window toward the dominant 2q peak just outside it; the windowed search returns a
  point near the window's edge, not a real spectral feature); at 980cm-1 it lands almost
  exactly on the true 2q peak itself (4.538 adopted vs 4.735 actual 2q peak) -- i.e. `fft`'s
  "q_p" there is 2q mislabeled as q. Not a uniform "off by exactly 2x" error across all
  three -- 860/900cm-1 are window-boundary artifacts with no clean multiplicative
  relationship to either q or 2q, 980cm-1 is a clean "wrong channel" mislabeling.
- Dropped the CHT `|T(k)|` overlay from the plot (tried, then removed): its near-k=0
  log-divergence (a real mathematical property of the Hankel-transform kernel, not a bug)
  dominates the curve, and away from that divergence the fitted q_re region only shows a
  faint shoulder, not a resolvable peak (CHT's q_p comes from a joint least-squares fit
  across the transform, not from peak-finding on `|T(k)|`, unlike FFT) -- showing it
  alongside FFT's peak-having spectrum invited exactly the wrong comparison. Replaced with:
  keep the raw FFT spectrum (still informative, shows the 2q/4q structure) and mark the
  identified dominant "2q" peak with a dashed vertical line + numeric label
  (`2q={value}`) instead of any CHT curve or the (now known unreliable) `fft`-adopted pkl
  marker. Re-saved to the same three `imrp_vs_q_{860,900,980}cm-1_MoO3_axis_compare.png`
  paths.
- **Tested and rejected a candidate fix (STATUS: UNRESOLVED, do not act on "peak/2" without
  re-reading this)**: checked whether pkl's `fft_q_p` should just be redefined as
  (global-max of the raw complex-FFT spectrum)/2, across all 15 lp1 wn, not just the 3
  tested above. Works well for 860/870/880/900/980cm-1 (peak/2 vs established
  CHT/Hankel/sqrtx average, ratio 0.96-1.11) but **fails badly for 920/930/941/960/1000cm-1**
  (ratio 0.04-0.07) -- for these wn the spectrum's global max sits at the very first
  non-zero frequency bin (q~0.19-0.21, essentially DC), a low-frequency leakage/background
  artifact unrelated to q or 2q (confirmed: amplitude=1.0 at the first bin then
  monotonically decreasing for several bins before any real structure). 911/891/970/950/991
  are intermediate/mixed (ratio 0.46-1.11, no clean pattern). Restricting the search to
  q>0.5 (1e5cm-1) doesn't fix it either -- 930/941cm-1 then land near the true "q" itself
  (not 2q) instead, 960/1000cm-1 still find noise, no single restriction rule works
  uniformly. **Conclusion: no blanket formula recovers a reliable fft q_p across all 15
  wn -- would need per-wn manual attention (re-tuned xr/q_guess and/or better background
  subtraction), comparable in effort to the original per-wn CHT/real-space tuning. Left
  pkl and all downstream figures untouched.** CHT/Hankel/1-sqrtx remain the trustworthy
  methods; fft's role should be treated with active suspicion, not silently "corrected",
  until someone does that per-wn work.

## 2026-07-13 Peak-spacing (skip-first-peak) real-space cross-check, wn=860-911cm-1

- User re-tuned real-space `xr_range`/`lam0_guess_um` for wn 860/870/880/890/900/911 in
  `fitting_pipeline_graphene_4x1_lp1.ipynb` (CHT params confirmed unchanged via git diff for
  all 6). Motivation: noticed a resolvable "q" (not just "2q") peak starting to appear in the
  FFT spectrum around 911cm-1 (section 7 of `rp_dispersion_plot.ipynb` -- below 900cm-1 the
  FFT is 2q-dominated). Wanted a real-space cross-check: find oscillation peaks directly,
  skip the first one (closest to the tip/edge, possibly contaminated by near-field
  artifacts), then compare two interval definitions: spacing(peak2,peak3) alone vs
  mean(spacing(peak1,peak2), spacing(peak2,peak3)).
- New `scripts/peak_spacing_crosscheck_graphene_4x1_lp1.py` (VERSION_SUFFIX pattern,
  prototype/investigation only, doesn't touch the committed pkl): refits CHT+Hankel+1sqrtx
  with the current (just-retuned) params, finds real-space peaks via `scipy.signal.find_peaks`
  within the tuned `xr_range` (3-4 peaks found per wn, enough to skip-first + compute both
  intervals), and overlays all 5 q sources (cht/hankel/sqrtx/peak2-3/avg) on one Im(r_p)
  MoO3(a) map using a **fixed reference E_F=0.771eV** (the already-established full-15-wn
  CHT fit -- not refit per source, 6 points is too thin). Output:
  `figures/rp_dispersion/peak_spacing_crosscheck_v1/` (CSV + combined rp png).
- **Result: skipping the first peak did NOT resolve a genuine "q" periodicity.** Both
  interval definitions land at ~2x the CHT/Hankel/sqrtx q_p at all 6 wn (e.g. 911cm-1:
  established~1.32-1.49, peak2-3=2.72, avg=3.03) -- visually confirmed on the combined rp
  plot: cht/hankel/sqrtx sit right on the ridge, peak2-3/avg sit far to the right of it,
  consistently ~2x. The 2q (tip-launched round-trip) channel still dominates the real-space
  oscillation throughout this window even after dropping the near-edge peak; 911cm-1's new
  FFT-visible "q" shoulder isn't yet strong enough to show up as the dominant real-space
  period via simple peak-counting.
- **Bug fixed during this work**: first version of the combined rp plot had all points
  collapsed near q=0 -- forgot to convert the peak-spacing/CHT/Hankel/sqrtx q values (in
  1e5 cm^-1) to raw cm^-1 before plotting on `rp_schematic`'s axes (which use raw cm^-1
  internally despite displaying 1e5 cm^-1-scaled tick labels). Fixed by multiplying by 1e5
  before the `ax.plot()` call; re-verified the resulting point positions visually track the
  ridge for cht/hankel/sqrtx as expected.
- **New unrelated finding (STATUS: UNRESOLVED, narrow scope)**: while refitting CHT for
  self-consistency, found the freshly-computed CHT q_p exactly matches the committed pkl for
  5 of these 6 wn (870/880/890/900/911, 0.00% diff) but **disagrees for 860cm-1 by -6.8%**
  (pkl=0.962, fresh=0.896) despite identical CHT parameters (x_start_cht/L_cutoff_cht/
  k_fit_range_cm/k_linked_guess_cm all confirmed unchanged via git diff against the original
  committed notebook cell). Reproducible/deterministic on the fresh side (re-ran twice,
  identical). Likely explanation: the committed pkl predates the 2026-07-08 CHT-engine
  migration to `snippet/dispersion_fitting.py` (pkl committed 2026-06-30, never regenerated
  since) and 860cm-1 was already flagged elsewhere this session as a wn where the CHT/
  real-space least-squares fit is unusually sensitive to a q/2q-degenerate local minimum --
  plausibly the migration shifted the optimizer's basin of attraction for this one marginal
  wn without affecting the other 14. Not investigated further -- flagging that the full
  lp1 (and possibly lp2/3x1_manual) pkl may be worth regenerating from current code and
  diffing against the committed version to check for other silent drifts, given this was
  found by accident rather than by a deliberate post-migration re-verification.
- **Follow-up on the same plot**: re-ran with (a) CHT dropped (redundant with Hankel/sqrtx
  here, down to 4 sources), (b) peak-spacing sources `/2`-corrected (confirmed 2q, now
  comparable scale to q), (c) E_F switched from the previous 0.771eV (full-15-wn CHT fit)
  to **0.65eV, the Raman-measured value for this dataset (graphene_4x1 lp1 = "edge")**.
  Result: Hankel/sqrtx track the E_F=0.65eV ridge fairly well across 860-911cm-1 (visually
  close); the /2-corrected peak-spacing sources land in the right ballpark (no longer 2x
  off) but sit consistently slightly to the right of the ridge (a small, systematic
  overestimate of q, not random noise) at every one of the 6 wn -- the /2 correction gets
  the scale right but doesn't fully reproduce the established fits or exactly hit the
  Raman-consistent ridge.
- **Cleanup per user feedback**: (1) legend label fixed from "peak2-3" to "peak3-peak2"
  (matches how the interval is actually computed, `xp[2]-xp[1]`, later minus earlier);
  (2) output folder renamed from the dataset/linecut-ambiguous
  `figures/rp_dispersion/peak_spacing_crosscheck_v1/` to
  `figures/rp_dispersion/graphene_4x1_manual_lp1_peak_spacing_crosscheck_v1/`; (3) now
  generates 3 panels, one per E_F (0.60eV Raman lower bound, 0.65eV Raman edge value,
  0.771eV previous full-15-wn CHT fit), same 4 sources each. Finding: the offset between
  hankel/sqrtx and the /2-corrected peak-spacing sources is essentially **E_F-independent**
  (changing E_F reshapes/shifts the ridge but doesn't change the already-extracted q
  values) -- hankel/sqrtx track the ridge reasonably at all 3 E_F choices in this narrow
  860-911cm-1 window (not much ridge curvature here vs. the full 15-wn range), while the
  peak-spacing/2 sources sit consistently to the right of the ridge at all 3, confirming
  the residual overestimate is intrinsic to the peak-spacing method itself, not an
  artifact of picking the wrong E_F reference.
- Added a 4th E_F panel, 0.67eV, per user request -- same pattern holds (hankel/sqrtx on
  the ridge, peak-spacing/2 sources consistently to its right).
- Added per-source E_F fit (reused rp_dispersion_plot.ipynb's "maximize sum Im(r_p) at
  data points" recipe, `scipy.optimize.minimize_scalar`) on the 6-point subset, to
  directly answer "which doping level is correct" instead of eyeballing fixed-E_F panels.
  Also fixed the `avg(1-2,2-3)` legend label to `avg(2-1,3-2)` (same "later-minus-earlier"
  convention as the peak3-peak2 fix).
- **Root cause found for the "870/900 look off-trend" complaint (bigger than the pixel-
  quantization issue already logged)**: peak-finding was restricted to `xr_range_rs`, but
  that window's lower bound already excludes the true first (tip/edge) oscillation peak --
  so what the script called "peak1/2/3" within that window were actually the physical
  2nd/3rd/4th peaks, and "peak3" (used for the interval) was the 4th real peak, deep in
  the damping tail and barely above the noise floor (e.g. 900cm-1's real peak4 amplitude
  ~0.09 vs peak1's ~4.2). **Fixed**: find peaks on the FULL `x_f`/`sig_f` (not restricted
  to `xr_range_rs`), so peak1 is the genuine tip/edge peak. Had to also fix the prominence
  threshold: scaling it to the *full* signal's peak-to-peak range let the huge first peak
  drown out real peak2/peak3 (860/900/911 dropped to only 2 detected peaks); fixed by
  scaling prominence to the ptp of the tail only (`x>0.25um`, past every wn's first peak
  here) while still searching the full range. Visually verified (2x3 diagnostic grid,
  scratchpad) peak1 now lands on the actual tallest early peak for all 6 wn, peak2/peak3
  on the next two real oscillation peaks.
- **Result improved substantially** after the fix: peak-spacing/2 points now sit close to
  the hankel/sqrtx points and the ridge (previously used the wrong, noisier 4th-peak
  interval and sat visibly right of the ridge at all 6 wn). Refitted per-source E_F:
  hankel=0.695eV, sqrtx=0.700eV, peak3-peak2/2=0.683eV, avg(2-1,3-2)/2=0.749eV -- all 4
  now cluster in 0.68-0.75eV (vs Raman ~0.60-0.65eV), with the peak-spacing variants
  bracketing hankel/sqrtx rather than sitting off to one side. Re-saved all 4 E_F panels
  + CSVs to the same `graphene_4x1_manual_lp1_peak_spacing_crosscheck_v1/` folder.
- **Simplified to peak2-peak3 only, sub-pixel refined**: per user feedback after visually
  checking 860/911cm-1 (peak2/peak3 confirmed genuine, but intervals were irregular --
  860cm-1: interval(1-2)=369nm vs interval(2-3)=292nm, 21% different; 911cm-1: 277nm vs
  185nm, 33% different, non-monotonic even out to peak4 -- consistent with peak position
  jitter from low SNR plus possible real q/2q beating, not a peak-finding bug), dropped the
  `avg(1-2,2-3)` source entirely (peak2-peak3 alone is the more reliable pair) and added
  3-point parabolic sub-pixel interpolation around each `find_peaks` index (`subpixel_peak_x`)
  instead of using the raw pixel-grid maximum -- this also removes the pixel-quantization
  artifact logged earlier (870/880cm-1 no longer land on the exact same interval value).
  Result: 880/890/900cm-1 now match hankel/sqrtx to within 3%; 860/911cm-1 (the two lowest-
  SNR/most irregular cases) still sit ~14% high, consistent with the low-SNR/beating
  explanation being a real signal limitation rather than a fixable peak-picking error.
  Fitted E_F now: hankel=0.695eV, sqrtx=0.700eV, peak3-peak2/2(subpixel)=0.683eV -- tighter
  cluster (0.68-0.70eV) than the pre-fix 3-source spread (0.68-0.75eV). Now only 3 sources
  (hankel/sqrtx/peak3-peak2 fit) plotted per E_F panel; same 4 panels re-saved.
- **Switched sub-pixel refinement from 3-point parabola to Lorentzian fit**, per user
  request ("should use Lorentzian, not quadratic"). `subpixel_peak_x` now fits
  `B + A/(1+((x-x0)/gamma)^2)` via `scipy.optimize.curve_fit` on a +-6-point window around
  each `find_peaks` index (falls back to the raw grid position if the fit doesn't
  converge). Visually verified fit quality first (all 6 wn x peak2/peak3, 12 panels,
  `scratchpad/lorentzian_peak_fit_preview.png`) before adopting -- curves track the local
  peak shape well, fitted center close to but distinct from the raw grid maximum in every
  panel. Result: even tighter than the parabolic version -- 870/880/890/900cm-1 now
  visually overlap hankel/sqrtx on the rp plot; 860/911cm-1 (lowest SNR) still sit
  slightly right. Fitted E_F: hankel=0.695eV, sqrtx=0.700eV,
  peak3-peak2/2(Lorentzian)=0.674eV -- cluster now 0.674-0.700eV. Re-saved all 4 E_F
  panels + CSV to the same folder.
- **Bug found and fixed in the Lorentzian fit itself**: user visually inspected a zoomed
  peak2 preview (`scratchpad/peak2_zoom_check.png`) and noticed raw data points visibly
  sticking out past the fitted curve. Root cause: the unbounded 4-parameter fit
  (B,A,x0,gamma) on a +-6-point window was under-constrained -- fitted gamma came out
  wildly inconsistent across wn (189.9nm at 860cm-1, but 4687nm/3496nm/6408nm at
  870/900/911cm-1), i.e. the optimizer found a nearly-flat degenerate solution at several
  wn instead of actually tracing the peak's real (fairly narrow, somewhat asymmetric)
  shape. Adding bounds alone just clamped gamma to the bound (still wrong, all wn pinned
  at the same 184.6nm ceiling). **Fixed** by narrowing the fit window to +-3 points (this
  peak shape isn't a sharp resonance, so a wider window invites the broad/degenerate
  solution) plus bounding gamma to a sane range relative to that smaller window span.
  Re-verified visually (`scratchpad/peak2_zoom_narrowwin.png`) before adopting -- gamma
  now 80-140nm across wn, curves visibly hug the actual peak data. Refitted E_F:
  hankel=0.695eV, sqrtx=0.700eV, peak3-peak2/2(Lorentzian, fixed)=0.688eV. Re-saved all 4
  E_F panels + CSV to the same folder (values shifted slightly from the buggy-gamma
  version: 860/870/880/890/900/911cm-1 q(peak2-3) now 1.998/2.107/2.341/2.541/2.743/3.203,
  vs 1.991/2.182/2.356/2.602/2.838/3.317 before the window/bounds fix).
- Added a real-space Lorentzian-fit preview to the script (`lorentzian_fit_realspace_preview_
  860to911.png`, same folder): 2x3 grid, one panel per wn, full bg-subtracted signal in gray
  with peak2/peak3's fitted Lorentzian curves overlaid (legend gives each fit's x0), and a
  text box giving both `2q` (the raw, pre-/2-correction value from the peak2-peak3 spacing)
  and `q` (`2q`/2) directly on each panel. 911cm-1's peak3 fit visibly diverges from a clean
  single-peak shape (broad/plateau-like region, x0 shifted further from the raw grid max
  than any other wn/peak) -- consistent with the low-SNR + possible q/2q beating explanation
  already logged, not a new issue.
