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
