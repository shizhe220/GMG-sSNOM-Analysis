# GMG — Project Instructions

nano-FTIR (sSNOM) graphene plasmon dispersion analysis: manually-extracted linecuts at
multiple wavenumbers -> CHT/Hankel/1-sqrtx/FFT fitting -> q_p vs wn dispersion -> E_F.
Three datasets: `graphene_3x1_manual` (Gr/SiO2/Si), `graphene_4x1_manual` lp1="edge"/
lp2="interior" (MoO3(4.2nm)/Gr/SiO2/Si). See `data/`/`figures/` per-dataset subfolders.

**Read `log_and_conventions.md`'s Progress Log before starting substantial work** — it has
the full dated history (bug fixes, architecture decisions, why things are the way they are).
This file only holds the stable structure/rules; the log is where the "why" lives.

## Folder Structure (as of 2026-07-10)

```
GMG/
├── CLAUDE.md                          This file (stable rules) -- see log_and_conventions.md for history.
├── log_and_conventions.md             Progress Log (dated entries) -- read for "why", not loaded automatically.
├── nanoftir_Shizhe.py                 Thin re-export shim + GMG-specific glue (load_aligned_wn_signal,
│                                        run_wn_comparison). The actual CHT/Hankel/1-sqrtx/FFT engine now
│                                        lives in /Users/shizhe/envsetting/snippet/dispersion_fitting.py
│                                        (migrated 2026-07-08; see envsetting/dispersion_fitting_migration_plan.md).
├── loadnpz.py                         Reader for single-frequency npz (graphene_3x1/3x1_manual).
├── loadnpz_combined.py                Reader for combined multi-frequency npz (graphene_4x1_manual).
│
├── manual_linecut_pipeline.ipynb              Manual linecut extraction, graphene_3x1_manual.
├── manual_linecut_pipeline_graphene_4x1.ipynb Manual linecut extraction, graphene_4x1_manual (lp1+lp2 both edges).
├── fitting_pipeline_graphene_3x1_manual.ipynb CHT/Hankel/1-sqrtx/FFT fitting, per wn, graphene_3x1_manual.
├── fitting_pipeline_graphene_4x1_lp1.ipynb    Same, graphene_4x1_manual lp1 (edge).
├── fitting_pipeline_graphene_4x1_lp2.ipynb    Same, graphene_4x1_manual lp2 (interior).
├── rp_dispersion_plot.ipynb                   Im(r_p) dispersion maps vs extracted q_p, per dataset.
│
├── data/
│   ├── processed_3x1um/                       Raw npz per wn, graphene_3x1(_manual).
│   ├── processed_4x1um/860-1000_AVG.npz       Raw npz, ALL 15 wn in one combined file, graphene_4x1_manual.
│   ├── graphene_3x1_manual/                   Exported linecut CSVs + linecut_geometry.csv.
│   ├── graphene_4x1_manual/                   Exported linecut CSVs (lp1+lp2) + linecut_geometry.csv.
│   └── fit_results_graphene_*.pkl             Per-wn lambda/q/damping, one pkl per dataset.
│
├── figures/<dataset>/{cht,realspace,fft,fit_comparison}/   Per-dataset fit figures.
├── figures/rp_dispersion/                     Im(r_p) maps.
│
├── scripts/                                   Standalone regen scripts, chdir to repo root, run from anywhere.
│   ├── save_fit_results_graphene_*.py           notebook-tuned params -> pkl + cht/realspace/fft figures
│   ├── save_fit_comparison_figures_graphene_*.py -> fit_comparison/ combined+panels figures
│   ├── export_comparison_csv_graphene_*.py      pkl -> comparison CSV
│   └── make_overview_pptx*.py                   figures -> slides/*.pptx (font Aptos)
│
├── docs/            Reference papers (Woessner et al. Nat. Mater. 2015 + SI).
├── slides/          Presentation decks (.pptx).
└── webtool/         Streamlit v1 interactive tool -- known rough edges, not the primary workflow.
```

**Shared code outside this repo**: `/Users/shizhe/envsetting/snippet/` (`dispersion_fitting.py`,
`linecut_extraction.py`, `extract_linecut.py`) -- see `/Users/shizhe/envsetting/PROJECTS.md` and
the global `~/.claude/CLAUDE.md` for the cross-project index/conventions.

## Conventions & Rules

### 1. Extensive Logging Policy
Every non-trivial change (bug fix, architecture decision, parameter re-tune) gets a dated
entry in `log_and_conventions.md`'s Progress Log — what changed, why, and where. This is the
project's institutional memory; don't let a decision live only in conversation history.

### 2. Code modification policy: non-destructive
New functions get appended; existing function signatures/behavior don't change without a
clear reason logged. When extending a widely-called shared function (e.g.
`compare_cavity_models`), grep every call site first to confirm a default-value change is
actually safe.

### 3. Plotting Style (APS Standard)
```python
plt.rcParams.update({
    'font.size': 12, 'font.family': 'Arial',
    'axes.linewidth': 1.5, 'xtick.major.width': 1.5, 'ytick.major.width': 1.5,
    'xtick.direction': 'in', 'ytick.direction': 'in', 'xtick.top': True,
    'ytick.right': True
})
```
- **Labels**: `fontweight='bold'` for `xlabel`/`ylabel`.
- **Legends**: `frameon=False`, `fontsize=12` (or smaller if the legend carries per-curve
  stats text, e.g. `fontsize=9`).
- **Math**: LaTeX in strings, e.g. `r'Frequency (cm$^{-1}$)'`.
- **Panel labels**: top-left inside the panel, `(a)`, `(b)`, `(c)`, bold, `fontsize=14`.
- **Fixed colors per method** (reuse everywhere, don't invent new ones): CHT `#b2182b`,
  Hankel `#1c7293`, 1/√x `#e08214`.
- **PPTX decks**: font Aptos.
- Before delivering a publication-facing figure, run it through
  `/Users/shizhe/envsetting/snippet/FIGURE_QA_CHECKLIST.md`.

### 4. Self-check discipline
After any code/data change, verify numerically — back up the pre-change output, re-run,
diff every value programmatically, report the diff count/magnitude (not just "ran without
error"). For PPTX/figures with no renderer available in this environment, verify via
python-pptx shape-bounds introspection instead of eyeballing.
