# GMG Linecut & Fit Tool (v1)

Interactive replacement for the manual notebook workflow (`manual_linecut_pipeline.ipynb`
-> `fitting_pipeline.ipynb`): load -> linecut -> align -> fit/FFT -> export, with
mouse-click interaction instead of typing every coordinate.

## Run

```bash
pip install -r webtool/requirements.txt
streamlit run webtool/app.py
```

## Workflow

1. **Mapping & Linecut** -- pick a wavenumber, view its amp/phase/Z map (Sky/bwr_r/gray,
   same convention as `snippet.linecut_extraction`). Click twice on the map to set the
   cut's start/end (or type the numbers directly); optionally flatten Z first
   (auto-suggest 3 flat points, or click your own on the Z map). "Extract & add to
   waterfall" stores the linecut and updates the running waterfall preview.
2. **Align** -- pick a wavenumber, click on its highlighted curve in the waterfall at
   the point that should become distance=0 (or type a number). Shifts are cumulative
   on top of whatever the align value already was.
3. **Fit / FFT** -- per wavenumber, three tabs (CHT / real-space Hankel+1/sqrtx / FFT)
   reusing `nanoftir_Shizhe.py`'s existing fitting+plotting functions. Click the
   |T(k)| or FFT-spectrum preview to set the momentum guess instead of typing it.
4. **Export & Summary** -- write linecut CSVs (same column format as the old data,
   `distance_nm, O3A, O3P, Z_nm, Z_nm_corrected`) plus a fit-results pickle/CSV, and
   plot q_p vs wavenumber across all four methods.

Every step reads/writes the same `st.session_state`, so you can jump back to any
wavenumber at any step (sidebar shows a per-wn progress status: pending / cut /
aligned / fitted) without losing what you've already done.

## Status / known limits (v1)

- Click-to-draw is two clicks (start, end), not a literal click-and-drag line --
  Streamlit's plotly event capture doesn't cleanly support drag gestures. Click 1 sets
  the start point, click 2 sets the end point and computes angle/radius from them.
- Z-plane correction is a single global correction per wn (not re-computed per linecut
  geometry change) -- click points on the Z map, or use the auto-suggest button.
- The interactive bits (click handlers) were tested by exercising the underlying
  functions directly with real data, not by an actual mouse click in a browser (no
  headless browser available in the dev environment) -- please sanity-check the click
  interactions yourself on first use and flag anything that feels off.
- Real-space fit (Hankel/1-sqrtx) doesn't have a clickable initial guess yet -- it
  still relies on `lam0_guess` + automatic peak-finding, same as the notebook.
