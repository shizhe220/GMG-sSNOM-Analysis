"""
Extract per-wn Hankel-fit and FFT figures for graphene_4x1_manual lp2 (interior) from
fitting_pipeline_graphene_4x1_lp2.ipynb's already-executed cell outputs, for the slide
deck built by make_hankel_fft_pptx_lp2.py. Doesn't recompute anything -- pulls the exact
images already produced by the notebook's live-tuned params (same source of truth as
scripts/export_shared_csv_graphene_4x1_lp2.py), matching how lp1's own
figures/graphene_4x1_manual/lp1/two_wave_hankel/ images were sourced.

860-890cm-1 (single-wave regime): the compare_cavity_models "sorted by AIC" real-space
comparison plot (Hankel + 1/sqrt(x) overlay) -> {wn}cm-1_singlewave.png
900-1000cm-1 (two-wave regime): the single-wave-vs-two-wave Hankel comparison plot
-> {wn}cm-1_hankel2wave.png
All 15 wn: the 4-panel spatial-data + amp/phase/complex FFT figure -> {wn}cm-1_fft.png

Cell-index map (fixed structure of the notebook, checked 2026-07-23):
  860/870/880/890 -> compare_cavity_models cell, FFT cell (9/11, 19/21, 27/29, 35/37)
  900-1000        -> FFT cell, two-wave cell (45/47, 55/57, ..., 145/147)
"""
import os
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import json
import base64

NB_PATH = 'fitting_pipeline_graphene_4x1_lp2.ipynb'
OUT_DIR = 'figures/graphene_4x1_manual/lp2/two_wave_hankel'
os.makedirs(OUT_DIR, exist_ok=True)

SINGLEWAVE_CELLS = {860: (9, 11), 870: (19, 21), 880: (27, 29), 890: (35, 37)}
TWOWAVE_CELLS = {
    900: (45, 47), 911: (55, 57), 920: (65, 67), 930: (75, 77), 941: (85, 87),
    950: (95, 97), 960: (105, 107), 970: (115, 117), 980: (125, 127),
    991: (135, 137), 1000: (145, 147),
}

nb = json.load(open(NB_PATH))


def last_png(cell_idx):
    outputs = nb['cells'][cell_idx]['outputs']
    pngs = [o['data']['image/png'] for o in outputs if 'image/png' in o.get('data', {})]
    if not pngs:
        raise RuntimeError(f'cell {cell_idx} has no image/png output')
    return base64.b64decode(pngs[-1])


for wn, (fit_cell, fft_cell) in SINGLEWAVE_CELLS.items():
    open(f'{OUT_DIR}/{wn}cm-1_singlewave.png', 'wb').write(last_png(fit_cell))
    open(f'{OUT_DIR}/{wn}cm-1_fft.png', 'wb').write(last_png(fft_cell))
    print(f'{wn}cm-1: singlewave (cell {fit_cell}) + fft (cell {fft_cell})')

for wn, (fft_cell, fit_cell) in TWOWAVE_CELLS.items():
    open(f'{OUT_DIR}/{wn}cm-1_hankel2wave.png', 'wb').write(last_png(fit_cell))
    open(f'{OUT_DIR}/{wn}cm-1_fft.png', 'wb').write(last_png(fft_cell))
    print(f'{wn}cm-1: hankel2wave (cell {fit_cell}) + fft (cell {fft_cell})')

print(f'\nSaved all images to {OUT_DIR}/')
