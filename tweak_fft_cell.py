import nbformat as nbf

filepath = 'fitting_pipeline.ipynb'
with open(filepath, 'r', encoding='utf-8') as f:
    nb = nbf.read(f, as_version=4)

code_fft = """# === Standard 1D FFT Analysis ===
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import importlib
import nanoftir_Shizhe as nanof
importlib.reload(nanof)

try:
    print("Running standard FFT comparison...")
    x_dist = amplp['distance_um']
    amp_raw = amplp[f'{target_wn}_O3A']
    phase_raw = phaselp[f'{target_wn}_O3P']
    
    # We subtract background here too to ensure FFT is analyzing pure oscillations
    from scipy.signal import savgol_filter
    window_len = min(41, len(amp_raw) if len(amp_raw) % 2 != 0 else len(amp_raw)-1)
    amp_bg = savgol_filter(amp_raw, window_length=window_len, polyorder=2)
    amp_osc = amp_raw - amp_bg
    
    # 🌟 q_guess is exposed here! Adjust [5, 30] to search for peaks in different k regions.
    nanof.plot_channel_fft(
        x_dist, amp_osc, phase_raw, 
        label='O3', wn=target_wn,
        xr=(0.05, 0.9), q_range=(0, 25),
        window='hann', pad_factor=3.0,
        q_guess=[3] 
    )
except Exception as e:
    import traceback
    print(f"Error running FFT: {e}")
    traceback.print_exc()
"""

count = 0
for cell in nb.cells:
    if cell.cell_type == 'code' and '=== Standard 1D FFT Analysis ===' in cell.source:
        cell.source = code_fft
        count += 1

with open(filepath, 'w', encoding='utf-8') as f:
    nbf.write(nb, f)

print(f"FFT cell successfully updated. Replaced {count} cells.")
