import nbformat as nbf
import os

filepath = 'fitting_pipeline.ipynb'
with open(filepath, 'r', encoding='utf-8') as f:
    nb = nbf.read(f, as_version=4)

code_unaligned = """import glob
import re
import matplotlib.cm as cm
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

data_dir = 'data/graphene_3x1'
file_paths = glob.glob(f'{data_dir}/*_AVG_lp1.csv')

# Visualization Settings
normalize_A = True   # Normalize Amplitude
normalize_P = False  # Normalize Phase

# Extract wavenumbers and sort
def get_wn(p):
    match = re.search(r'(\d+)cm-1', p)
    return int(match.group(1)) if match else 0

file_paths = sorted(file_paths, key=get_wn)
wn_list = [f"{get_wn(p)}cm-1" for p in file_paths]

# Load raw data and apply custom Z correction
raw_data_dict = {}
for wn, path in zip(wn_list, file_paths):
    df = pd.read_csv(path)
    # The 'Z_nm_corrected' column has the linear slope removed.
    # We subtract its minimum so the substrate sits nicely at 0.
    df['Z_nm_custom'] = df['Z_nm_corrected'] - df['Z_nm_corrected'].min()
    raw_data_dict[wn] = df

# Plot 1x3 Unaligned Waterfall
fig, axs = plt.subplots(1, 3, figsize=(15, 10))
# Using a nice coolwarm colormap similar to the reference image
colors = cm.coolwarm(np.linspace(0, 1, len(wn_list)))

# Calculate clean offsets so curves don't overlap too much
max_A = max([raw_data_dict[wn]['O3A'].max() for wn in wn_list])
max_P = max([raw_data_dict[wn]['O3P'].max() for wn in wn_list])
max_Z = max([raw_data_dict[wn]['Z_nm_custom'].max() for wn in wn_list])

offset_A = 1.2 if normalize_A else max_A * 0.8
offset_P = 1.2 if normalize_P else max_P * 0.8
offset_Z = max_Z * 1.1

for i, wn in enumerate(wn_list):
    df = raw_data_dict[wn]
    x = df['distance_nm'] / 1000.0  # to um
    
    a_data = df['O3A']
    p_data = df['O3P']
    
    if normalize_A: a_data = a_data / a_data.max()
    if normalize_P: p_data = p_data / p_data.max()
    
    axs[0].plot(x, a_data - offset_A * i, color=colors[i], lw=1.5)
    axs[1].plot(x, p_data - offset_P * i, color=colors[i], lw=1.5)
    axs[2].plot(x, df['Z_nm_custom'] - offset_Z * i, color=colors[i], lw=1.5)

axs[0].set_title('O3A (Unaligned' + (', Norm)' if normalize_A else ')'))
axs[1].set_title('O3P (Unaligned' + (', Norm)' if normalize_P else ')'))
axs[2].set_title('Z_nm (Levelled)')

for ax in axs:
    ax.set_xlabel('Raw Distance (μm)')
    ax.set_yticks([]) # Hide Y-ticks for a much cleaner waterfall look
    ax.axvline(0.35, color='gray', linestyle='--', alpha=0.5)
plt.tight_layout()
"""

code_aligned = """# 🌟 Interactive Alignment Dictionary 🌟
# Modify the shift value (in nm) for any wavenumber to manually align its edge to X=0.
align_dict = {wn: 350.0 for wn in wn_list}

# Apply alignment
aligned_data_dict = {}
for wn in wn_list:
    df_aligned = raw_data_dict[wn].copy()
    shift_nm = align_dict[wn]
    # Shift X-axis so edge is exactly at 0
    df_aligned['distance_um'] = (df_aligned['distance_nm'] - shift_nm) / 1000.0
    aligned_data_dict[wn] = df_aligned

# Plot 1x3 Aligned Waterfall
fig, axs = plt.subplots(1, 3, figsize=(15, 10))

for i, wn in enumerate(wn_list):
    df = aligned_data_dict[wn]
    x = df['distance_um']
    
    a_data = df['O3A']
    p_data = df['O3P']
    
    if normalize_A: a_data = a_data / a_data.max()
    if normalize_P: p_data = p_data / p_data.max()
    
    axs[0].plot(x, a_data - offset_A * i, color=colors[i], lw=1.5)
    axs[1].plot(x, p_data - offset_P * i, color=colors[i], lw=1.5)
    axs[2].plot(x, df['Z_nm_custom'] - offset_Z * i, color=colors[i], lw=1.5)

axs[0].set_title('O3A (Aligned' + (', Norm)' if normalize_A else ')'))
axs[1].set_title('O3P (Aligned' + (', Norm)' if normalize_P else ')'))
axs[2].set_title('Z_nm (Levelled & Aligned)')

for ax in axs:
    ax.set_xlabel('Aligned Distance (μm)')
    ax.set_yticks([]) # Clean aesthetics
    ax.axvline(0.0, color='red', linestyle='--', alpha=0.8, lw=2) # The aligned edge

plt.tight_layout()
"""

count = 0
for cell in nb.cells:
    if cell.cell_type == 'code':
        if 'Plot 1x3 Unaligned Waterfall' in cell.source:
            cell.source = code_unaligned
            count += 1
        elif 'Plot 1x3 Aligned Waterfall' in cell.source:
            cell.source = code_aligned
            count += 1

with open(filepath, 'w', encoding='utf-8') as f:
    nbf.write(nb, f)

print(f"Notebook successfully updated. Replaced {count} cells.")
