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

# Visualization Setting
normalize_max = True  # Set to False to view absolute amplitudes in the waterfall

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
    # Custom Z correction: subtract min to ensure Z >= 0
    df['Z_nm_custom'] = df['Z_nm'] - df['Z_nm'].min()
    raw_data_dict[wn] = df

# Plot 1x3 Unaligned Waterfall
fig, axs = plt.subplots(1, 3, figsize=(15, 6))
colors = cm.viridis(np.linspace(0, 1, len(wn_list)))

# Offsets for waterfall plotting
offset_A = 0.5 if normalize_max else (raw_data_dict[wn_list[0]]['O3A'].max() * 0.5)
offset_P = 0.5 if normalize_max else (raw_data_dict[wn_list[0]]['O3P'].max() * 0.5)
offset_Z = 2.0  

for i, wn in enumerate(wn_list):
    df = raw_data_dict[wn]
    x = df['distance_nm'] / 1000.0  # to um
    
    a_data = df['O3A']
    p_data = df['O3P']
    
    if normalize_max:
        a_data = a_data / a_data.max()
        p_data = p_data / p_data.max()
    
    axs[0].plot(x, a_data - offset_A * i, color=colors[i])
    axs[1].plot(x, p_data - offset_P * i, color=colors[i])
    axs[2].plot(x, df['Z_nm_custom'] - offset_Z * i, color=colors[i])

axs[0].set_title('O3A (Unaligned' + (', Normalized)' if normalize_max else ')'))
axs[1].set_title('O3P (Unaligned' + (', Normalized)' if normalize_max else ')'))
axs[2].set_title('Z_nm_custom (Unaligned)')

for ax in axs:
    ax.set_xlabel('Raw Distance (μm)')
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
fig, axs = plt.subplots(1, 3, figsize=(15, 6))

for i, wn in enumerate(wn_list):
    df = aligned_data_dict[wn]
    x = df['distance_um']
    
    a_data = df['O3A']
    p_data = df['O3P']
    
    if normalize_max:
        a_data = a_data / a_data.max()
        p_data = p_data / p_data.max()
    
    axs[0].plot(x, a_data - offset_A * i, color=colors[i])
    axs[1].plot(x, p_data - offset_P * i, color=colors[i])
    axs[2].plot(x, df['Z_nm_custom'] - offset_Z * i, color=colors[i])

axs[0].set_title('O3A (Aligned' + (', Normalized)' if normalize_max else ')'))
axs[1].set_title('O3P (Aligned' + (', Normalized)' if normalize_max else ')'))
axs[2].set_title('Z_nm_custom (Aligned)')

for ax in axs:
    ax.set_xlabel('Aligned Distance (μm)')
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
