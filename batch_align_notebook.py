import nbformat as nbf
import os

filepath = 'fitting_pipeline.ipynb'
try:
    with open(filepath, 'r', encoding='utf-8') as f:
        nb = nbf.read(f, as_version=4)
except Exception as e:
    print(f"Error loading nb: {e}")
    exit(1)

md_cell = nbf.v4.new_markdown_cell("""## 1. Batch Data Alignment & Initialization
Here we load all 15 wavenumbers, visualize them in an unaligned state, apply a customizable edge shift via `align_dict`, and visualize the aligned result. Once aligned, we select one target wavenumber to pass to the downstream fitting pipelines.""")

code_unaligned = nbf.v4.new_code_cell("""import glob
import re
import matplotlib.cm as cm

data_dir = 'data/graphene_3x1'
file_paths = glob.glob(f'{data_dir}/*_AVG_lp1.csv')

# Extract wavenumbers and sort
def get_wn(p):
    match = re.search(r'(\d+)cm-1', p)
    return int(match.group(1)) if match else 0

file_paths = sorted(file_paths, key=get_wn)
wn_list = [f"{get_wn(p)}cm-1" for p in file_paths]

# Load raw data
raw_data_dict = {}
for wn, path in zip(wn_list, file_paths):
    raw_data_dict[wn] = pd.read_csv(path)

# Plot 1x3 Unaligned Waterfall
fig, axs = plt.subplots(1, 3, figsize=(15, 6))
colors = cm.viridis(np.linspace(0, 1, len(wn_list)))

offset_A, offset_P, offset_Z = 0.5, 0.5, 2.0  # visual offsets

for i, wn in enumerate(wn_list):
    df = raw_data_dict[wn]
    x = df['distance_nm'] / 1000.0  # to um
    
    axs[0].plot(x, df['O3A'] - offset_A * i, color=colors[i])
    axs[1].plot(x, df['O3P'] - offset_P * i, color=colors[i])
    axs[2].plot(x, df['Z_nm_corrected'] - offset_Z * i, color=colors[i])

axs[0].set_title('O3A (Unaligned)')
axs[1].set_title('O3P (Unaligned)')
axs[2].set_title('Z_nm_corrected (Unaligned)')

for ax in axs:
    ax.set_xlabel('Raw Distance (μm)')
    ax.axvline(0.35, color='gray', linestyle='--', alpha=0.5)
plt.tight_layout()
""")

code_aligned = nbf.v4.new_code_cell("""# 🌟 Interactive Alignment Dictionary 🌟
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
    
    axs[0].plot(x, df['O3A'] - offset_A * i, color=colors[i])
    axs[1].plot(x, df['O3P'] - offset_P * i, color=colors[i])
    axs[2].plot(x, df['Z_nm_corrected'] - offset_Z * i, color=colors[i])

axs[0].set_title('O3A (Aligned)')
axs[1].set_title('O3P (Aligned)')
axs[2].set_title('Z_nm_corrected (Aligned)')

for ax in axs:
    ax.set_xlabel('Aligned Distance (μm)')
    ax.axvline(0.0, color='red', linestyle='--', alpha=0.8, lw=2) # The aligned edge

plt.tight_layout()
""")

code_select = nbf.v4.new_code_cell("""# === Select Target Wavenumber for Downstream Fitting ===
target_wn = '860cm-1'
L_cutoff = 0.9

df_target = aligned_data_dict[target_wn]
x_mat = df_target['distance_um'].values
y_mat = df_target['O3A'].values

# Phase data preparation (used by downstream FFT cell)
# Because some downstream scripts rely on pandas Series mapped directly,
# we construct amplp and phaselp dummy DataFrames to preserve compatibility
amplp = pd.DataFrame({'distance_um': x_mat, f'{target_wn}_O3A': y_mat})
phaselp = pd.DataFrame({f'{target_wn}_O3P': df_target['O3P'].values})

# Keep only positive x (on the material) up to L_cutoff for CHT/RealSpace fits
mask_fit = (x_mat >= 0) & (x_mat <= L_cutoff)
x_f = x_mat[mask_fit]
sig_f = y_mat[mask_fit]

print(f"Loaded {target_wn} for fitting. Edge aligned at x=0. Total fitting points: {len(x_f)}")
""")

new_cells = []
for idx, cell in enumerate(nb.cells):
    if idx == 2:
        new_cells.extend([md_cell, code_unaligned, code_aligned, code_select])
    elif idx == 3:
        pass # replace old cell 3
    else:
        new_cells.append(cell)

nb.cells = new_cells

with open(filepath, 'w', encoding='utf-8') as f:
    nbf.write(nb, f)

print("Notebook batched aligned structure successfully generated!")
