import nbformat as nbf

filepath = 'fitting_pipeline.ipynb'
with open(filepath, 'r', encoding='utf-8') as f:
    nb = nbf.read(f, as_version=4)

for cell in nb.cells:
    if cell.cell_type == 'code' and '# === Select Target Wavenumber' in cell.source:
        if '# Background Subtraction' not in cell.source:
            lines = cell.source.split('\n')
            new_lines = []
            for line in lines:
                new_lines.append(line)
                if "y_mat = df_target['O3A'].values" in line:
                    new_lines.append("")
                    new_lines.append("from scipy.signal import savgol_filter")
                    new_lines.append("# Background Subtraction")
                    new_lines.append("window_len = min(41, len(y_mat) if len(y_mat) % 2 != 0 else len(y_mat)-1)")
                    new_lines.append("y_bg = savgol_filter(y_mat, window_length=window_len, polyorder=2)")
                    new_lines.append("y_mat = y_mat - y_bg")
            cell.source = '\n'.join(new_lines)

with open(filepath, 'w', encoding='utf-8') as f:
    nbf.write(nb, f)

print("Notebook patched to restore background subtraction!")
