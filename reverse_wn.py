import nbformat as nbf

filepath = 'fitting_pipeline.ipynb'
with open(filepath, 'r', encoding='utf-8') as f:
    nb = nbf.read(f, as_version=4)

count = 0
for cell in nb.cells:
    if cell.cell_type == 'code':
        if 'file_paths = sorted(file_paths, key=get_wn)' in cell.source:
            cell.source = cell.source.replace('file_paths = sorted(file_paths, key=get_wn)', 'file_paths = sorted(file_paths, key=get_wn, reverse=True)')
            count += 1
            
with open(filepath, 'w', encoding='utf-8') as f:
    nbf.write(nb, f)

print(f"Reversed wn sorting in {count} cells.")
