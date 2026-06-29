"""Folder loading for the web tool. Thin wrapper around loadnpz.py."""
import glob
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append('/Users/shizhe/envsetting')

from loadnpz import loadnpz  # noqa: E402


def get_wn(path_or_name):
    m = re.search(r'(\d+(?:\.\d+)?)\s*cm-1', path_or_name)
    return float(m.group(1)) if m else 0.0


def load_folder(data_dir):
    """Load every *_AVG.npz (or *.npz) file in data_dir. Returns (scans, wn_list)
    where wn_list is a list of wn label strings sorted high -> low."""
    paths = sorted(glob.glob(os.path.join(data_dir, '*.npz')))
    if not paths:
        raise FileNotFoundError(f"No .npz files found in {data_dir}")

    scans = {}
    for p in paths:
        wn_val = get_wn(os.path.basename(p))
        label = f"{wn_val:g}cm-1"
        scans[label] = loadnpz(p)

    wn_list = sorted(scans.keys(), key=get_wn, reverse=True)
    return scans, wn_list
