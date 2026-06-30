import io
import json
import re
from pathlib import Path

import numpy as np


def _decode_npz_text(arr):
    if hasattr(arr, "tobytes"):
        return arr.tobytes().decode("utf-8")
    return bytes(arr).decode("utf-8")


def _lp_sort_key(key):
    m = re.fullmatch(r"lp(\d+)", key)
    return int(m.group(1)) if m else 10**9


_ROI_COLORS = ["#ff5c8a", "#5cc8ff", "#8aff5c", "#ffd05c", "#c08aff"]
_LP_COLORS = ["#5cc8ff", "#ff5c8a", "#8aff5c", "#ffd05c", "#c08aff"]


class Scan(dict):
    """A loaded FwdBwd .npz summary: a dict with built-in plotting.

    Dict access is unchanged / backward compatible:
        scan["O3A"]                 # 2D channel array, shape (ny, nx)
        scan["info"]                # parsed info.json
        scan["stats"]               # ROI stats DataFrame (or None)
        scan["point_spectra"]       # ROI ratio spectra DataFrame (or None)
        scan["lineprofile"]["lp1"]  # line-profile DataFrame (or None)
        scan["pixelsize_nm"]

    Plotting helpers need no fig/ax setup and never ask you for pixels — they
    use the coordinates stored in the file (which are this image's own native
    pixels), so the boxes/lines land in the right place automatically even when
    images have different pixel sizes:
        scan.show()              # image + crop box + every ROI box + every line
        scan.show_roi("P1")      # image with ROI box(es) overlaid
        scan.show_lp("L1")       # image with line-profile path(s) overlaid
        scan.plot_lp("L1")       # the extracted profile curves vs distance (nm)
    Omit the name to draw all ROIs / all lines.

    Combined multi-frequency files (saved when several frequencies are loaded;
    info["combined"] is True) instead expose:
        scan["frequencies"]                     # list of frequencies (cm^-1)
        scan["maps"][freq][channel]             # 2D averaged map per frequency
        scan["lineprofile"]["lp1"]["O3A"]        # DataFrame: distance_nm, <freq1>, <freq2>, ...
        scan["lineprofile"]["lp1"]["O3P"]        # same distance_nm grid & freq columns
        scan["lineprofile"]["lp1"]["Z_nm"]       # corrected topography, same grid/columns
        scan["fft"]["fft1"]["O3A"]               # DataFrame: q_cm1, k_um1, <freq1>, ...
        scan["fft"]["fft1"]["O3P"]               # same q grid & freq columns
        scan["fft"]["fft1"]["O3complex"]         # FFT of amp*e^(i*phase), same grid/columns
        scan["point_spectra"]                    # ROI ratios, one row per frequency
        scan.show_map(freq, channel)             # one frequency's map
        scan.plot_lp_freq("O3A")                 # a channel vs distance, all frequencies
        scan.plot_fft_freq("O3A")                # a channel's FFT vs q, all frequencies
    """

    # ----- combined multi-frequency files -------------------------------
    def is_combined(self):
        return bool(self.get("combined"))

    def freqs(self):
        return list(self.get("frequencies") or sorted((self.get("maps") or {}).keys()))

    def get_map(self, freq, channel=None):
        maps = self.get("maps") or {}
        if freq not in maps:
            for k in maps:
                if str(k) == str(freq):
                    freq = k
                    break
        chans = maps.get(freq)
        if not chans:
            raise KeyError("frequency %r not found; have %s" % (freq, self.freqs()))
        if channel is None:
            for pref in ("O3A", "O2A", "O4A", "Z"):
                if pref in chans:
                    channel = pref
                    break
            else:
                channel = next(iter(chans))
        if channel not in chans:
            raise KeyError("channel %r not at %r; have %s" % (channel, freq, list(chans)))
        return channel, np.asarray(chans[channel])

    def show_map(self, freq=None, channel=None, ax=None, **imshow_kw):
        import matplotlib.pyplot as plt
        if freq is None:
            freq = self.freqs()[0]
        ch, img = self.get_map(freq, channel)
        if ax is None:
            _, ax = plt.subplots()
        cmap = imshow_kw.pop("cmap", "gray" if ch.upper().startswith("Z") else "viridis")
        ax.imshow(img, origin="lower", cmap=cmap, **imshow_kw)
        ax.set_title("%s  %s" % (freq, ch), fontsize=9)
        return ax

    def _resolve_key(self, container, name):
        keys = list(container.keys())
        if name is None:
            return keys[0]
        s = str(name).strip().lower()
        digits = re.sub(r"\D", "", s)
        for k in keys:
            if s == k.lower() or (digits and re.sub(r"\D", "", k) == digits):
                return k
        return keys[0]

    def _resolve_role(self, chans, name):
        """Pick a channel sub-table (e.g. "O3A"/"O3P"/"Z_nm") from a {role: table} dict."""
        keys = list(chans.keys())
        if name is None:
            return keys[0]
        s = str(name).strip().lower()
        for k in keys:
            if s == k.lower():
                return k
        for k in keys:
            if k.lower().startswith(s) or s.startswith(k.lower()):
                return k
        return keys[0]

    def _freq_cols(self, table, exclude):
        """cols, get(col)->ndarray, and the subset of columns that are frequency values."""
        cols = get = None
        try:
            import pandas as pd
            if isinstance(table, pd.DataFrame):
                cols = list(table.columns)
                get = lambda c: table[c].to_numpy()
        except Exception:
            pass
        if cols is None:
            arr = np.asarray(table)
            cols = list(range(arr.shape[1]))
            get = lambda c: arr[:, cols.index(c)]
        freq_cols = [c for c in cols if str(c) not in exclude]
        return cols, get, freq_cols

    def plot_lp_freq(self, channel="O3A", name=None, ax=None):
        """Plot one line-profile channel (default O3A) vs distance for every frequency."""
        import matplotlib.pyplot as plt
        lps = self.get("lineprofile") or {}
        if not lps:
            raise ValueError("this file has no line profiles")
        key = self._resolve_key(lps, name)
        chans = lps[key]
        if not isinstance(chans, dict):
            raise ValueError("%r is not a combined multi-frequency line profile" % key)
        role = self._resolve_role(chans, channel)
        cols, get, freq_cols = self._freq_cols(chans[role], ("distance_nm",))
        dist = get("distance_nm")
        if ax is None:
            _, ax = plt.subplots(figsize=(6, 3.2))
        for fl in freq_cols:
            ax.plot(dist, get(fl), lw=1, label=str(fl))
        ax.set_xlabel("distance (nm)"); ax.set_ylabel(role)
        ax.legend(fontsize=7, ncol=2, title="cm$^{-1}$")
        ax.set_title("%s  %s vs frequency" % (key, role), fontsize=9)
        return ax

    def plot_fft_freq(self, channel="O3A", name=None, ax=None):
        """Plot one FFT channel (default O3A) vs q for every frequency."""
        import matplotlib.pyplot as plt
        ffts = self.get("fft") or {}
        if not ffts:
            raise ValueError("this file has no FFT tables")
        key = self._resolve_key(ffts, name)
        chans = ffts[key]
        if not isinstance(chans, dict):
            raise ValueError("%r is not a combined multi-frequency FFT table" % key)
        role = self._resolve_role(chans, channel)
        cols, get, freq_cols = self._freq_cols(chans[role], ("q_cm1", "k_um1"))
        q = get("q_cm1")
        if ax is None:
            _, ax = plt.subplots(figsize=(6, 3.2))
        for fl in freq_cols:
            ax.plot(q, get(fl), lw=1, label=str(fl))
        ax.set_xlabel("q (cm$^{-1}$)"); ax.set_ylabel("%s FFT (norm.)" % role)
        ax.legend(fontsize=7, ncol=2, title="cm$^{-1}$")
        ax.set_title("%s  %s FFT vs frequency" % (key, role), fontsize=9)
        return ax

    # ----- channels -----------------------------------------------------
    def channel_names(self):
        return [k for k, v in self.items()
                if isinstance(v, np.ndarray) and v.ndim == 2]

    def default_channel(self):
        names = self.channel_names()
        for pref in ("O3A", "O2A", "O4A", "O1A", "O5A"):
            if pref in names:
                return pref
        for n in names:
            if re.fullmatch(r"O\d+A", n):
                return n
        for pref in ("Z", "Z C"):
            if pref in names:
                return pref
        return names[0] if names else None

    def image(self, channel=None):
        if self.is_combined():
            raise TypeError("combined multi-frequency file — use show_map(freq, channel), "
                            "plot_lp_freq(channel), plot_fft_freq(channel), or scan['maps'][freq][channel]")
        ch = channel or self.default_channel()
        if ch is None or not isinstance(self.get(ch), np.ndarray):
            raise KeyError("no 2D channel found; available: %s" % self.channel_names())
        return ch, np.asarray(self[ch])

    # ----- coordinate lookups (from info.json) --------------------------
    def roi_points(self):
        return ((self.get("info") or {}).get("roi") or {}).get("points") or []

    def line_meta(self):
        return ((self.get("info") or {}).get("lineprofile") or {}).get("profiles") or []

    def _pick_roi(self, name):
        pts = self.roi_points()
        if not pts:
            raise ValueError("this file has no ROI points")
        if name is None:
            return list(enumerate(pts))
        key = str(name).strip().lower()
        for i, p in enumerate(pts):
            label = str(p.get("point", "P%d" % (i + 1))).lower()
            if key in (label, label.lstrip("p"), "p%d" % (i + 1), str(i + 1)):
                return [(i, p)]
        raise KeyError("ROI %r not found; have %s"
                       % (name, [p.get("point") for p in pts]))

    def _pick_lp(self, name):
        meta = self.line_meta()
        if not meta:
            raise ValueError("this file has no line profiles")
        if name is None:
            return list(enumerate(meta))
        key = str(name).strip().lower()
        for i, m in enumerate(meta):
            lpkey = str(m.get("key", "lp%d" % (i + 1))).lower()   # e.g. 'lp1'
            num = lpkey.replace("lp", "")
            if key in (lpkey, "lp" + num, "l" + num, num, "l%d" % (i + 1), str(i + 1)):
                return [(i, m)]
        raise KeyError("line %r not found; have %s"
                       % (name, [m.get("key") for m in meta]))

    # ----- plotting -----------------------------------------------------
    def show_image(self, channel=None, ax=None, **imshow_kw):
        import matplotlib.pyplot as plt
        ch, img = self.image(channel)
        if ax is None:
            ny, nx = img.shape
            _, ax = plt.subplots(figsize=(5.0, 5.0 * ny / max(1, nx)))
        topo = ch in ("Z", "Z C", "Z OFFSET")
        cmap = imshow_kw.pop("cmap", "gray" if topo else "inferno")
        ax.imshow(img, origin="upper", cmap=cmap, **imshow_kw)
        prefix = (self.get("info") or {}).get("dataset_prefix", "")
        ax.set_title(("%s\n%s" % (prefix, ch)).strip(), fontsize=8)
        ax.set_xlabel("x (px)")
        ax.set_ylabel("y (px)")
        return ax

    def _draw_crop(self, ax):
        from matplotlib.patches import Rectangle
        r = ((self.get("info") or {}).get("crop") or {}).get("rect_native_px")
        if not r:
            return
        ax.add_patch(Rectangle((r["x0"], r["y0"]), r["x1"] - r["x0"], r["y1"] - r["y0"],
                               fill=False, edgecolor="#ff3048", lw=1.4))
        ax.text(r["x0"] + 3, r["y0"] + 12, "crop", color="#ff3048", fontsize=7)

    def show_roi(self, name=None, channel=None, ax=None, crop=True):
        from matplotlib.patches import Rectangle
        ax = self.show_image(channel, ax)
        if crop:
            self._draw_crop(ax)
        roi = (self.get("info") or {}).get("roi") or {}
        box = int(roi.get("roi_px") or (2 * int(roi.get("roi_half_px", 5)) + 1))
        half = box / 2.0
        for i, p in self._pick_roi(name):
            x, y = p["x_px"], p["y_px"]
            c = _ROI_COLORS[i % len(_ROI_COLORS)]
            ax.add_patch(Rectangle((x - half, y - half), box, box,
                                   fill=False, edgecolor=c, lw=1.6))
            ax.plot([x], [y], "+", color=c, ms=6)
            ax.text(x + half + 1, y - half, p.get("point", "P%d" % (i + 1)),
                    color=c, fontsize=8)
        return ax

    def show_lp(self, name=None, channel=None, ax=None, crop=True):
        ax = self.show_image(channel, ax)
        if crop:
            self._draw_crop(ax)
        for i, m in self._pick_lp(name):
            x0, y0 = m["native_start_px"]
            x1, y1 = m["native_end_px"]
            c = _LP_COLORS[i % len(_LP_COLORS)]
            ax.plot([x0, x1], [y0, y1], "-", color=c, lw=1.8)
            ax.plot([x0, x1], [y0, y1], ".", color=c, ms=3)
            w = float(m.get("width_px", 1) or 1)
            if w > 1:
                dx, dy = x1 - x0, y1 - y0
                length = (dx * dx + dy * dy) ** 0.5 or 1.0
                ux, uy = -dy / length, dx / length
                for s in (w / 2.0, -w / 2.0):
                    ax.plot([x0 + ux * s, x1 + ux * s], [y0 + uy * s, y1 + uy * s],
                            "--", color=c, lw=0.7, alpha=0.6)
            label = "L%s" % str(m.get("key", "lp%d" % (i + 1))).replace("lp", "")
            ax.text(x0 + 3, y0 - 5, label, color=c, fontsize=8)
        return ax

    def show(self, channel=None, ax=None):
        ax = self.show_image(channel, ax)
        self._draw_crop(ax)
        if self.roi_points():
            self.show_roi(None, channel=channel, ax=ax, crop=False)
        if self.line_meta():
            self.show_lp(None, channel=channel, ax=ax, crop=False)
        return ax

    def _lp_arrays(self, m):
        key = m.get("key")
        data = (self.get("lineprofile") or {}).get(key)
        if data is None:
            raise KeyError("no profile data stored for %s" % key)
        try:
            import pandas as pd
            if isinstance(data, pd.DataFrame):
                dist = data["distance_nm"].to_numpy()
                cols = {c: data[c].to_numpy() for c in data.columns if c != "distance_nm"}
                return dist, cols
        except Exception:
            pass
        arr = np.asarray(data)
        names = m.get("columns") or (["distance_nm"]
                                     + ["c%d" % j for j in range(arr.shape[1] - 1)])
        dist = arr[:, 0]
        cols = {names[j]: arr[:, j] for j in range(1, arr.shape[1])}
        return dist, cols

    def plot_lp(self, name=None, axes=None):
        if self.is_combined():
            return self.plot_lp_freq(name=name)
        import matplotlib.pyplot as plt
        i, m = self._pick_lp(name)[0]
        dist, cols = self._lp_arrays(m)
        amp = next((k for k in cols if k.endswith("A")), None)
        ph = next((k for k in cols if k.endswith("P")), None)
        # prefer the step-aware slope-corrected Z (matches the GUI), else raw Z
        topo = (next((k for k in cols if "Z" in k.upper() and "CORRECT" in k.upper()), None)
                or next((k for k in cols if "Z" in k.upper()), None))
        if axes is None:
            _, axes = plt.subplots(2, 1, figsize=(6, 4), sharex=True,
                                   gridspec_kw={"height_ratios": [1.6, 1.0]})
        ax_top, ax_bot = axes
        if amp:
            ax_top.plot(dist, cols[amp], color="#1f77b4")
            ax_top.set_ylabel(amp, color="#1f77b4")
        if ph:
            axp = ax_top.twinx()
            axp.plot(dist, cols[ph], color="#e377c2")
            axp.set_ylabel(ph, color="#e377c2")
        if topo:
            ax_bot.plot(dist, cols[topo], color="#ff7f0e")
        ax_bot.set_ylabel((topo or "Z") + " (nm)", color="#ff7f0e")
        ax_bot.set_xlabel("distance (nm)")
        ax_top.set_title("L%s profile" % str(m.get("key", "lp%d" % (i + 1))).replace("lp", ""),
                         fontsize=9)
        return axes


def _load_combined(z, info, as_df):
    """Load a combined multi-frequency npz (info['combined'] is True).

    Layout:
      scan["maps"][freq][channel]         -> 2D averaged map (ny, nx)
      scan["frequencies"]                 -> list of frequencies (cm^-1)
      scan["lineprofile"]["lp1"]["O3A"]   -> DataFrame: distance_nm, <freq1>, <freq2>, ...
      scan["lineprofile"]["lp1"]["O3P"]   -> same distance_nm grid & freq columns
      scan["lineprofile"]["lp1"]["Z_nm"]  -> corrected topography, same grid/columns
      scan["fft"]["fft1"]["O3A"]          -> DataFrame: q_cm1, k_um1, <freq1>, ...
      scan["fft"]["fft1"]["O3P"]          -> same q grid & freq columns
      scan["fft"]["fft1"]["O3complex"]    -> FFT of amp*e^(i*phase), same grid/columns
      scan["point_spectra"]               -> ROI ratio DataFrame (rows = frequencies)
      scan["info"]                        -> parsed info.json
    """
    out = {"info": info, "combined": True, "frequencies": info.get("frequencies", [])}
    ds_meta = info.get("datasets", [])
    if ds_meta:
        out["pixelsize_nm"] = ds_meta[0].get("pixelsize_nm")

    def _table(key, cols):
        arr = z[key]
        if as_df and cols:
            try:
                import pandas as pd
                return pd.DataFrame(arr, columns=cols)
            except Exception:
                return arr
        return arr

    # per-frequency 2D maps from keys like "f860__O3A"
    maps = {}
    for key in z.files:
        m = re.fullmatch(r"f(?P<f>.+?)__(?P<ch>.+)", key)
        if not m:
            continue
        fl = m.group("f")
        try:
            fkey = float(fl)
        except ValueError:
            fkey = fl
        maps.setdefault(fkey, {})[m.group("ch")] = z[key]
    out["maps"] = maps or None

    def _channel_tables(profiles):
        out_d = {}
        for item in profiles:
            base = item.get("base") or item.get("key")
            keys = item.get("keys") or ({item["key"]: item["key"]} if item.get("key") else {})
            cols = item.get("columns")
            chans = {}
            for role, k in keys.items():
                if k in z.files:
                    chans[role] = _table(k, cols)
            if base and chans:
                out_d[base] = chans
        return out_d or None

    out["lineprofile"] = _channel_tables(info.get("lineprofile", {}).get("profiles", []))
    out["fft"] = _channel_tables(info.get("fft", {}).get("profiles", []))

    if "point_spectra" in z.files:
        out["point_spectra"] = _table("point_spectra", info.get("point_spectra", {}).get("columns"))
    else:
        out["point_spectra"] = None
    out["stats"] = None
    return Scan(out)


def loadnpz(path, *, include_channels_dict=False, lineprofile_as_dataframe=True):
    """
    Load a FwdBwd .npz summary file.

    Z and Z C arrays are stored in nm.  Line profiles are returned as
    scan["lineprofile"]["lp1"], scan["lineprofile"]["lp2"], ... .
    By default each lp# entry is a pandas DataFrame with columns from
    scan["info"]["lineprofile"]["profiles"], typically:
        distance_nm, O3A, O3P, Z_nm
    If pandas is unavailable or lineprofile_as_dataframe=False, entries
    remain NumPy arrays.

    Returns a Scan (a dict subclass) with built-in plotting that uses the pixel
    coordinates stored in the file — you never specify pixels:
        scan.show()              # image + crop + all ROI boxes + all lines
        scan.show_roi("P1")      # ROI box(es) drawn on the image
        scan.show_lp("L1")       # line-profile path(s) drawn on the image
        scan.plot_lp("L1")       # the profile curves vs distance (nm)

    Example
    -------
    scan = loadnpz("880cm-1_AVG.npz")
    O3A  = scan["O3A"]
    info = scan["info"]
    lp1  = scan["lineprofile"]["lp1"]   # DataFrame or ndarray
    px   = scan["pixelsize_nm"]
    scan.show()                          # quick overview of ROIs + line cuts
    """
    path = Path(path)
    z = np.load(path, allow_pickle=False)
    info0 = json.loads(_decode_npz_text(z["info"])) if "info" in z.files else {}
    if info0.get("combined"):
        return _load_combined(z, info0, lineprofile_as_dataframe)
    out = {}
    raw_lps = {}
    special = {"info", "stats", "point_spectra"}

    for key in z.files:
        if re.fullmatch(r"lp\d+", key):
            raw_lps[key] = z[key]
        elif key not in special:
            out[key] = z[key]

    out["info"] = info0
    info = out["info"]
    # 'pixelsize_nm' is the current key; fall back to the older 'pixel_size_nm'.
    px = info.get("pixelsize_nm", info.get("pixel_size_nm"))
    out["pixelsize_nm"] = px
    out["x_pixelsize_nm"] = info.get("x_pixelsize_nm", info.get("x_pixel_size_nm", px))
    out["y_pixelsize_nm"] = info.get("y_pixelsize_nm", info.get("y_pixel_size_nm", out["x_pixelsize_nm"]))

    if "stats" in z.files:
        text = _decode_npz_text(z["stats"]).strip()
        if text:
            import pandas as pd
            out["stats"] = pd.read_csv(io.StringIO(text))
        else:
            out["stats"] = None
    else:
        out["stats"] = None

    if "point_spectra" in z.files:
        arr = z["point_spectra"]
        cols = out["info"].get("point_spectra", {}).get("columns")
        if lineprofile_as_dataframe and cols:
            try:
                import pandas as pd
                out["point_spectra"] = pd.DataFrame(arr, columns=cols)
            except Exception:
                out["point_spectra"] = arr
        else:
            out["point_spectra"] = arr
    else:
        out["point_spectra"] = None

    lp_meta = {}
    for item in out["info"].get("lineprofile", {}).get("profiles", []):
        key = item.get("key")
        if key:
            lp_meta[key] = item

    lineprofiles = {}
    for key in sorted(raw_lps, key=_lp_sort_key):
        arr = raw_lps[key]
        cols = lp_meta.get(key, {}).get("columns")
        if lineprofile_as_dataframe and cols:
            try:
                import pandas as pd
                lineprofiles[key] = pd.DataFrame(arr, columns=cols)
            except Exception:
                lineprofiles[key] = arr
        else:
            lineprofiles[key] = arr

    out["lineprofile"] = lineprofiles or None

    if include_channels_dict:
        reserved = {"info", "stats", "lineprofile", "point_spectra",
                    "pixelsize_nm", "x_pixelsize_nm", "y_pixelsize_nm", "channels"}
        out["channels"] = {k: out[k] for k in out.keys() if k not in reserved}

    return Scan(out)


__all__ = ["loadnpz", "Scan"]
