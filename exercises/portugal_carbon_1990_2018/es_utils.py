"""Shared helpers for the Portugal carbon exercise.

Companion resource for Pedro Cabral, *Monitoring Dynamics in Ecosystem Service
Indicators with Remote Sensing* (Elsevier). Kept in one small module so the
notebooks can concentrate on the analytical decisions rather than plumbing.
"""

from __future__ import annotations

import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Patch

DATA_DIR = Path(__file__).parent / "data"
OUT_DIR = Path(__file__).parent / "outputs"

CLC_NODATA_CLASS = 48  # "NODATA" class inside the CLC raster
POOLS = ["c_above", "c_below", "c_soil", "c_dead"]

# ---------------------------------------------------------------------------
# CLC grid codes grouped into ten classes that matter for carbon.
# Grid codes (1-44) are what the CLC raster stores; they are not the
# three-digit CLC nomenclature codes (111, 112, ...).
# ---------------------------------------------------------------------------
GROUPS = {
    "Urban & artificial": range(1, 12),
    "Arable land": range(12, 15),
    "Permanent crops": range(15, 18),
    "Pastures & mixed agriculture": range(18, 22),
    "Agro-forestry (montado)": range(22, 23),
    "Forest": range(23, 26),
    "Shrubland & grassland": range(26, 30),
    "Open & burnt": range(30, 35),
    "Wetlands": range(35, 40),
    "Water": range(40, 45),
}
GROUP_ORDER = list(GROUPS)
GROUP_COLORS = ["#b2182b", "#f6e8a6", "#e0a33a", "#d9c96b", "#a6b35a",
                "#1b7837", "#7fbf7b", "#9e9e9e", "#8c6bb1", "#4393c3"]
_GROUP_OF_CODE = np.full(256, -1, dtype=np.int16)
for gi, codes in enumerate(GROUPS.values()):
    for c in codes:
        _GROUP_OF_CODE[c] = gi


def read_clc(path: Path) -> tuple[np.ndarray, dict]:
    """Read a clipped CLC raster; invalid pixels (nodata, class 48) become -1."""
    with rasterio.open(path) as src:
        arr = src.read(1).astype(np.int16)
        profile = src.profile.copy()
        if src.nodata is not None:
            arr[arr == int(src.nodata)] = -1
    arr[(arr == CLC_NODATA_CLASS) | (arr < 0) | (arr > 255)] = -1
    return arr, profile


def pixel_area_ha(profile: dict) -> float:
    t = profile["transform"]
    return abs(t.a * t.e) / 10_000


def to_groups(codes: np.ndarray) -> np.ndarray:
    """Map CLC grid codes to group index (0-9); invalid pixels stay -1."""
    out = np.full(codes.shape, -1, dtype=np.int16)
    valid = codes >= 0
    out[valid] = _GROUP_OF_CODE[codes[valid]]
    return out


def lookup(codes: np.ndarray, table: pd.DataFrame, column: str) -> np.ndarray:
    """Assign a table value to every pixel by class code (what InVEST does)."""
    lut = np.full(256, np.nan, dtype=np.float32)
    lut[table["lucode"].to_numpy(int)] = table[column].to_numpy(np.float32)
    out = np.full(codes.shape, np.nan, dtype=np.float32)
    valid = codes >= 0
    out[valid] = lut[codes[valid]]
    return out


def extent(profile: dict) -> list[float]:
    """Matplotlib extent [left, right, bottom, top] in km."""
    t, w, h = profile["transform"], profile["width"], profile["height"]
    return [t.c / 1e3, (t.c + w * t.a) / 1e3, (t.f + h * t.e) / 1e3, t.f / 1e3]


def plot_groups(ax, groups: np.ndarray, profile: dict, title: str = "", legend: bool = True):
    cmap = ListedColormap(GROUP_COLORS)
    cmap.set_bad("white")
    data = np.ma.masked_less(groups, 0)
    ax.imshow(data, cmap=cmap, norm=BoundaryNorm(np.arange(-0.5, 10.5), 10),
              extent=extent(profile), interpolation="nearest")
    ax.set_title(title)
    ax.set_xlabel("Easting (km, EPSG:3035)")
    ax.set_ylabel("Northing (km)")
    if legend:
        handles = [Patch(color=c, label=g) for g, c in zip(GROUP_ORDER, GROUP_COLORS)]
        ax.legend(handles=handles, loc="lower left", fontsize=7, frameon=True)


def plot_continuous(ax, data: np.ndarray, profile: dict, title: str, label: str,
                    cmap: str = "YlGn", vmin=None, vmax=None, diverging: bool = False):
    if diverging:
        lim = vmax if vmax is not None else np.nanpercentile(np.abs(data), 99)
        vmin, vmax = -lim, lim
    im = ax.imshow(data, cmap=cmap, vmin=vmin, vmax=vmax, extent=extent(profile),
                   interpolation="nearest")
    ax.set_title(title)
    ax.set_xticks([])
    ax.set_yticks([])
    plt.colorbar(im, ax=ax, shrink=0.7, label=label)
    return im


def area_by_group(groups: np.ndarray, area_ha: float) -> pd.Series:
    counts = np.bincount(groups[groups >= 0], minlength=len(GROUP_ORDER))
    return pd.Series(counts * area_ha, index=GROUP_ORDER, name="area_ha")


def transition_matrix(g_from: np.ndarray, g_to: np.ndarray, area_ha: float) -> pd.DataFrame:
    """Area (ha) moving from each group in year 1 (rows) to each group in year 2 (columns)."""
    n = len(GROUP_ORDER)
    valid = (g_from >= 0) & (g_to >= 0)
    idx = g_from[valid].astype(np.int64) * n + g_to[valid]
    m = np.bincount(idx, minlength=n * n).reshape(n, n) * area_ha
    return pd.DataFrame(m, index=GROUP_ORDER, columns=GROUP_ORDER)


# ---------------------------------------------------------------------------
# Trend statistics for Notebook 03 (Section 5.3.1)
# ---------------------------------------------------------------------------
def theil_sen_mann_kendall(cube: np.ndarray, years: np.ndarray, block_rows: int = 64):
    """Per-pixel Theil-Sen slope and Mann-Kendall Z for a (time, y, x) array.

    Pixels with fewer than 8 valid years return NaN. Processed in row blocks to
    keep memory low on a laptop.
    """
    t, h, w = cube.shape
    i, j = np.triu_indices(t, k=1)
    dt = (years[j] - years[i]).astype(np.float32)
    slope = np.full((h, w), np.nan, np.float32)
    z = np.full((h, w), np.nan, np.float32)
    for r0 in range(0, h, block_rows):
        blk = cube[:, r0:r0 + block_rows, :].astype(np.float32)  # (t, b, w)
        diff = blk[j] - blk[i]                                     # (pairs, b, w)
        slopes = diff / dt[:, None, None]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)  # all-NaN pixels (sea, no data)
            s_med = np.nanmedian(slopes, axis=0)
        sign_sum = np.nansum(np.sign(diff), axis=0)
        n = np.sum(~np.isnan(blk), axis=0).astype(np.float32)
        var = n * (n - 1) * (2 * n + 5) / 18.0
        with np.errstate(invalid="ignore", divide="ignore"):
            zz = np.where(sign_sum > 0, (sign_sum - 1) / np.sqrt(var),
                          np.where(sign_sum < 0, (sign_sum + 1) / np.sqrt(var), 0.0))
        enough = n >= 8
        slope[r0:r0 + block_rows] = np.where(enough, s_med, np.nan)
        z[r0:r0 + block_rows] = np.where(enough, zz, np.nan)
    return slope, z
