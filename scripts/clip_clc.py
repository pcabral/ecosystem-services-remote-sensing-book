"""Clip CORINE Land Cover rasters to mainland Portugal.

Companion resource for Pedro Cabral, *Monitoring Dynamics in Ecosystem Service
Indicators with Remote Sensing* (Elsevier).

The full-Europe CLC rasters (100 m, EPSG:3035) are too large to keep in a
teaching repository. This script cuts them to mainland Portugal and writes small,
compressed GeoTIFFs that the exercise notebooks use.

Book link (Section 3.3): CLC is a categorical product, so the script never
reprojects or resamples the raster. Instead, the boundary *vector* is
reprojected to the raster's CRS and the raster is cut on its original grid.
The class codes, pixel size and pixel alignment therefore stay exactly as
published, and the 1990 and 2018 outputs share one grid.

Usage (from the repository root)::

    python scripts/clip_clc.py \
        --clc1990 path/to/U2000_CLC1990_V2020_20u1.tif \
        --clc2018 path/to/U2018_CLC2018_V2020_20u1.tif \
        --boundary path/to/Portugal.shp \
        --caop path/to/Continente_CAOP2025.gpkg

``--caop`` is optional. When given, the municipality layer is extracted,
simplified and saved so that Notebook 02 can aggregate results by municipality.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from rasterio.mask import mask

DEFAULT_OUT = Path("exercises/portugal_carbon_1990_2018/data")


def sha256(path: Path, chunk: int = 2**20) -> str:
    """Checksum of an input file, recorded for reproducibility (Section 3.6)."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def fallback_nodata(dtype: str):
    """Value used for pixels outside the boundary if the source declares none."""
    return {"int8": -128, "uint8": 128, "int16": -32768, "uint16": 65535}.get(dtype, 0)


def clip_raster(src_path: Path, boundary: gpd.GeoDataFrame, out_path: Path) -> dict:
    """Cut one categorical raster to the boundary, keeping its native grid."""
    with rasterio.open(src_path) as src:
        geoms = boundary.to_crs(src.crs).geometry  # reproject the vector, not the raster
        nodata = src.nodata if src.nodata is not None else fallback_nodata(src.dtypes[0])
        data, transform = mask(src, geoms, crop=True, nodata=nodata, filled=True,
                               all_touched=False)
        profile = {
            "driver": "GTiff",
            "dtype": src.dtypes[0],
            "count": 1,
            "crs": src.crs,
            "transform": transform,
            "width": data.shape[2],
            "height": data.shape[1],
            "nodata": nodata,
            "compress": "deflate",
            "tiled": True,
            "blockxsize": 256,
            "blockysize": 256,
        }
        try:
            colormap = src.colormap(1)  # keep the official CLC colour table if present
        except ValueError:
            colormap = None

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with rasterio.open(out_path, "w", **profile) as dst:
        dst.write(data[0], 1)
        if colormap:
            dst.write_colormap(1, colormap)
        dst.update_tags(source=Path(src_path).name,
                        note="Clipped to mainland Portugal on the native CLC grid; no resampling.")

    values, counts = np.unique(data[0][data[0] != nodata], return_counts=True)
    return {
        "source": Path(src_path).name,
        "source_sha256": sha256(Path(src_path)),
        "output": out_path.name,
        "crs": str(profile["crs"]),
        "pixel_size_m": [transform.a, -transform.e],
        "shape": [profile["height"], profile["width"]],
        "nodata": nodata,
        "class_pixel_counts": {int(v): int(c) for v, c in zip(values, counts)},
    }


def check_codes(reports: list[dict], carbon_table: Path) -> dict:
    """Compare raster class codes with the carbon table's lucode column (Section 4.6.3)."""
    table_codes = set(pd.read_csv(carbon_table)["lucode"].astype(int))
    result = {}
    for r in reports:
        raster_codes = set(r["class_pixel_counts"])
        result[r["output"]] = {
            "codes_missing_from_table": sorted(raster_codes - table_codes),
            "table_codes_absent_from_raster": sorted(table_codes - raster_codes),
        }
    return result


def extract_municipalities(caop_path: Path, out_path: Path, layer: str | None,
                           tolerance_m: float = 50.0) -> dict:
    """Extract and simplify the CAOP municipality layer for Notebook 02 (Section 3.5)."""
    import pyogrio

    layers = [name for name, _ in pyogrio.list_layers(caop_path)]
    if layer is None:
        candidates = [n for n in layers if "munic" in n.lower()]
        if not candidates:
            raise SystemExit(f"No municipality layer found in {caop_path.name}. "
                             f"Layers available: {layers}. Pass one with --caop-layer.")
        layer = candidates[0]

    gdf = gpd.read_file(caop_path, layer=layer)
    name_cols = [c for c in gdf.columns if c.lower() in ("municipio", "concelho", "nome", "name")]
    code_cols = [c for c in gdf.columns if c.lower() in ("dtmn", "dico", "codigo", "code")]
    keep = (code_cols[:1] + name_cols[:1]) or [c for c in gdf.columns if c != "geometry"][:2]
    gdf = gdf[keep + ["geometry"]].rename(columns=dict(zip(keep, ["code", "name"][: len(keep)])))
    if "name" not in gdf.columns:
        gdf["name"] = gdf.iloc[:, 0].astype(str)
    gdf = gdf.dissolve(by="name", as_index=False)  # one feature per municipality
    gdf["geometry"] = gdf.geometry.simplify(tolerance_m, preserve_topology=True)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    gdf.to_file(out_path, driver="GPKG", layer="municipalities")
    return {"source": caop_path.name, "layer": layer, "n_municipalities": len(gdf),
            "simplify_tolerance_m": tolerance_m, "output": out_path.name}


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--clc1990", type=Path, required=True)
    p.add_argument("--clc2018", type=Path, required=True)
    p.add_argument("--boundary", type=Path, required=True,
                   help="Mainland Portugal outline (any vector format, any CRS).")
    p.add_argument("--caop", type=Path, help="Optional CAOP GeoPackage for municipalities.")
    p.add_argument("--caop-layer", help="Municipality layer name if auto-detection fails.")
    p.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = p.parse_args()

    boundary = gpd.read_file(args.boundary)
    boundary = boundary.dissolve()  # a single mainland outline
    args.out.mkdir(parents=True, exist_ok=True)
    boundary.to_file(args.out / "portugal_mainland.gpkg", driver="GPKG", layer="mainland")

    reports = [
        clip_raster(args.clc1990, boundary, args.out / "clc1990_portugal.tif"),
        clip_raster(args.clc2018, boundary, args.out / "clc2018_portugal.tif"),
    ]
    a, b = reports
    if a["shape"] != b["shape"] or a["crs"] != b["crs"]:
        raise SystemExit("The two clipped rasters do not share a grid; check the inputs.")

    report = {"rasters": reports}
    table = args.out / "carbon_pools.csv"
    if table.exists():
        report["code_check"] = check_codes(reports, table)
    if args.caop:
        report["municipalities"] = extract_municipalities(
            args.caop, args.out / "municipalities.gpkg", args.caop_layer)

    (args.out / "clip_report.json").write_text(json.dumps(report, indent=2))
    for r in reports:
        size_mb = (args.out / r["output"]).stat().st_size / 1e6
        print(f"{r['output']}: {r['shape'][0]} x {r['shape'][1]} pixels, {size_mb:.1f} MB")
    if "code_check" in report:
        for name, chk in report["code_check"].items():
            print(f"{name}: codes missing from carbon table -> {chk['codes_missing_from_table']}")
    print(f"Report written to {args.out / 'clip_report.json'}")


if __name__ == "__main__":
    main()
