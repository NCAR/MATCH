#!/usr/bin/env python3
"""
Monthly-mean AOD maps from MATCH output (CORe and/or CDAS).

Plots total and species-resolved aerosol optical depth on Robinson projection
using NCAR brand styling (adapted from CDAS-CORe / DAVINCI-MONET).

Usage:
    python plot_aod_core.py                                    # CORe only
    python plot_aod_core.py --ref-dir ~/Data/MATCH/.../cdas     # + CDAS & diffs
    python plot_aod_core.py -o ~/plots --label "Feb 2026"
"""

import argparse
import logging
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.util import add_cyclic_point
import xarray as xr

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

# ── NCAR style (adapted from DAVINCI-MONET / CDAS-CORe) ────────────────

NCAR_COLORS = {
    "space": "#011837",
    "dark_blue": "#00357A",
    "ncar_blue": "#0A5DDA",
    "aqua": "#00A2B4",
    "light_blue": "#CEDFF8",
    "light_gray": "#F1F0EE",
    "orange": "#FF8C00",
    "yellow": "#FFDD31",
    "gray": "#58595B",
    "red": "#D62839",
    "green": "#2E8B57",
    "purple": "#7B68EE",
}

NCAR_PALETTE = [
    NCAR_COLORS["ncar_blue"], NCAR_COLORS["aqua"], NCAR_COLORS["orange"],
    NCAR_COLORS["purple"], NCAR_COLORS["green"], NCAR_COLORS["red"],
    NCAR_COLORS["yellow"], NCAR_COLORS["dark_blue"],
]


def apply_ncar_style():
    """Apply NCAR brand styling globally."""
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Poppins", "Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"],
        "mathtext.fontset": "dejavusans",
        "axes.labelsize": 12,
        "axes.titlesize": 14,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.fontsize": 10,
        "figure.titlesize": 16,
        "axes.prop_cycle": plt.cycler(color=NCAR_PALETTE),
        "lines.linewidth": 1.5,
        "axes.grid": True,
        "grid.alpha": 0.3,
        "grid.linestyle": "-",
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
    })


def save_figure(fig, name, output_dir):
    """Save figure as PNG (300 DPI) and PDF."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        path = output_dir / f"{name}.{ext}"
        fig.savefig(path, dpi=300, bbox_inches="tight")
    log.info(f"Saved: {output_dir / name}.png/pdf")


# ── Plot functions ──────────────────────────────────────────────────────

def _aod_levels(vmax):
    """Non-uniform AOD contour levels (CERES-SARB / DAVINCI convention).

    Fine resolution at low values, then coarse steps above 0.04, capped at *vmax*.
    For small vmax (< 0.25) use 0.01 coarse steps to keep enough contour detail.
    """
    fine = np.array([0.00, 0.005, 0.01, 0.02, 0.03, 0.04])
    step = 0.01 if vmax < 0.25 else 0.05
    coarse = np.arange(0.05, vmax + step / 2, step)
    levels = np.concatenate([fine, coarse])
    levels = levels[levels <= vmax + 1e-9]
    return levels


def plot_aod_map(data, lats, lons, title, vmin, vmax, output_dir, filename,
                 cmap="turbo", nlevels=20):
    """Global map of a 2-D AOD field on Robinson projection."""
    data_cyc, lons_cyc = add_cyclic_point(data, coord=lons)

    fig = plt.figure(figsize=(12, 5), constrained_layout=True)
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson())

    levels = _aod_levels(vmax)
    cf = ax.contourf(lons_cyc, lats, data_cyc, levels=levels, cmap=cmap,
                     extend="max", transform=ccrs.PlateCarree())
    cf.set_rasterized(True)

    ax.coastlines(linewidth=0.5, color=NCAR_COLORS["gray"])
    ax.add_feature(cfeature.BORDERS, linewidth=0.3, edgecolor=NCAR_COLORS["gray"])
    ax.set_global()
    ax.set_title(title)

    cbar = fig.colorbar(cf, ax=ax, orientation="horizontal", pad=0.05,
                        shrink=0.7, aspect=35)
    cbar.set_label("Optical Depth")

    save_figure(fig, filename, output_dir)
    plt.close(fig)


def plot_diff_map(data, lats, lons, title, vmax_diff, output_dir, filename,
                  nlevels=20):
    """Global map of a difference field (diverging colormap)."""
    data_cyc, lons_cyc = add_cyclic_point(data, coord=lons)

    fig = plt.figure(figsize=(12, 5), constrained_layout=True)
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson())

    levels = np.arange(-vmax_diff, vmax_diff + 0.025, 0.05)
    cf = ax.contourf(lons_cyc, lats, data_cyc, levels=levels, cmap="RdBu_r",
                     extend="both", transform=ccrs.PlateCarree())
    cf.set_rasterized(True)

    ax.coastlines(linewidth=0.5, color=NCAR_COLORS["gray"])
    ax.add_feature(cfeature.BORDERS, linewidth=0.3, edgecolor=NCAR_COLORS["gray"])
    ax.set_global()
    ax.set_title(title)

    cbar = fig.colorbar(cf, ax=ax, orientation="horizontal", pad=0.05,
                        shrink=0.7, aspect=35)
    cbar.set_label("Δ Optical Depth")

    save_figure(fig, filename, output_dir)
    plt.close(fig)


# ── AOD field definitions ───────────────────────────────────────────────

AOD_FIELDS = [
    {
        "name": "Total AOD",
        "varname": "AEROD",
        "filename": "aod_total",
        "vmin": 0.0,
        "vmax": 0.7,
        "vmax_diff": 0.3,
    },
    {
        "name": "Sulfate AOD",
        "varname": "SO4OD",
        "filename": "aod_sulfate",
        "vmin": 0.0,
        "vmax": 0.3,
        "vmax_diff": 0.1,
    },
    {
        "name": "Dust AOD",
        "varname": "DSTODXC",
        "filename": "aod_dust",
        "vmin": 0.0,
        "vmax": 0.8,
        "vmax_diff": 0.3,
    },
    {
        "name": "Sea Salt AOD",
        "varname": "SSLTOD",
        "filename": "aod_seasalt",
        "vmin": 0.0,
        "vmax": 0.15,
        "vmax_diff": 0.05,
    },
    {
        "name": "Carbonaceous AOD",
        "varname": None,  # computed as sum of OC + BC
        "components": ["OCPHOOD", "OCPHIOD", "BCPHOOD", "BCPHIOD"],
        "filename": "aod_carbon",
        "vmin": 0.0,
        "vmax": 0.2,
        "vmax_diff": 0.1,
    },
]


# ── Data loading ────────────────────────────────────────────────────────

def load_aod_means(input_dir):
    """Accumulate running means from h*.nc files. Returns (means_dict, lats, lons)."""
    input_dir = Path(input_dir)
    nc_files = sorted(input_dir.glob("h[0-9][0-9][0-9][0-9].nc"))
    if not nc_files:
        log.error(f"No netCDF files found in {input_dir}")
        return None, None, None

    aod_varnames = set()
    for f in AOD_FIELDS:
        if f["varname"] is not None:
            aod_varnames.add(f["varname"])
        else:
            aod_varnames.update(f["components"])

    log.info(f"Accumulating means from {len(nc_files)} netCDF files in {input_dir}")

    sums = {v: None for v in aod_varnames}
    count = 0

    for ncf in nc_files:
        with xr.open_dataset(ncf) as ds1:
            if count == 0:
                lats = ds1.lat.values
                lons = ds1.lon.values
            for v in aod_varnames:
                vals = ds1[v].values.squeeze()
                sums[v] = vals if sums[v] is None else sums[v] + vals
            count += 1
        if count % 100 == 0:
            log.info(f"  ... processed {count}/{len(nc_files)} files")

    means = {v: sums[v] / count for v in aod_varnames}
    log.info(f"Computed monthly means from {count} timesteps")
    return means, lats, lons


def get_field_data(means, field):
    """Extract or compute the AOD data for a field definition."""
    if field["varname"] is not None:
        return means[field["varname"]]
    return sum(means[v] for v in field["components"])


# ── Main ────────────────────────────────────────────────────────────────

def generate_all_plots(input_dir, ref_dir, output_dir, label,
                       run_name="CORe", ref_name="CDAS", ref_label=None):
    """Generate AOD maps, and optionally reference + difference plots."""
    if ref_label is None:
        ref_label = label

    means_test, lats, lons = load_aod_means(input_dir)
    if means_test is None:
        return

    for field in AOD_FIELDS:
        log.info(f"Plotting {run_name} {field['name']}...")
        data = get_field_data(means_test, field)
        plot_aod_map(data, lats, lons,
                     title=f"{field['name']} — MATCH {run_name} {label}",
                     vmin=field["vmin"], vmax=field["vmax"],
                     output_dir=output_dir, filename=field["filename"])

    # Load reference data if provided
    if ref_dir is None:
        return

    means_ref, lats_ref, lons_ref = load_aod_means(ref_dir)
    if means_ref is None:
        return

    for field in AOD_FIELDS:
        name = field["name"]
        data_ref = get_field_data(means_ref, field)

        # Reference absolute plot
        log.info(f"Plotting {ref_name} {name}...")
        plot_aod_map(data_ref, lats_ref, lons_ref,
                     title=f"{name} — MATCH {ref_name} {ref_label}",
                     vmin=field["vmin"], vmax=field["vmax"],
                     output_dir=output_dir,
                     filename=f"{field['filename']}_{ref_name.lower()}")

        # Difference plot
        log.info(f"Plotting Δ{name} ({run_name} − {ref_name})...")
        data_test = get_field_data(means_test, field)
        diff = data_test - data_ref
        plot_diff_map(diff, lats, lons,
                      title=f"Δ{name} ({run_name} − {ref_name}) — {label}",
                      vmax_diff=field["vmax_diff"],
                      output_dir=output_dir,
                      filename=f"{field['filename']}_diff")

    log.info(f"All plots saved to {output_dir}")


def main():
    parser = argparse.ArgumentParser(description="MATCH AOD maps (CORe vs CDAS)")
    parser.add_argument("--input-dir",
                        default=str(Path.home() / "Data/MATCH/2026/sarb4_hourly_202602_core"),
                        help="CORe run directory containing h*.nc files")
    parser.add_argument("--ref-dir",
                        default=None,
                        help="CDAS reference run directory containing h*.nc files")
    parser.add_argument("-o", "--output",
                        default=None,
                        help="Output directory (default: <input-dir>/plots)")
    parser.add_argument("--label", default="Feb 2026",
                        help="Label for plot titles")
    parser.add_argument("--run-name", default="CORe",
                        help="Run name for plot titles (e.g. CORe, CDAS)")
    parser.add_argument("--ref-name", default="CDAS",
                        help="Reference run name for plot titles")
    parser.add_argument("--ref-label", default=None,
                        help="Reference label (default: same as --label)")
    args = parser.parse_args()

    output_dir = args.output or str(Path(args.input_dir) / "plots")
    ref_label = args.ref_label or args.label

    apply_ncar_style()
    generate_all_plots(args.input_dir, args.ref_dir, output_dir, args.label,
                       run_name=args.run_name, ref_name=args.ref_name,
                       ref_label=ref_label)


if __name__ == "__main__":
    main()
