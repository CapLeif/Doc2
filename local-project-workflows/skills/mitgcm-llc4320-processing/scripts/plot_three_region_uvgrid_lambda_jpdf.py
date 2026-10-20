from __future__ import annotations

import csv
import os
from dataclasses import dataclass
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np


ROOT = Path(os.environ.get("MITGCM_PROJECT_ROOT", Path(__file__).resolve().parents[3]))
OUT = ROOT / "UVgrid_lambda_JPDF_masked"

X_LIMITS = (-10.0, 10.0)
Y_LIMITS = (-10.0, 10.0)
BIN_WIDTH = 0.05
LOG_MIN = np.log10(0.001)
OMEGA = 7.2921159e-5


@dataclass(frozen=True)
class RegionConfig:
    key: str
    title: str
    diagnostic_path: Path
    bbox: tuple[float, float, float, float] | None
    source_note: str
    local_f_source_path: Path | None = None


REGIONS = [
    RegionConfig(
        key="Kuroshio_Extension",
        title="Kuroshio Extension",
        diagnostic_path=ROOT
        / "KuroshioExt"
        / "lambda_theta_surface_20120101T050000"
        / "Kuroshio_surface_diagnostics.npz",
        bbox=(155.0, 165.0, 36.0, 40.0),
        source_note="KuroshioExt UV-grid diagnostics; mask 36-40N, 155-165E",
    ),
    RegionConfig(
        key="South_China_Sea",
        title="South China Sea",
        diagnostic_path=ROOT
        / "SouthChinaSea3"
        / "lambda_theta_surface_0001492560_omega_positive_gray_v2"
        / "South_China_Sea_surface_diagnostics.npz",
        bbox=(115.0, 125.0, 17.0, 23.0),
        source_note="SouthChinaSea3 UV-grid diagnostics; mask 17-23N, 115-125E",
    ),
    RegionConfig(
        key="West_Atlantic",
        title="West Atlantic",
        diagnostic_path=ROOT / "Westatlantic" / "llc4320_flow_diagnostics_20120913_00UTC.npz",
        bbox=None,
        source_note="Westatlantic UV-grid diagnostics; full-scale mask; f recomputed from local latitude",
        local_f_source_path=Path("Westatlantic") / "LLC4320_pre-SWOT_WestAtlantic_20120913.nc",
    ),
]


def as_float_array(value) -> np.ndarray:
    if np.ma.isMaskedArray(value):
        return value.filled(np.nan).astype(float)
    return np.asarray(value, dtype=float)


def geographic_gradient(q: np.ndarray, lon: np.ndarray, lat: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    radius = 6371000.0
    lon_r = np.deg2rad(lon)
    lat_r = np.deg2rad(lat)
    dlon_i = np.gradient(lon_r, axis=0)
    dlon_j = np.gradient(lon_r, axis=1)
    dlat_i = np.gradient(lat_r, axis=0)
    dlat_j = np.gradient(lat_r, axis=1)
    coslat = np.cos(lat_r)
    dx_i = radius * coslat * dlon_i
    dx_j = radius * coslat * dlon_j
    dy_i = radius * dlat_i
    dy_j = radius * dlat_j
    dq_i = np.gradient(q, axis=0)
    dq_j = np.gradient(q, axis=1)
    det = dx_i * dy_j - dx_j * dy_i
    good = np.isfinite(det) & (np.abs(det) > 1e-3)
    qx = np.divide(dq_i * dy_j - dq_j * dy_i, det, out=np.full_like(q, np.nan), where=good)
    qy = np.divide(dx_i * dq_j - dx_j * dq_i, det, out=np.full_like(q, np.nan), where=good)
    return qx, qy


def load_west_atlantic_local_f(config: RegionConfig) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    if config.local_f_source_path is None:
        raise ValueError("West Atlantic local-f source path is required.")
    source_path = config.local_f_source_path
    if not source_path.is_absolute():
        source_path = ROOT / source_path
    with nc.Dataset(source_path) as ds:
        lon = as_float_array(ds.variables["XC"][:-1, :-1])
        lat = as_float_array(ds.variables["YC"][:-1, :-1])
        angle_cs = as_float_array(ds.variables["AngleCS"][:-1, :-1])
        angle_sn = as_float_array(ds.variables["AngleSN"][:-1, :-1])
        depth = as_float_array(ds.variables["Depth"][:-1, :-1])
        u_face = as_float_array(ds.variables["U"][0, 0, :-1, :])
        v_face = as_float_array(ds.variables["V"][0, 0, :, :-1])

    u_grid = 0.5 * (u_face[:, :-1] + u_face[:, 1:])
    v_grid = 0.5 * (v_face[:-1, :] + v_face[1:, :])
    angle_ok = np.isfinite(angle_cs) & np.isfinite(angle_sn) & (np.abs(np.hypot(angle_cs, angle_sn) - 1.0) < 1e-3)
    u = np.where(angle_ok, angle_cs * u_grid - angle_sn * v_grid, np.nan)
    v = np.where(angle_ok, angle_sn * u_grid + angle_cs * v_grid, np.nan)
    wet = (
        (depth > 0)
        & np.isfinite(lon)
        & np.isfinite(lat)
        & np.isfinite(u)
        & np.isfinite(v)
        & (np.abs(u) < 20)
        & (np.abs(v) < 20)
    )
    u = np.where(wet, u, np.nan)
    v = np.where(wet, v, np.nan)
    ux, uy = geographic_gradient(u, lon, lat)
    vx, vy = geographic_gradient(v, lon, lat)
    zeta = vx - uy
    sigma = np.hypot(ux - vy, vx + uy)
    delta = ux + vy
    f = 2.0 * OMEGA * np.sin(np.deg2rad(lat))
    good_f = np.isfinite(f) & (np.abs(f) > 1e-8)
    div_f = np.divide(delta, f, out=np.full_like(delta, np.nan), where=good_f)
    sigma_f = np.divide(sigma, np.abs(f), out=np.full_like(sigma, np.nan), where=good_f)
    ro = np.divide(zeta, f, out=np.full_like(zeta, np.nan), where=good_f)
    land = ~wet
    return lon, lat, div_f, sigma_f, ro, land


def load_jpdf_arrays(config: RegionConfig) -> dict[str, np.ndarray]:
    if config.local_f_source_path is not None and config.local_f_source_path.suffix.lower() == ".nc":
        lon, lat, div_f, sigma_f, ro, land = load_west_atlantic_local_f(config)
    else:
        with np.load(config.diagnostic_path) as z:
            land = z["land"].astype(bool) if "land" in z.files else None
            lon = z["lon"].astype(float)
            lat = z["lat"].astype(float)
            div_f = z["div_f"].astype(float)
            sigma_f = z["sigma_f"].astype(float)
            ro = z["Ro"].astype(float)
        if land is None:
            land = np.zeros(lon.shape, dtype=bool)

    coord = np.isfinite(lon) & np.isfinite(lat) & (np.abs(lon) <= 180) & (np.abs(lat) <= 90)
    if config.bbox is None:
        roi = coord
    else:
        lon_min, lon_max, lat_min, lat_max = config.bbox
        roi = coord & (lon >= lon_min) & (lon <= lon_max) & (lat >= lat_min) & (lat <= lat_max)

    wet = roi & ~land & np.isfinite(div_f) & np.isfinite(sigma_f) & np.isfinite(ro)
    omega_f2 = sigma_f * sigma_f - ro * ro
    y = np.full(div_f.shape, np.nan)
    real = wet & (omega_f2 >= 0)
    complex_pair = wet & (omega_f2 < 0)
    y[real] = -np.sqrt(omega_f2[real])
    y[complex_pair] = np.sqrt(-omega_f2[complex_pair])
    valid = wet & np.isfinite(div_f) & np.isfinite(y)
    return {
        "x": div_f[valid],
        "y": y[valid],
        "real": real,
        "complex": complex_pair,
        "roi": roi,
        "wet": wet,
    }


def make_lambda_jpdf(x: np.ndarray, y: np.ndarray) -> dict[str, np.ndarray | float | int]:
    x_edges = np.arange(X_LIMITS[0], X_LIMITS[1] + BIN_WIDTH * 0.5, BIN_WIDTH)
    y_edges = np.arange(Y_LIMITS[0], Y_LIMITS[1] + BIN_WIDTH * 0.5, BIN_WIDTH)
    counts, _, _ = np.histogram2d(x, y, bins=(x_edges, y_edges))
    in_range_count = int(counts.sum())
    pdf = counts / (in_range_count * BIN_WIDTH * BIN_WIDTH) if in_range_count else counts

    positive = pdf > 0
    log_pdf = np.full_like(pdf, np.nan, dtype=float)
    log_pdf[positive] = np.log10(pdf[positive])
    if np.any(positive):
        pdf_max = float(pdf[positive].max())
        log_max = float(np.log10(pdf_max))
    else:
        pdf_max = np.nan
        log_max = 0.0
    zero_level = LOG_MIN - max(0.01, (log_max - LOG_MIN) / 255.0)
    plot_values = log_pdf.copy()
    plot_values[positive] = np.maximum(plot_values[positive], LOG_MIN)
    plot_values[~positive] = zero_level
    return {
        "x_edges": x_edges,
        "y_edges": y_edges,
        "x_centers": 0.5 * (x_edges[:-1] + x_edges[1:]),
        "y_centers": 0.5 * (y_edges[:-1] + y_edges[1:]),
        "counts": counts,
        "pdf": pdf,
        "plot_values": plot_values,
        "in_range_count": in_range_count,
        "pdf_max": pdf_max,
        "log_max": log_max,
        "zero_level": zero_level,
    }


def draw_jpdf(config: RegionConfig, arrays: dict[str, np.ndarray], out_png: Path) -> dict[str, float | int | str]:
    x = arrays["x"]
    y = arrays["y"]
    finite = np.isfinite(x) & np.isfinite(y)
    x = x[finite]
    y = y[finite]
    jpdf = make_lambda_jpdf(x, y)

    fig, ax = plt.subplots(figsize=(6.0, 4.5), dpi=200)
    hm = ax.imshow(
        np.asarray(jpdf["plot_values"]).T,
        origin="lower",
        extent=(X_LIMITS[0], X_LIMITS[1], Y_LIMITS[0], Y_LIMITS[1]),
        interpolation="nearest",
        aspect="equal",
    )

    raw_cmap = mpl.colormaps["turbo"](np.linspace(0, 1, 255))
    cmap = mpl.colors.ListedColormap(np.vstack(([1, 1, 1, 1], raw_cmap)))
    hm.set_cmap(cmap)
    hm.set_clim(float(jpdf["zero_level"]), float(jpdf["log_max"]))
    cb = fig.colorbar(hm, ax=ax, pad=0.018)
    cb.set_ticks(np.log10([0.001, 0.01, 0.1, 1.0]))
    cb.set_ticklabels(["0.001", "0.01", "0.1", "1"])

    grid = np.linspace(X_LIMITS[0], X_LIMITS[1], 401)
    xx, yy = np.meshgrid(grid, grid)
    lambda_minus = 0.5 * (xx + yy)
    lambda_minus = np.where(yy < 0, lambda_minus, np.nan)
    lambda_minus_levels = np.arange(-10.0, 1.0, 1.0)
    ax.contour(
        xx,
        yy,
        lambda_minus,
        lambda_minus_levels,
        colors=[(0.5, 0.5, 0.5)],
        linewidths=0.65,
        linestyles="--",
    )
    for level in lambda_minus_levels:
        label_x = 2.0 * level + 9.0
        if X_LIMITS[0] <= label_x <= X_LIMITS[1]:
            ax.text(
                label_x,
                -9.0,
                f"{level:g}",
                color=(0.35, 0.35, 0.35),
                fontsize=6.5,
                ha="center",
                va="center",
                bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.5},
            )

    ax.plot([X_LIMITS[0], 0], [Y_LIMITS[0], 0], "k--", linewidth=0.8)
    ax.plot([0, X_LIMITS[1]], [0, Y_LIMITS[0]], "k--", linewidth=0.8)

    lambda_plus = 0.5 * (xx - yy)
    lambda_plus = np.where((xx < 0) & (yy < 0) & (yy > xx), lambda_plus, np.nan)
    lambda_plus_levels = np.arange(-10.0, 1.0, 1.0)
    ax.contour(
        xx,
        yy,
        lambda_plus,
        lambda_plus_levels,
        colors=[(0.85, 0.55, 0.55)],
        linewidths=0.65,
        linestyles="--",
    )
    for level in lambda_plus_levels:
        label_y = -9.0 - 2.0 * level
        if Y_LIMITS[0] < label_y < 0 and label_y > -9.0:
            ax.text(
                -9.0,
                label_y,
                f"{level:g}",
                color=(0.75, 0.35, 0.35),
                fontsize=6.5,
                ha="center",
                va="center",
                bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.5},
            )

    ax.axvline(0, color="k", linestyle="--", linewidth=0.8)
    ax.axhline(0, color="k", linestyle="--", linewidth=0.8)
    ax.set_xlim(*X_LIMITS)
    ax.set_ylim(*Y_LIMITS)
    ax.set_xlabel(r"$\Delta/f$", fontsize=9)
    ax.set_ylabel(r"$-\sqrt{\Omega}/|f|$", fontsize=9)
    ax.set_title(r"$\lambda$-JPDF | " + config.title, fontsize=9.5, fontweight="normal")
    ax.tick_params(direction="out", labelsize=8)
    ax.set_facecolor("white")
    for spine in ax.spines.values():
        spine.set_linewidth(0.8)

    fig.savefig(out_png, dpi=600, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    return {
        "region": config.title,
        "source": str(config.diagnostic_path),
        "mask": config.source_note,
        "roi_cells": int(np.count_nonzero(arrays["roi"])),
        "wet_cells": int(np.count_nonzero(arrays["wet"])),
        "valid_jpdf_points": int(x.size),
        "in_range_points": int(jpdf["in_range_count"]),
        "real_pair_points": int(np.count_nonzero(arrays["real"])),
        "complex_pair_points": int(np.count_nonzero(arrays["complex"])),
        "pdf_max": float(jpdf["pdf_max"]),
        "output_png": str(out_png),
    }


def main() -> None:
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "axes.linewidth": 0.8,
            "mathtext.fontset": "dejavusans",
        }
    )
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for config in REGIONS:
        arrays = load_jpdf_arrays(config)
        out_png = OUT / f"{config.key}_UVgrid_masked_lambda_JPDF.png"
        rows.append(draw_jpdf(config, arrays, out_png))
        print(f"{config.title}: {out_png}")
        print(
            f"  valid={rows[-1]['valid_jpdf_points']:,}, "
            f"in_range={rows[-1]['in_range_points']:,}, "
            f"real={rows[-1]['real_pair_points']:,}, "
            f"complex={rows[-1]['complex_pair_points']:,}"
        )

    csv_path = OUT / "UVgrid_masked_lambda_JPDF_summary.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Summary: {csv_path}")


if __name__ == "__main__":
    main()

