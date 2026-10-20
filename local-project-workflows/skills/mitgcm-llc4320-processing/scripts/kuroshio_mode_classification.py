from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from kuroshio_fl_theta_gamma_epsilon import (
    Front,
    extract_fronts,
    load_surface_fields,
)
from regional_fl_jpdf import load_binary_region, load_westatlantic_nc


MODE_NAMES = ("Mode1", "Mode2", "Mode3", "Mode4", "Unclassified")
MODE_COLORS = {
    "Mode1": "#3f65b7",
    "Mode2": "#c95047",
    "Mode3": "#4b9460",
    "Mode4": "#c1842c",
    "Unclassified": "#8f8f8f",
}


def eigenvector_for_lambda(a: np.ndarray, b: np.ndarray, c: np.ndarray, d: np.ndarray, lam: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    vx = b.copy()
    vy = lam - a
    n = np.hypot(vx, vy)
    alt = n <= 1e-14
    vx[alt] = lam[alt] - d[alt]
    vy[alt] = c[alt]
    n = np.hypot(vx, vy)
    good = n > 0
    out_x = np.full_like(vx, np.nan, dtype=float)
    out_y = np.full_like(vy, np.nan, dtype=float)
    out_x[good] = vx[good] / n[good]
    out_y[good] = vy[good] / n[good]
    return out_x, out_y


def acute_axis_angle_deg(ax: np.ndarray, ay: np.ndarray, bx: np.ndarray, by: np.ndarray) -> np.ndarray:
    dot = np.abs(ax * bx + ay * by)
    dot = np.clip(dot, 0.0, 1.0)
    return np.rad2deg(np.arccos(dot))


def signed_axis_angle_deg(
    ref_x: np.ndarray,
    ref_y: np.ndarray,
    axis_x: np.ndarray,
    axis_y: np.ndarray,
    zeta: np.ndarray,
) -> np.ndarray:
    cross = ref_x * axis_y - ref_y * axis_x
    dot = ref_x * axis_x + ref_y * axis_y
    ccw = np.rad2deg(np.arctan2(cross, dot))
    ccw_acute = (ccw + 90.0) % 180.0 - 90.0
    out = -np.sign(zeta) * ccw_acute
    out[np.sign(zeta) == 0] = np.nan
    return out


def compute_full_field_angles(data: dict[str, np.ndarray], bbox: tuple[float, float, float, float] | None) -> dict[str, np.ndarray]:
    lon = data["lon"]
    lat = data["lat"]
    if bbox is None:
        region = np.isfinite(lon) & np.isfinite(lat)
    else:
        lon_min, lon_max, lat_min, lat_max = bbox
        region = (lon >= lon_min) & (lon <= lon_max) & (lat >= lat_min) & (lat <= lat_max)
    base = data["valid"] & region

    ux = data["ux"]
    uy = data["uy"]
    vx = data["vx"]
    vy = data["vy"]
    bx = data["bx"]
    by = data["by"]
    f = data["f"]

    shape = ux.shape
    out = {
        "theta_deg": np.full(shape, np.nan, dtype=np.float32),
        "theta_acute_deg": np.full(shape, np.nan, dtype=np.float32),
        "alpha_deg": np.full(shape, np.nan, dtype=np.float32),
        "gamma_signed_deg": np.full(shape, np.nan, dtype=np.float32),
        "epsilon_signed_deg": np.full(shape, np.nan, dtype=np.float32),
        "lambda_at_minus_over_absf": np.full(shape, np.nan, dtype=np.float32),
        "lambda_at_plus_over_absf": np.full(shape, np.nan, dtype=np.float32),
        "zeta": np.full(shape, np.nan, dtype=np.float32),
        "valid_angle": np.zeros(shape, dtype=bool),
        "region": region,
    }

    sigma = np.hypot(ux - vy, uy + vx)
    zeta = vx - uy
    delta = ux + vy
    det = ux * vy - uy * vx
    disc = delta * delta - 4.0 * det
    grad_norm = np.hypot(bx, by)

    valid = (
        base
        & np.isfinite(ux)
        & np.isfinite(uy)
        & np.isfinite(vx)
        & np.isfinite(vy)
        & np.isfinite(bx)
        & np.isfinite(by)
        & np.isfinite(f)
        & (np.abs(f) > 0)
        & np.isfinite(sigma)
        & (sigma > 0)
        & np.isfinite(zeta)
        & (np.sign(zeta) != 0)
        & np.isfinite(disc)
        & (disc > 0)
        & np.isfinite(grad_norm)
        & (grad_norm > 0)
    )

    if not np.any(valid):
        return out

    idx = np.where(valid)
    uxv = ux[idx]
    uyv = uy[idx]
    vxv = vx[idx]
    vyv = vy[idx]
    bxv = bx[idx]
    byv = by[idx]
    fv = f[idx]
    sigv = sigma[idx]
    zetav = zeta[idx]
    deltav = delta[idx]
    discv = disc[idx]
    gradv = grad_norm[idx]

    sqrt_disc = np.sqrt(discv)
    lam_at_minus = 0.5 * (deltav - sqrt_disc)
    lam_at_plus = 0.5 * (deltav + sqrt_disc)

    # A.T eigenvector for lambda_minus: blue direction in the existing project.
    blue_x, blue_y = eigenvector_for_lambda(uxv, vxv, uyv, vyv, lam_at_minus)

    # S compression eigendirection.
    s_a = uxv
    s_b = 0.5 * (uyv + vxv)
    s_d = vyv
    s_trace = s_a + s_d
    s_disc = np.sqrt((s_a - s_d) ** 2 + 4.0 * s_b ** 2)
    lam_s_minus = 0.5 * (s_trace - s_disc)
    lam_s_plus = 0.5 * (s_trace + s_disc)
    comp_x, comp_y = eigenvector_for_lambda(s_a, s_b, s_b, s_d, lam_s_minus)

    grad_x = bxv / gradv
    grad_y = byv / gradv

    finite_vectors = (
        np.isfinite(blue_x)
        & np.isfinite(blue_y)
        & np.isfinite(comp_x)
        & np.isfinite(comp_y)
        & (lam_s_minus < 0)
        & (lam_s_plus > 0)
    )
    theta = np.full_like(uxv, np.nan, dtype=float)
    theta[finite_vectors] = np.rad2deg(np.arccos(np.clip(-zetav[finite_vectors] / sigv[finite_vectors], -1.0, 1.0)))
    theta_acute = np.minimum(theta, 180.0 - theta)
    alpha = acute_axis_angle_deg(blue_x, blue_y, comp_x, comp_y)
    gamma = signed_axis_angle_deg(grad_x, grad_y, comp_x, comp_y, zetav)
    epsilon = signed_axis_angle_deg(grad_x, grad_y, blue_x, blue_y, zetav)

    valid_angle = finite_vectors & np.isfinite(theta) & np.isfinite(alpha) & np.isfinite(gamma) & np.isfinite(epsilon)
    write_idx = (idx[0][valid_angle], idx[1][valid_angle])
    absf = np.abs(fv[valid_angle])
    out["theta_deg"][write_idx] = theta[valid_angle].astype(np.float32)
    out["theta_acute_deg"][write_idx] = theta_acute[valid_angle].astype(np.float32)
    out["alpha_deg"][write_idx] = alpha[valid_angle].astype(np.float32)
    out["gamma_signed_deg"][write_idx] = gamma[valid_angle].astype(np.float32)
    out["epsilon_signed_deg"][write_idx] = epsilon[valid_angle].astype(np.float32)
    out["lambda_at_minus_over_absf"][write_idx] = (lam_at_minus[valid_angle] / absf).astype(np.float32)
    out["lambda_at_plus_over_absf"][write_idx] = (lam_at_plus[valid_angle] / absf).astype(np.float32)
    out["zeta"][write_idx] = zetav[valid_angle].astype(np.float32)
    out["valid_angle"][write_idx] = True
    return out


def classify_modes(angles: dict[str, np.ndarray]) -> np.ndarray:
    theta_acute = angles["theta_acute_deg"]
    alpha = angles["alpha_deg"]
    gamma = angles["gamma_signed_deg"]
    epsilon = angles["epsilon_signed_deg"]
    valid = angles["valid_angle"]
    boundary = theta_acute + alpha
    mode = np.zeros(valid.shape, dtype=np.uint8)

    mode[(epsilon < 0.0) & (gamma < 0.0) & valid] = 1
    mode[(epsilon >= 0.0) & (gamma > 0.0) & (gamma > boundary) & valid] = 2
    mode[(epsilon >= 0.0) & (gamma > epsilon) & (gamma < boundary) & valid] = 3
    mode[(epsilon < 0.0) & (gamma > 0.0) & valid] = 4
    return mode


def mask_data_to_bbox(data: dict[str, np.ndarray], bbox: tuple[float, float, float, float] | None) -> dict[str, np.ndarray]:
    if bbox is None:
        return data
    lon_min, lon_max, lat_min, lat_max = bbox
    region = (
        (data["lon"] >= lon_min)
        & (data["lon"] <= lon_max)
        & (data["lat"] >= lat_min)
        & (data["lat"] <= lat_max)
    )
    masked: dict[str, np.ndarray] = {}
    for key, value in data.items():
        if isinstance(value, np.ndarray) and value.shape == region.shape and key not in {"lon", "lat", "x", "y", "dxc", "dyc", "f"}:
            fill = False if value.dtype == bool else np.nan
            masked[key] = np.where(region, value, fill)
        else:
            masked[key] = value.copy() if isinstance(value, np.ndarray) else value
    masked["valid"] = data["valid"] & region
    return masked


def collect_fl_mask(data: dict[str, np.ndarray], args: argparse.Namespace) -> tuple[np.ndarray, list[Front]]:
    fronts = extract_fronts(
        data,
        args.top_percent,
        args.min_major_axis_km,
        10000,
        args.max_tangent_angle_deg,
    )
    selected = [front for front in fronts if front.length_km > args.stats_min_length_km]
    fl_mask = np.zeros(data["valid"].shape, dtype=bool)
    for front in selected:
        fl_mask[front.rows, front.cols] = True
    return fl_mask, selected


def mode_masks(mode: np.ndarray, valid: np.ndarray) -> dict[str, np.ndarray]:
    masks = {f"Mode{i}": valid & (mode == i) for i in range(1, 5)}
    masks["Unclassified"] = valid & (mode == 0)
    return masks


def write_counts_csv(out_file: Path, rows: list[dict[str, str | int | float]]) -> None:
    out_file.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["scope", "mode", "count", "percent_of_valid"]
    with out_file.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_points_csv(out_file: Path, angles: dict[str, np.ndarray], mode: np.ndarray, mask: np.ndarray) -> None:
    out_file.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "row",
        "col",
        "lon",
        "lat",
        "theta_deg",
        "theta_acute_deg",
        "alpha_deg",
        "epsilon_signed_deg",
        "gamma_signed_deg",
        "lambda_at_minus_over_absf",
        "lambda_at_plus_over_absf",
        "zeta",
        "mode",
    ]
    rows, cols = np.where(mask)
    with out_file.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for r, c in zip(rows, cols):
            m = int(mode[r, c])
            writer.writerow({
                "row": int(r),
                "col": int(c),
                "lon": float(angles["lon"][r, c]),
                "lat": float(angles["lat"][r, c]),
                "theta_deg": float(angles["theta_deg"][r, c]),
                "theta_acute_deg": float(angles["theta_acute_deg"][r, c]),
                "alpha_deg": float(angles["alpha_deg"][r, c]),
                "epsilon_signed_deg": float(angles["epsilon_signed_deg"][r, c]),
                "gamma_signed_deg": float(angles["gamma_signed_deg"][r, c]),
                "lambda_at_minus_over_absf": float(angles["lambda_at_minus_over_absf"][r, c]),
                "lambda_at_plus_over_absf": float(angles["lambda_at_plus_over_absf"][r, c]),
                "zeta": float(angles["zeta"][r, c]),
                "mode": f"Mode{m}" if m else "Unclassified",
            })


def finite_values(angles: dict[str, np.ndarray], mask: np.ndarray) -> dict[str, np.ndarray]:
    return {
        "theta": angles["theta_deg"][mask].astype(float),
        "epsilon": angles["epsilon_signed_deg"][mask].astype(float),
        "gamma": angles["gamma_signed_deg"][mask].astype(float),
        "lm": angles["lambda_at_minus_over_absf"][mask].astype(float),
        "lp": angles["lambda_at_plus_over_absf"][mask].astype(float),
    }


def bin_stat(x: np.ndarray, y: np.ndarray, edges: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    centers = 0.5 * (edges[:-1] + edges[1:])
    mean = np.full(centers.shape, np.nan)
    median = np.full(centers.shape, np.nan)
    count = np.zeros(centers.shape, dtype=int)
    for i, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
        q = (x >= lo) & ((x < hi) | ((i == len(centers) - 1) & (x <= hi)))
        vals = y[q]
        vals = vals[np.isfinite(vals)]
        count[i] = vals.size
        if vals.size:
            mean[i] = float(np.mean(vals))
            median[i] = float(np.median(vals))
    return centers, mean, median, count


def plot_theta_curves(
    angles: dict[str, np.ndarray],
    mode: np.ndarray,
    scope_mask: np.ndarray,
    scope_name: str,
    out_file: Path,
    theta_bin_deg: float,
) -> None:
    valid = angles["valid_angle"] & scope_mask
    masks = mode_masks(mode, valid)
    total = int(np.sum(valid))
    edges = np.arange(0.0, 180.0 + theta_bin_deg, theta_bin_deg)
    fig, axes = plt.subplots(2, 2, figsize=(12, 8.6), sharex=True, sharey=True)
    axes = axes.ravel()
    for ax, name in zip(axes, MODE_NAMES[:4]):
        q = masks[name]
        values = finite_values(angles, q)
        pct = 100.0 * values["theta"].size / total if total else 0.0
        centers, eps_mean, eps_med, cnt = bin_stat(values["theta"], values["epsilon"], edges)
        _, gam_mean, gam_med, _ = bin_stat(values["theta"], values["gamma"], edges)
        ax.axhline(0.0, color="0.55", lw=0.8)
        ax.plot(centers, eps_mean, "-", color="#2364aa", lw=1.8, label="epsilon mean")
        ax.plot(centers, gam_mean, "-", color="#2f8f46", lw=1.8, label="gamma mean")
        ax.plot(centers, eps_med, "--", color="#2364aa", lw=1.1, label="epsilon median")
        ax.plot(centers, gam_med, "--", color="#2f8f46", lw=1.1, label="gamma median")
        ax.set_title(f"{name}: n={values['theta'].size:,}, {pct:.2f}%")
        ax.set_xlim(0, 180)
        ax.set_ylim(-90, 90)
        ax.grid(True, color="0.88", lw=0.6)
        ax.set_xlabel("theta (deg, original definition)")
        ax.set_ylabel("signed acute angle (deg)")
        ax2 = ax.twinx()
        ax2.fill_between(centers, 0, cnt, color="0.82", alpha=0.22, step="mid")
        ax2.set_ylim(0, max(1, int(np.nanmax(cnt)) if cnt.size else 1) * 1.2)
        ax2.set_yticks([])
    axes[0].legend(loc="upper right", fontsize=8)
    fig.suptitle(f"{scope_name}: Mode-conditioned epsilon/gamma statistics by theta", fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out_file.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_file, dpi=300, facecolor="white")
    plt.close(fig)


def conditional_grid(
    x: np.ndarray,
    y: np.ndarray,
    value: np.ndarray,
    x_edges: np.ndarray,
    y_edges: np.ndarray,
    statistic: str,
) -> tuple[np.ndarray, np.ndarray]:
    grid = np.full((len(y_edges) - 1, len(x_edges) - 1), np.nan)
    counts = np.zeros_like(grid, dtype=int)
    xi = np.digitize(x, x_edges) - 1
    yi = np.digitize(y, y_edges) - 1
    ok = (
        np.isfinite(x)
        & np.isfinite(y)
        & np.isfinite(value)
        & (xi >= 0)
        & (xi < len(x_edges) - 1)
        & (yi >= 0)
        & (yi < len(y_edges) - 1)
    )
    for j, i in sorted(set(zip(yi[ok], xi[ok]))):
        vals = value[ok & (xi == i) & (yi == j)]
        counts[j, i] = vals.size
        if vals.size:
            grid[j, i] = float(np.mean(vals) if statistic == "mean" else np.median(vals))
    return grid, counts


def plot_lambda_conditionals(
    angles: dict[str, np.ndarray],
    mode: np.ndarray,
    scope_mask: np.ndarray,
    scope_name: str,
    out_file: Path,
    statistic: str,
    lambda_bin: float,
) -> None:
    valid = angles["valid_angle"] & scope_mask
    masks = mode_masks(mode, valid)
    total = int(np.sum(valid))
    x_edges = np.arange(-5.0, 5.0 + lambda_bin, lambda_bin)
    y_edges = np.arange(-5.0, 5.0 + lambda_bin, lambda_bin)
    variables = [
        ("theta", "theta (deg)", (0, 180), "turbo"),
        ("epsilon", "signed epsilon (deg)", (-90, 90), "RdBu_r"),
        ("gamma", "signed gamma (deg)", (-90, 90), "RdBu_r"),
    ]
    fig, axes = plt.subplots(4, 3, figsize=(12.5, 13.2), sharex=True, sharey=True)
    for row, name in enumerate(MODE_NAMES[:4]):
        values = finite_values(angles, masks[name])
        pct = 100.0 * values["theta"].size / total if total else 0.0
        for col, (key, label, clim, cmap_name) in enumerate(variables):
            ax = axes[row, col]
            grid, counts = conditional_grid(values["lm"], values["lp"], values[key], x_edges, y_edges, statistic)
            cmap = mpl.colormaps[cmap_name].copy()
            cmap.set_bad("#efefef")
            im = ax.pcolormesh(x_edges, y_edges, grid, shading="flat", cmap=cmap, vmin=clim[0], vmax=clim[1])
            ax.axhline(0, color="k", lw=0.7, alpha=0.55)
            ax.axvline(0, color="k", lw=0.7, alpha=0.55)
            ax.set_title(f"{name} {label}\nn={values['theta'].size:,}, {pct:.2f}%")
            ax.set_xlim(x_edges[0], x_edges[-1])
            ax.set_ylim(y_edges[0], y_edges[-1])
            ax.set_aspect("equal", adjustable="box")
            if row == 3:
                ax.set_xlabel("lambda_AT_minus / |f|")
            if col == 0:
                ax.set_ylabel("lambda_AT_plus / |f|")
            cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.025)
            cb.ax.tick_params(labelsize=7)
            if np.any(counts):
                ax.contour(
                    0.5 * (x_edges[:-1] + x_edges[1:]),
                    0.5 * (y_edges[:-1] + y_edges[1:]),
                    counts,
                    levels=[1],
                    colors="0.35",
                    linewidths=0.35,
                )
    fig.suptitle(f"{scope_name}: lambda-JPDF conditional {statistic}", fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.975))
    out_file.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_file, dpi=300, facecolor="white")
    plt.close(fig)


def write_conditional_csv(
    angles: dict[str, np.ndarray],
    mode: np.ndarray,
    scope_mask: np.ndarray,
    out_file: Path,
    lambda_bin: float,
) -> None:
    valid = angles["valid_angle"] & scope_mask
    x_edges = np.arange(-5.0, 5.0 + lambda_bin, lambda_bin)
    y_edges = np.arange(-5.0, 5.0 + lambda_bin, lambda_bin)
    fields = [
        "mode",
        "lambda_at_minus_over_absf_center",
        "lambda_at_plus_over_absf_center",
        "count",
        "theta_mean_deg",
        "theta_median_deg",
        "epsilon_mean_deg",
        "epsilon_median_deg",
        "gamma_mean_deg",
        "gamma_median_deg",
    ]
    out_file.parent.mkdir(parents=True, exist_ok=True)
    with out_file.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for name, q in mode_masks(mode, valid).items():
            if name == "Unclassified":
                continue
            values = finite_values(angles, q)
            xi = np.digitize(values["lm"], x_edges) - 1
            yi = np.digitize(values["lp"], y_edges) - 1
            ok = (
                np.isfinite(values["lm"])
                & np.isfinite(values["lp"])
                & (xi >= 0)
                & (xi < len(x_edges) - 1)
                & (yi >= 0)
                & (yi < len(y_edges) - 1)
            )
            for j, i in sorted(set(zip(yi[ok], xi[ok]))):
                qq = ok & (xi == i) & (yi == j)
                row = {
                    "mode": name,
                    "lambda_at_minus_over_absf_center": float(0.5 * (x_edges[i] + x_edges[i + 1])),
                    "lambda_at_plus_over_absf_center": float(0.5 * (y_edges[j] + y_edges[j + 1])),
                    "count": int(np.sum(qq)),
                }
                for key in ("theta", "epsilon", "gamma"):
                    vals = values[key][qq]
                    row[f"{key}_mean_deg"] = float(np.mean(vals))
                    row[f"{key}_median_deg"] = float(np.median(vals))
                writer.writerow(row)


def load_region_data(args: argparse.Namespace) -> tuple[dict[str, np.ndarray], str]:
    root = args.root
    if args.input_mode == "kuroshio":
        return load_surface_fields(root), "20120101T060000"
    if args.input_mode == "binary":
        return load_binary_region(root, args.use_temperature_gradient)
    if args.input_mode == "westatlantic-nc":
        return load_westatlantic_nc(root / args.nc_name)
    raise ValueError(f"Unsupported input mode: {args.input_mode}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Classify regional full-field and FL points into four angle modes.")
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1] / "KuroshioExt")
    parser.add_argument("--input-mode", choices=("kuroshio", "binary", "westatlantic-nc"), default="kuroshio")
    parser.add_argument("--region-name", default="Kuroshio Extension")
    parser.add_argument("--nc-name", default="LLC4320_pre-SWOT_WestAtlantic_20120913.nc")
    parser.add_argument("--out", type=Path)
    parser.add_argument("--bbox", nargs=4, type=float, metavar=("LON_MIN", "LON_MAX", "LAT_MIN", "LAT_MAX"))
    parser.add_argument("--top-percent", type=float, default=10.0)
    parser.add_argument("--min-major-axis-km", type=float, default=10.0)
    parser.add_argument("--stats-min-length-km", type=float, default=20.0)
    parser.add_argument("--max-tangent-angle-deg", type=float, default=45.0)
    parser.add_argument("--theta-bin-deg", type=float, default=5.0)
    parser.add_argument("--lambda-bin", type=float, default=0.1)
    parser.add_argument("--write-point-csv", action="store_true")
    parser.add_argument("--use-temperature-gradient", action="store_true")
    args = parser.parse_args()

    root = args.root
    data, stamp = load_region_data(args)
    bbox = tuple(args.bbox) if args.bbox is not None else ((155.0, 165.0, 36.0, 40.0) if args.input_mode == "kuroshio" else None)
    out_root = args.out or (root / f"mode_classification_{stamp}")
    out_root.mkdir(parents=True, exist_ok=True)

    angles = compute_full_field_angles(data, bbox)
    angles["lon"] = data["lon"]
    angles["lat"] = data["lat"]
    mode = classify_modes(angles)

    masked_data = mask_data_to_bbox(data, bbox)
    fl_mask, fronts = collect_fl_mask(masked_data, args)
    valid_all = angles["valid_angle"]
    valid_fl = valid_all & fl_mask

    counts: list[dict[str, str | int | float]] = []
    for scope_name, scope_mask in (("AllField", valid_all), ("FL", valid_fl)):
        total = int(np.sum(scope_mask))
        for name, q in mode_masks(mode, scope_mask).items():
            count = int(np.sum(q))
            counts.append({
                "scope": scope_name,
                "mode": name,
                "count": count,
                "percent_of_valid": 100.0 * count / total if total else np.nan,
            })
    write_counts_csv(out_root / "mode_counts_allfield_and_fl.csv", counts)

    np.savez_compressed(
        out_root / "mode_classification_fields.npz",
        lon=data["lon"].astype(np.float32),
        lat=data["lat"].astype(np.float32),
        region=angles["region"],
        valid_angle=angles["valid_angle"],
        fl_mask=fl_mask,
        mode=mode,
        theta_deg=angles["theta_deg"],
        theta_acute_deg=angles["theta_acute_deg"],
        alpha_deg=angles["alpha_deg"],
        epsilon_signed_deg=angles["epsilon_signed_deg"],
        gamma_signed_deg=angles["gamma_signed_deg"],
        lambda_at_minus_over_absf=angles["lambda_at_minus_over_absf"],
        lambda_at_plus_over_absf=angles["lambda_at_plus_over_absf"],
        zeta=angles["zeta"],
    )

    if args.write_point_csv:
        write_points_csv(out_root / "allfield_valid_mode_points.csv", angles, mode, valid_all)
        write_points_csv(out_root / "fl_valid_mode_points.csv", angles, mode, valid_fl)

    plot_theta_curves(angles, mode, valid_all, "All computable field points", out_root / "allfield_mode_theta_epsilon_gamma_curves.png", args.theta_bin_deg)
    plot_theta_curves(angles, mode, valid_fl, "FL grid points", out_root / "fl_mode_theta_epsilon_gamma_curves.png", args.theta_bin_deg)
    for scope_name, scope_mask, stem in (
        ("All computable field points", valid_all, "allfield"),
        ("FL grid points", valid_fl, "fl"),
    ):
        plot_lambda_conditionals(angles, mode, scope_mask, scope_name, out_root / f"{stem}_lambda_jpdf_conditional_mean.png", "mean", args.lambda_bin)
        plot_lambda_conditionals(angles, mode, scope_mask, scope_name, out_root / f"{stem}_lambda_jpdf_conditional_median.png", "median", args.lambda_bin)
        write_conditional_csv(angles, mode, scope_mask, out_root / f"{stem}_lambda_jpdf_conditional_stats.csv", args.lambda_bin)

    with (out_root / "task_parameters.txt").open("w", encoding="utf-8") as f:
        f.write(f"region_name={args.region_name}\n")
        f.write(f"input_mode={args.input_mode}\n")
        f.write(f"stamp={stamp}\n")
        f.write("bbox_lon_lat=" + ("full_domain\n" if bbox is None else f"{bbox[0]},{bbox[1]},{bbox[2]},{bbox[3]}\n"))
        f.write("sign_rule=signed acute epsilon/gamma = -sign(relative_vorticity) * ccw_acute_angle_from_buoyancy_gradient\n")
        f.write("theta_rule=theta keeps original 0-180 definition; theta_acute=min(theta,180-theta) is used only in Mode2/Mode3 boundary\n")
        f.write("mode_partition=mutually exclusive; Mode2/Mode3 require epsilon>=0 to avoid overlap with Mode4\n")
        f.write(f"top_percent={args.top_percent}\n")
        f.write(f"min_major_axis_km={args.min_major_axis_km}\n")
        f.write(f"stats_min_length_km={args.stats_min_length_km}\n")
        f.write(f"max_tangent_angle_deg={args.max_tangent_angle_deg}\n")
        f.write(f"theta_bin_deg={args.theta_bin_deg}\n")
        f.write(f"lambda_bin={args.lambda_bin}\n")
        f.write(f"fl_count={len(fronts)}\n")
        f.write(f"allfield_valid_angle_points={int(np.sum(valid_all))}\n")
        f.write(f"fl_valid_angle_points={int(np.sum(valid_fl))}\n")

    print(f"wrote {out_root}")
    print(f"allfield valid angle points: {int(np.sum(valid_all))}")
    print(f"FL count: {len(fronts)}")
    print(f"FL valid angle points: {int(np.sum(valid_fl))}")


if __name__ == "__main__":
    main()

