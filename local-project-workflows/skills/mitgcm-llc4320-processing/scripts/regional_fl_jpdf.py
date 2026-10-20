from __future__ import annotations

import argparse
from pathlib import Path

import netCDF4 as nc
import numpy as np

from kuroshio_fl_theta_gamma_epsilon import (
    ALPHA_T,
    BETA_S,
    GRAVITY,
    OMEGA,
    compute_angles,
    extract_fronts,
    point_rows,
    projected_xy_km,
    read_be_f32,
    select_valid_segment,
    theta_lambda_drawable,
    write_rows_csv,
    write_teg_jpdf,
)
from process_surface_diagnostics import geographic_gradient


def gradient_xy(q: np.ndarray, dx_m: np.ndarray, dy_m: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    return (
        np.divide(np.gradient(q, axis=0), dx_m, out=np.full_like(q, np.nan), where=dx_m > 0),
        np.divide(np.gradient(q, axis=1), dy_m, out=np.full_like(q, np.nan), where=dy_m > 0),
    )


def locate_root_file(root: Path, stem: str) -> Path:
    patterns = (f"{stem}_*", f"*_{stem}_*", stem)
    hits: list[Path] = []
    for pattern in patterns:
        hits = [p for p in root.glob(pattern) if p.is_file()]
        if hits:
            break
    if len(hits) != 1:
        raise FileNotFoundError(f"Expected one root-level {stem} file under {root}, found {len(hits)}")
    return hits[0]


def locate_optional_root_file(root: Path, stem: str) -> Path | None:
    patterns = (f"{stem}_*", f"*_{stem}_*", stem)
    for pattern in patterns:
        hits = [p for p in root.glob(pattern) if p.is_file()]
        if hits:
            return hits[0] if len(hits) == 1 else None
    return None


def infer_root_shape(root: Path) -> tuple[int, int]:
    xc = locate_root_file(root, "XC")
    dims = xc.name.split("_")[-1].split("x")
    return int(dims[0]), int(dims[1])


def load_binary_region(root: Path, use_temperature_gradient: bool) -> tuple[dict[str, np.ndarray], str]:
    nx, ny = infer_root_shape(root)
    stamp_file = root / "timestamp.txt"
    stamp = stamp_file.read_text(encoding="utf-8").strip() if stamp_file.exists() else locate_root_file(root, "U").name.split("_", 2)[0]
    lon = read_be_f32(locate_root_file(root, "XC"), (nx, ny))
    lat = read_be_f32(locate_root_file(root, "YC"), (nx, ny))
    dxc_path = locate_optional_root_file(root, "DXC")
    dxc = read_be_f32(dxc_path, (nx, ny)) if dxc_path is not None else None
    dyc = read_be_f32(locate_root_file(root, "DYC"), (nx, ny))
    hfac = read_be_f32(locate_root_file(root, "hFacC"), (nx, ny))
    theta = read_be_f32(locate_root_file(root, "Theta"), (nx, ny))
    salt_path = locate_optional_root_file(root, "Salt")
    salt = read_be_f32(salt_path, (nx, ny)) if salt_path is not None else None
    u_face = read_be_f32(locate_root_file(root, "U"), (nx, ny))
    v_face = read_be_f32(locate_root_file(root, "V"), (nx, ny))
    u_face = np.where(np.isfinite(u_face) & (np.abs(u_face) < 20), u_face, np.nan)
    v_face = np.where(np.isfinite(v_face) & (np.abs(v_face) < 20), v_face, np.nan)
    u = np.full_like(u_face, np.nan)
    v = np.full_like(v_face, np.nan)
    u[:-1, :] = 0.5 * (u_face[:-1, :] + u_face[1:, :])
    v[:, :-1] = 0.5 * (v_face[:, :-1] + v_face[:, 1:])

    base_ok = (
        (hfac > 0)
        & np.isfinite(lon)
        & np.isfinite(lat)
        & (np.abs(lon) <= 180)
        & (np.abs(lat) <= 90)
        & np.isfinite(theta)
        & np.isfinite(u)
        & np.isfinite(v)
    )
    if salt is not None and not use_temperature_gradient:
        scalar = GRAVITY * (ALPHA_T * theta - BETA_S * salt)
        scalar_note = "buoyancy"
        base_ok &= np.isfinite(salt)
    else:
        scalar = GRAVITY * ALPHA_T * theta
        scalar_note = "temperature_gradient_proxy"

    u = np.where(base_ok, u, np.nan)
    v = np.where(base_ok, v, np.nan)
    scalar = np.where(base_ok, scalar, np.nan)
    if dxc is None:
        bx, by = geographic_gradient(scalar, lon, lat)
        ux, uy = geographic_gradient(u, lon, lat)
        vx, vy = geographic_gradient(v, lon, lat)
        dxc, dyc = metric_from_lonlat(lon, lat)
    else:
        bx, by = gradient_xy(scalar, dxc, dyc)
        ux, uy = gradient_xy(u, dxc, dyc)
        vx, vy = gradient_xy(v, dxc, dyc)
    f = 2 * OMEGA * np.sin(np.deg2rad(lat))
    x, y = projected_xy_km(lon, lat)
    finite = base_ok & np.isfinite(bx) & np.isfinite(by) & np.isfinite(ux) & np.isfinite(uy) & np.isfinite(vx) & np.isfinite(vy)
    return make_data(lon, lat, x, y, dxc, dyc, theta, u, v, bx, by, ux, uy, vx, vy, f, finite), f"{stamp}_{scalar_note}"


def metric_from_lonlat(lon: np.ndarray, lat: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    radius = 6371000.0
    lon_r = np.deg2rad(lon)
    lat_r = np.deg2rad(lat)
    dlon_i = np.gradient(lon_r, axis=0)
    dlon_j = np.gradient(lon_r, axis=1)
    dlat_i = np.gradient(lat_r, axis=0)
    dlat_j = np.gradient(lat_r, axis=1)
    coslat = np.cos(lat_r)
    dx = radius * np.hypot(coslat * dlon_i, dlat_i)
    dy = radius * np.hypot(coslat * dlon_j, dlat_j)
    return dx, dy


def load_westatlantic_nc(path: Path) -> tuple[dict[str, np.ndarray], str]:
    with nc.Dataset(path) as ds:
        lon = np.asarray(ds.variables["XC"][:-1, :-1], dtype=float)
        lat = np.asarray(ds.variables["YC"][:-1, :-1], dtype=float)
        angle_cs = np.asarray(ds.variables["AngleCS"][:-1, :-1], dtype=float)
        angle_sn = np.asarray(ds.variables["AngleSN"][:-1, :-1], dtype=float)
        depth = np.asarray(ds.variables["Depth"][:-1, :-1], dtype=float)
        theta = np.asarray(ds.variables["Theta"][0, 0, :-1, :-1], dtype=float)
        salt = np.asarray(ds.variables["Salt"][0, 0, :-1, :-1], dtype=float)
        u_face = np.asarray(ds.variables["U"][0, 0, :-1, :], dtype=float)
        v_face = np.asarray(ds.variables["V"][0, 0, :, :-1], dtype=float)
    u_grid = 0.5 * (u_face[:, :-1] + u_face[:, 1:])
    v_grid = 0.5 * (v_face[:-1, :] + v_face[1:, :])
    angle_ok = np.isfinite(angle_cs) & np.isfinite(angle_sn) & (np.abs(np.hypot(angle_cs, angle_sn) - 1.0) < 1e-3)
    u = np.where(angle_ok, angle_cs * u_grid - angle_sn * v_grid, np.nan)
    v = np.where(angle_ok, angle_sn * u_grid + angle_cs * v_grid, np.nan)
    base_ok = (
        (depth > 0)
        & np.isfinite(lon)
        & np.isfinite(lat)
        & (np.abs(lon) <= 180)
        & (np.abs(lat) <= 90)
        & np.isfinite(theta)
        & np.isfinite(salt)
        & np.isfinite(u)
        & np.isfinite(v)
        & (np.abs(u) < 20)
        & (np.abs(v) < 20)
    )
    u = np.where(base_ok, u, np.nan)
    v = np.where(base_ok, v, np.nan)
    buoyancy = np.where(base_ok, GRAVITY * (ALPHA_T * theta - BETA_S * salt), np.nan)
    bx, by = geographic_gradient(buoyancy, lon, lat)
    ux, uy = geographic_gradient(u, lon, lat)
    vx, vy = geographic_gradient(v, lon, lat)
    dxc, dyc = metric_from_lonlat(lon, lat)
    f = 2 * OMEGA * np.sin(np.deg2rad(lat))
    x, y = projected_xy_km(lon, lat)
    finite = base_ok & np.isfinite(bx) & np.isfinite(by) & np.isfinite(ux) & np.isfinite(uy) & np.isfinite(vx) & np.isfinite(vy)
    return make_data(lon, lat, x, y, dxc, dyc, theta, u, v, bx, by, ux, uy, vx, vy, f, finite), "20120913_00UTC_buoyancy"


def make_data(
    lon: np.ndarray,
    lat: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    dxc: np.ndarray,
    dyc: np.ndarray,
    theta: np.ndarray,
    u: np.ndarray,
    v: np.ndarray,
    bx: np.ndarray,
    by: np.ndarray,
    ux: np.ndarray,
    uy: np.ndarray,
    vx: np.ndarray,
    vy: np.ndarray,
    f: np.ndarray,
    finite: np.ndarray,
) -> dict[str, np.ndarray]:
    return {
        "lon": lon,
        "lat": lat,
        "x": x,
        "y": y,
        "dxc": dxc,
        "dyc": dyc,
        "theta": np.where(finite, theta, np.nan),
        "u": np.where(finite, u, np.nan),
        "v": np.where(finite, v, np.nan),
        "bx": np.where(finite, bx, np.nan),
        "by": np.where(finite, by, np.nan),
        "ux": np.where(finite, ux, np.nan),
        "uy": np.where(finite, uy, np.nan),
        "vx": np.where(finite, vx, np.nan),
        "vy": np.where(finite, vy, np.nan),
        "M2": np.where(finite, np.hypot(bx, by), np.nan),
        "f": np.where(finite, f, np.nan),
        "valid": finite,
    }


def write_region_jpdf(
    data: dict[str, np.ndarray],
    out_root: Path,
    bin_deg: float,
    top_percent: float,
    min_major_axis_km: float,
    stats_min_length_km: float,
    max_tangent_angle_deg: float,
    stamp: str,
) -> None:
    out_root.mkdir(parents=True, exist_ok=True)
    table_dir = out_root / "tables"
    fronts = extract_fronts(data, top_percent, min_major_axis_km, 10000, max_tangent_angle_deg)
    selected_raw = [front for front in fronts if front.length_km > stats_min_length_km]
    all_rows: list[dict[str, float | int | bool]] = []
    inventory = []
    drawable_count = 0
    for new_id, front in enumerate(selected_raw, start=1):
        full_angles = compute_angles(front, data)
        front.fl_id = new_id
        try:
            segment, original_arc = select_valid_segment(front, full_angles, None)
        except RuntimeError:
            all_rows.extend(point_rows(front, full_angles))
            inventory.append({
                "fl_id": new_id,
                "raw_length_km": front.length_km,
                "length_km": np.nan,
                "raw_original_arc_start_km": np.nan,
                "raw_original_arc_end_km": np.nan,
                "n_points": int(front.rows.size),
                "n_valid_angles": int(np.sum(full_angles["valid"])),
                "n_theta_lambda_drawable": int(np.sum(theta_lambda_drawable(full_angles))),
                "status": "no_continuous_theta_lambda_segment",
            })
            continue
        segment.fl_id = new_id
        angles = compute_angles(segment, data)
        all_rows.extend(point_rows(segment, angles))
        drawable_count += 1
        inventory.append({
            "fl_id": new_id,
            "raw_length_km": front.length_km,
            "length_km": segment.length_km,
            "raw_original_arc_start_km": original_arc[0],
            "raw_original_arc_end_km": original_arc[1],
            "n_points": int(segment.rows.size),
            "n_valid_angles": int(np.sum(angles["valid"])),
            "n_theta_lambda_drawable": int(np.sum(theta_lambda_drawable(angles))),
            "status": "drawable_segment_selected",
        })
    if not all_rows:
        raise RuntimeError("No FL angle rows were produced")
    write_rows_csv(all_rows, table_dir / "all_gt20_FL_theta_gamma_epsilon_points.csv")
    write_rows_csv(inventory, table_dir / "all_gt20_FL_inventory.csv")
    write_teg_jpdf(
        all_rows,
        bin_deg,
        table_dir / f"theta_epsilon_gamma_jpdf_bin{bin_deg:g}deg.csv",
        out_root,
    )
    with (out_root / "task_parameters.txt").open("w", encoding="utf-8") as f:
        f.write(f"stamp={stamp}\n")
        f.write(f"top_percent={top_percent}\n")
        f.write(f"min_major_axis_km={min_major_axis_km}\n")
        f.write(f"stats_min_length_km={stats_min_length_km}\n")
        f.write(f"raw_candidates={len(fronts)}\n")
        f.write(f"selected_raw_gt_length={len(selected_raw)}\n")
        f.write(f"drawable_segments={drawable_count}\n")
        f.write(f"angle_rows={len(all_rows)}\n")
        f.write(f"valid_angle_rows={sum(bool(r['valid']) for r in all_rows)}\n")
        f.write(f"jpdf_bin_deg={bin_deg}\n")
    print(f"{out_root}: candidates={len(fronts)}, selected={len(selected_raw)}, drawable={drawable_count}, rows={len(all_rows)}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("region_dir", type=Path)
    parser.add_argument("--mode", choices=("binary", "westatlantic-nc"), required=True)
    parser.add_argument("--out-name", default="")
    parser.add_argument("--bin-deg", type=float, default=4.0)
    parser.add_argument("--top-percent", type=float, default=10.0)
    parser.add_argument("--min-major-axis-km", type=float, default=10.0)
    parser.add_argument("--stats-min-length-km", type=float, default=20.0)
    parser.add_argument("--max-tangent-angle-deg", type=float, default=45.0)
    parser.add_argument("--use-temperature-gradient", action="store_true")
    args = parser.parse_args()

    root = args.region_dir
    if args.mode == "westatlantic-nc":
        data, stamp = load_westatlantic_nc(root / "LLC4320_pre-SWOT_WestAtlantic_20120913.nc")
    else:
        data, stamp = load_binary_region(root, args.use_temperature_gradient)
    name = args.out_name or f"FL_theta_epsilon_gamma_JPDF_{stamp}"
    write_region_jpdf(
        data,
        root / name,
        args.bin_deg,
        args.top_percent,
        args.min_major_axis_km,
        args.stats_min_length_km,
        args.max_tangent_angle_deg,
        stamp,
    )


if __name__ == "__main__":
    main()

