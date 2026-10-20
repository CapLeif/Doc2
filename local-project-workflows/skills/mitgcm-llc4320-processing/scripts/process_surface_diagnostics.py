from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage


THETA_COLORS = np.array([
    [100, 0, 100], [175, 0, 175], [220, 0, 220], [50, 50, 200],
    [0, 100, 255], [0, 150, 150], [0, 200, 50], [100, 255, 0],
    [128, 255, 0], [148, 255, 0], [200, 255, 50], [255, 255, 0],
    [255, 200, 0], [255, 160, 0], [255, 125, 0], [225, 25, 0],
], dtype=float) / 255.0
THETA_CMAP = mpl.colors.ListedColormap(THETA_COLORS, name="nclCM353_16")
THETA_CMAP.set_bad("#b3b3b3")
THETA_NORM = mpl.colors.BoundaryNorm(np.linspace(0, 180, 17), 16)


def read_be_f32(path: Path, nx: int, ny: int) -> np.ndarray:
    return np.fromfile(path, dtype=">f4", count=nx * ny).reshape((nx, ny), order="F").astype(np.float32, copy=False)


def locate(root: Path, stem: str) -> Path:
    hits = list(root.rglob(f"{stem}_*"))
    if not hits:
        hits = list(root.rglob(f"*_{stem}_*"))
    if not hits:
        hits = list(root.rglob(stem))
    if len(hits) != 1:
        raise FileNotFoundError(f"Expected one {stem} file under {root}, found {len(hits)}")
    return hits[0]


def locate_optional(root: Path, stem: str) -> Path | None:
    hits = list(root.rglob(f"{stem}_*"))
    if not hits:
        hits = list(root.rglob(f"*_{stem}_*"))
    return hits[0] if len(hits) == 1 else None


def geographic_gradient(q: np.ndarray, lon: np.ndarray, lat: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Gradient of a scalar on a locally curvilinear lon/lat grid."""
    radius = np.float32(6371000.0)
    lon_r = np.deg2rad(lon)
    lat_r = np.deg2rad(lat)
    dlon_i = np.gradient(lon_r, axis=0)
    dlon_j = np.gradient(lon_r, axis=1)
    dlat_i = np.gradient(lat_r, axis=0)
    dlat_j = np.gradient(lat_r, axis=1)
    coslat = np.cos(lat_r)
    dx_i, dx_j = radius * coslat * dlon_i, radius * coslat * dlon_j
    dy_i, dy_j = radius * dlat_i, radius * dlat_j
    dq_i, dq_j = np.gradient(q, axis=0), np.gradient(q, axis=1)
    det = dx_i * dy_j - dx_j * dy_i
    good = np.isfinite(det) & (np.abs(det) > 1e-3)
    qx = np.divide(dq_i * dy_j - dq_j * dy_i, det, out=np.full_like(q, np.nan), where=good)
    qy = np.divide(dx_i * dq_j - dx_j * dq_i, det, out=np.full_like(q, np.nan), where=good)
    return qx, qy


def infer_shape(root: Path) -> tuple[int, int]:
    xc = locate(root, "XC")
    dims = xc.name.split("_")[-1].split("x")
    return int(dims[0]), int(dims[1])


def finite_limits(a: np.ndarray, symmetric: bool = False, q=(1.0, 99.0)) -> tuple[float, float]:
    v = a[np.isfinite(a)]
    lo, hi = np.percentile(v, q)
    if symmetric:
        z = max(abs(lo), abs(hi))
        return -z, z
    return float(lo), float(hi)


def style_axis(ax):
    ax.set_xlabel("Longitude (°E)")
    ax.set_ylabel("Latitude (°N)")
    ax.tick_params(direction="out")


def finite_plot_coordinates(lon, lat):
    """Fill only invalid plotting coordinates from their nearest valid cell."""
    good = np.isfinite(lon) & np.isfinite(lat) & (np.abs(lon) <= 180) & (np.abs(lat) <= 90)
    if good.all():
        return lon, lat
    nearest = ndimage.distance_transform_edt(~good, return_distances=False, return_indices=True)
    return lon[tuple(nearest)], lat[tuple(nearest)]


def surface_plot(lon, lat, field, title, label, cmap, norm, out, figsize=(11, 7.7)):
    fig, ax = plt.subplots(figsize=figsize, dpi=200)
    plot_cmap = mpl.colormaps[cmap].copy() if isinstance(cmap, str) else cmap.copy()
    plot_cmap.set_bad("#b3b3b3")
    pc = ax.pcolormesh(lon, lat, field, shading="nearest", cmap=plot_cmap, norm=norm, rasterized=True)
    style_axis(ax)
    ax.set_title(title)
    cb = fig.colorbar(pc, ax=ax, pad=0.018)
    cb.set_label(label)
    fig.savefig(out, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def median_curve(x, y, bins=36):
    edges = np.linspace(np.nanmin(x), np.nanmax(x), bins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    med = np.full(bins, np.nan)
    ib = np.digitize(x, edges) - 1
    for k in range(bins):
        vals = y[ib == k]
        if vals.size:
            med[k] = np.nanmedian(vals)
    return centers, med


def stats_figure(theta, fields, labels, region, tag, mask, out):
    x = theta[mask]
    fig, axes = plt.subplots(2, 2, figsize=(19.5, 13.5), dpi=200)
    cmaps = ["RdBu_r", "RdBu_r", "magma", "RdBu_r"]
    for ax, field, label, cmap in zip(axes.flat, fields, labels, cmaps):
        y = field[mask]
        hb = ax.hexbin(x, y, gridsize=(72, 55), bins="log", mincnt=1, cmap=cmap, linewidths=0)
        xc, ym = median_curve(x, y)
        ax.plot(xc, ym, "k-", lw=2.2, label="median")
        ax.legend(loc="upper right", frameon=False)
        ax.set_xlabel(r"$\theta$ (degrees)")
        ax.set_ylabel(label)
        ax.set_title(f"{tag}: " + r"$\theta$" + f" vs {label}")
        ax.grid(True, alpha=0.22)
        cb = fig.colorbar(hb, ax=ax, pad=0.015)
        cb.set_label("log10(count)")
    fig.suptitle(f"{region} surface statistics | {tag} | N={x.size:,}", y=0.995)
    fig.savefig(out, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def write_stats(path, theta, fields, mask):
    names = ["lambda_minus_over_f", "Ro", "strain", "divergence"]
    x = theta[mask]
    lines = [f"N={x.size}"]
    for name, a in zip(names, fields):
        lines.append(f"corr(theta,{name})={np.corrcoef(x, a[mask])[0, 1]:.6f}")
    for p in (1, 25, 50, 75, 99):
        lines.append(f"theta_p{p}={np.percentile(x, p):.6f}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("region_dir", type=Path)
    ap.add_argument("--name", required=True)
    ap.add_argument("--output-suffix", default="")
    ap.add_argument("--bbox", nargs=4, type=float, metavar=("LON_MIN", "LON_MAX", "LAT_MIN", "LAT_MAX"))
    args = ap.parse_args()
    root = args.region_dir.resolve()
    nx, ny = infer_shape(root)
    stamp_file = root / "timestamp.txt"
    if stamp_file.exists():
        stamp = stamp_file.read_text(encoding="utf-8").strip()
    else:
        stamp = locate(root, "U").name.split(".", 1)[-1]
    suffix = f"_{args.output_suffix}" if args.output_suffix else ""
    out = root / f"lambda_theta_surface_{stamp}{suffix}"
    out.mkdir(exist_ok=True)

    # All derived diagnostics are collocated at tracer (C) points.
    lon = read_be_f32(locate(root, "XC"), nx, ny)
    lat = read_be_f32(locate(root, "YC"), nx, ny)
    dxc_path = locate_optional(root, "DXC")
    dxc = read_be_f32(dxc_path, nx, ny) if dxc_path is not None else None
    dyc = read_be_f32(locate(root, "DYC"), nx, ny)
    hfac = read_be_f32(locate(root, "hFacC"), nx, ny)
    u_face = read_be_f32(locate(root, "U"), nx, ny)
    v_face = read_be_f32(locate(root, "V"), nx, ny)
    u_face = np.where(np.isfinite(u_face) & (np.abs(u_face) < 20), u_face, np.nan)
    v_face = np.where(np.isfinite(v_face) & (np.abs(v_face) < 20), v_face, np.nan)

    # MITgcm C-grid: U is on west/east faces and V on south/north faces.
    # Average the two faces bounding each tracer cell before combining gradients.
    u = np.full_like(u_face, np.nan)
    v = np.full_like(v_face, np.nan)
    u[:-1, :] = 0.5 * (u_face[:-1, :] + u_face[1:, :])
    v[:, :-1] = 0.5 * (v_face[:, :-1] + v_face[:, 1:])
    base_ok = (hfac > 0) & np.isfinite(lon) & np.isfinite(lat) & (np.abs(lon) <= 180) & (np.abs(lat) <= 90)
    if args.bbox is not None:
        lon_min, lon_max, lat_min, lat_max = args.bbox
        base_ok &= (lon >= lon_min) & (lon <= lon_max) & (lat >= lat_min) & (lat <= lat_max)
    u = np.where(base_ok, u, np.nan)
    v = np.where(base_ok, v, np.nan)

    # Rotate model-grid components to geographic east/north when LLC angle
    # fields are supplied. For unrotated regional tiles AngleCS=1, AngleSN=0.
    angle_cs_path = locate_optional(root, "AngleCS")
    angle_sn_path = locate_optional(root, "AngleSN")
    if angle_cs_path is not None and angle_sn_path is not None:
        angle_cs = read_be_f32(angle_cs_path, nx, ny)
        angle_sn = read_be_f32(angle_sn_path, nx, ny)
        finite_angle = np.isfinite(angle_cs) & np.isfinite(angle_sn) & (np.abs(angle_cs) <= 1.01) & (np.abs(angle_sn) <= 1.01)
        angle_cs = np.where(finite_angle, angle_cs, np.nan)
        angle_sn = np.where(finite_angle, angle_sn, np.nan)
        angle_ok = finite_angle & (np.abs(np.hypot(angle_cs, angle_sn) - 1) < 1e-3)
        u_grid, v_grid = u.copy(), v.copy()
        u = np.where(angle_ok, angle_cs * u_grid - angle_sn * v_grid, np.nan)
        v = np.where(angle_ok, angle_sn * u_grid + angle_cs * v_grid, np.nan)

    if dxc is None:
        lon = np.where(base_ok, lon, np.nan)
        lat = np.where(base_ok, lat, np.nan)
        ux, uy = geographic_gradient(u, lon, lat)
        vx, vy = geographic_gradient(v, lon, lat)
    else:
        ux = np.divide(np.gradient(u, axis=0), dxc, out=np.full_like(u, np.nan), where=dxc > 0)
        uy = np.divide(np.gradient(u, axis=1), dyc, out=np.full_like(u, np.nan), where=dyc > 0)
        vx = np.divide(np.gradient(v, axis=0), dxc, out=np.full_like(v, np.nan), where=dxc > 0)
        vy = np.divide(np.gradient(v, axis=1), dyc, out=np.full_like(v, np.nan), where=dyc > 0)
    zeta = vx - uy
    sigma = np.hypot(ux - vy, vx + uy)
    delta = ux + vy
    f = 2 * 7.2921159e-5 * np.sin(np.deg2rad(lat))
    disc = sigma * sigma - zeta * zeta
    wet = (hfac > 0) & np.isfinite(u) & np.isfinite(v) & np.isfinite(f) & (np.abs(f) > 1e-8)
    # Omega > 0 is required for a real, distinct eigenvalue pair and theta.
    real = wet & (disc > 0)
    theta = np.full((nx, ny), np.nan)
    theta[real & (sigma > 0)] = np.rad2deg(np.arccos(np.clip(-zeta[real & (sigma > 0)] / sigma[real & (sigma > 0)], -1, 1)))
    lambda_minus = np.full((nx, ny), np.nan)
    lambda_minus[real] = 0.5 * (delta[real] - np.sqrt(disc[real]))
    lm_f = lambda_minus / np.abs(f)
    ro = zeta / f
    sigma_f = sigma / np.abs(f)
    div_f = delta / f
    valid = real & np.isfinite(theta) & np.isfinite(lm_f) & np.isfinite(ro) & np.isfinite(sigma_f) & np.isfinite(div_f)
    strong = valid & (lm_f < -1)
    plot_lon, plot_lat = finite_plot_coordinates(lon, lat)

    np.savez_compressed(out / f"{args.name}_surface_diagnostics.npz", lon=lon.astype(">f4"), lat=lat.astype(">f4"),
        lambda_minus_f=lm_f.astype("f4"), theta=theta.astype("f4"), Ro=ro.astype("f4"),
        sigma_f=sigma_f.astype("f4"), div_f=div_f.astype("f4"), valid=valid, land=~wet, sd=strong)
    np.savez_compressed(out / "surface_lambda_theta.npz", lon=lon.astype(">f4"), lat=lat.astype(">f4"),
        theta=theta.astype("f4"), lambda_minus=lambda_minus.astype("f4"), lambda_minus_f=lm_f.astype("f4"), mask=valid, sd=strong)

    name = args.name.replace("_", " ")
    surface_plot(plot_lon, plot_lat, theta, f"{name} surface $\\theta$", r"$\theta$ (degrees)", THETA_CMAP, THETA_NORM, out / "theta_surface.png")
    surface_plot(plot_lon, plot_lat, theta, rf"$\theta$ | {name} surface", r"$\theta$ (degrees)", THETA_CMAP, THETA_NORM, out / "theta.png", (12.5, 8.75))
    specs = [
        (lm_f, r"$\lambda_-/|f|$", r"$\lambda_-/|f|$", "RdBu_r", True, "lambda_minus_over_f"),
        (ro, r"$Ro=\zeta/f$", r"$Ro=\zeta/f$", "RdBu_r", True, "Ro"),
        (sigma_f, r"$\sigma/|f|$", r"$\sigma/|f|$", "magma", False, "strain_rate_over_f"),
        (div_f, r"$\delta/f$", r"$\delta/f$", "RdBu_r", True, "divergence_over_f"),
    ]
    for a, ttl, lab, cmap, sym, fn in specs:
        plot_mask = valid if fn == "lambda_minus_over_f" else wet
        shown = np.where(plot_mask, a, np.nan)
        lim = finite_limits(shown, sym)
        norm = mpl.colors.Normalize(*lim)
        surface_plot(plot_lon, plot_lat, shown, f"{ttl} | {name} surface", lab, cmap, norm, out / f"{fn}.png", (12.5, 8.75))

    suite_stride = max(1, int(np.ceil(max(nx, ny) / 1200)))
    ss = (slice(None, None, suite_stride), slice(None, None, suite_stride))
    fig, axs = plt.subplots(2, 3, figsize=(20, 11.25), dpi=120)
    suite = [(lm_f, specs[0]), (theta, None), (ro, specs[1]), (sigma_f, specs[2]), (div_f, specs[3])]
    for ax, (a, sp) in zip(axs.flat, suite):
        if sp is None:
            pc = ax.pcolormesh(plot_lon[ss], plot_lat[ss], a[ss], shading="nearest", cmap=THETA_CMAP, norm=THETA_NORM, rasterized=True); ttl=r"$\theta$"; lab=r"$\theta$ (degrees)"
        else:
            _, ttl, lab, cmap, sym, fn = sp; plot_mask=valid if fn == "lambda_minus_over_f" else wet; shown=np.where(plot_mask,a,np.nan); norm=mpl.colors.Normalize(*finite_limits(shown,sym)); plot_cmap=mpl.colormaps[cmap].copy(); plot_cmap.set_bad("#b3b3b3"); pc=ax.pcolormesh(plot_lon[ss],plot_lat[ss],shown[ss],shading="nearest",cmap=plot_cmap,norm=norm,rasterized=True)
        style_axis(ax); ax.set_title(f"{ttl} | {name} surface"); fig.colorbar(pc, ax=ax, pad=.015).set_label(lab)
    axs.flat[-1].axis("off")
    fig.savefig(out / f"{args.name}_surface_diagnostic_suite.png", dpi=120, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    fields = [lm_f, ro, sigma_f, div_f]
    labels = [r"$\lambda_-/|f|$", r"$Ro=\zeta/f$", r"$\sigma/|f|$", r"$\delta/f$"]
    stats_figure(theta, fields, labels, name, "all_valid", valid, out / "theta_parameter_statistics_all_valid.png")
    write_stats(out / "theta_parameter_statistics_all_valid.txt", theta, fields, valid)
    if strong.sum() > 0:
        stats_figure(theta, fields, labels, name, "high_negative_lambda", strong, out / "theta_parameter_statistics_high_negative_lambda.png")
        write_stats(out / "theta_parameter_statistics_high_negative_lambda.txt", theta, fields, strong)

    summary = [f"grid={nx}x{ny}", f"lon={np.nanmin(lon):.6f},{np.nanmax(lon):.6f}", f"lat={np.nanmin(lat):.6f},{np.nanmax(lat):.6f}",
        f"valid={valid.sum()} ({100*valid.mean():.2f}%)", f"SD_lambda_minus_over_absf_lt_minus1={strong.sum()} ({100*strong.sum()/max(valid.sum(),1):.3f}% of valid)",
        f"theta_valid_deg={np.percentile(theta[valid],[1,50,99])}", f"lambda_minus_over_absf_valid={np.percentile(lm_f[valid],[1,50,99])}"]
    (out / "summary.txt").write_text("\n".join(summary) + "\n", encoding="utf-8")
    print(out)


if __name__ == "__main__":
    main()

