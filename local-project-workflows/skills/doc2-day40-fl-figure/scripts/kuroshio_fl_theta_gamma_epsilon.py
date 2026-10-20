from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components, dijkstra


NX, NY, NZ = 960, 931, 40
STAMP = "20120101T060000"
OMEGA = 7.2921159e-5
EARTH_RADIUS_KM = 6371.0
GRAVITY = 9.81
ALPHA_T = 2.0e-4
BETA_S = 7.6e-4


@dataclass
class Front:
    fl_id: int
    rows: np.ndarray
    cols: np.ndarray
    x: np.ndarray
    y: np.ndarray
    lon: np.ndarray
    lat: np.ndarray
    arc_km: np.ndarray
    length_km: float


@dataclass(frozen=True)
class ArrowScale:
    global_max: float
    half_length_km: float


def read_be_f32(path: Path, shape: tuple[int, ...]) -> np.ndarray:
    count = int(np.prod(shape))
    data = np.fromfile(path, dtype=">f4", count=count)
    if data.size != count:
        raise ValueError(f"{path} has {data.size} floats, expected {count}")
    return data.reshape(shape, order="F").astype(np.float64, copy=False)


def fill_nearest(a: np.ndarray) -> np.ndarray:
    good = np.isfinite(a)
    if good.all():
        return a
    nearest = ndimage.distance_transform_edt(~good, return_distances=False, return_indices=True)
    return a[tuple(nearest)]


def gradient_xy(q: np.ndarray, dx_m: np.ndarray, dy_m: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    return (
        np.divide(np.gradient(q, axis=0), dx_m, out=np.full_like(q, np.nan), where=dx_m > 0),
        np.divide(np.gradient(q, axis=1), dy_m, out=np.full_like(q, np.nan), where=dy_m > 0),
    )


def projected_xy_km(lon: np.ndarray, lat: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    lon0 = np.nanmin(lon)
    lat0 = np.nanmin(lat)
    lat_ref = np.deg2rad(np.nanmedian(lat))
    x = EARTH_RADIUS_KM * np.cos(lat_ref) * np.deg2rad(lon - lon0)
    y = EARTH_RADIUS_KM * np.deg2rad(lat - lat0)
    return x, y


def load_surface_fields(root: Path) -> dict[str, np.ndarray]:
    lon = read_be_f32(root / "XG_960x931", (NX, NY))
    lat = read_be_f32(root / "YG_960x931", (NX, NY))
    dxc = read_be_f32(root / "DXC_960x931", (NX, NY))
    dyc = read_be_f32(root / "DYC_960x931", (NX, NY))
    hfac = read_be_f32(root / "hFacC_960x931x40", (NX, NY))
    theta = read_be_f32(root / f"Theta_960x931x40.{STAMP}", (NX, NY, NZ))[:, :, 0]
    salt = read_be_f32(root / f"Salt_960x931x40.{STAMP}", (NX, NY, NZ))[:, :, 0]
    u_face = read_be_f32(root / f"U_960x931x40.{STAMP}", (NX, NY, NZ))[:, :, 0]
    v_face = read_be_f32(root / f"V_960x931x40.{STAMP}", (NX, NY, NZ))[:, :, 0]

    base_ok = (
        (hfac > 0)
        & np.isfinite(lon)
        & np.isfinite(lat)
        & (np.abs(lon) <= 180)
        & (np.abs(lat) <= 90)
        & np.isfinite(theta)
        & np.isfinite(salt)
        & np.isfinite(u_face)
        & np.isfinite(v_face)
        & (np.abs(u_face) < 20)
        & (np.abs(v_face) < 20)
    )

    u = np.full_like(u_face, np.nan)
    v = np.full_like(v_face, np.nan)
    u[:-1, :] = 0.5 * (u_face[:-1, :] + u_face[1:, :])
    v[:, :-1] = 0.5 * (v_face[:, :-1] + v_face[:, 1:])
    u = np.where(base_ok, u, np.nan)
    v = np.where(base_ok, v, np.nan)

    buoyancy = GRAVITY * (ALPHA_T * theta - BETA_S * salt)
    bx, by = gradient_xy(np.where(base_ok, buoyancy, np.nan), dxc, dyc)
    ux, uy = gradient_xy(u, dxc, dyc)
    vx, vy = gradient_xy(v, dxc, dyc)
    m2 = np.hypot(bx, by)
    f = 2 * OMEGA * np.sin(np.deg2rad(lat))
    x, y = projected_xy_km(lon, lat)

    finite = base_ok & np.isfinite(bx) & np.isfinite(by) & np.isfinite(ux) & np.isfinite(uy) & np.isfinite(vx) & np.isfinite(vy)
    for key, value in {"theta": theta, "salt": salt, "u": u, "v": v, "buoyancy": buoyancy}.items():
        _ = key
        value[~base_ok] = np.nan
    return {
        "lon": lon,
        "lat": lat,
        "x": x,
        "y": y,
        "dxc": dxc,
        "dyc": dyc,
        "theta": theta,
        "salt": salt,
        "u": u,
        "v": v,
        "bx": np.where(finite, bx, np.nan),
        "by": np.where(finite, by, np.nan),
        "ux": np.where(finite, ux, np.nan),
        "uy": np.where(finite, uy, np.nan),
        "vx": np.where(finite, vx, np.nan),
        "vy": np.where(finite, vy, np.nan),
        "M2": np.where(finite, m2, np.nan),
        "f": np.where(base_ok, f, np.nan),
        "valid": finite,
    }


def component_major_axis_km(rows: np.ndarray, cols: np.ndarray, x: np.ndarray, y: np.ndarray) -> float:
    pts = np.column_stack((x[rows, cols], y[rows, cols]))
    pts = pts[np.all(np.isfinite(pts), axis=1)]
    if pts.shape[0] < 2:
        return 0.0
    centered = pts - pts.mean(axis=0, keepdims=True)
    _, s, _ = np.linalg.svd(centered, full_matrices=False)
    return float(4.0 * s[0] / np.sqrt(max(pts.shape[0], 1)))


def diameter_path(graph, nodes: np.ndarray) -> np.ndarray:
    subgraph = graph[nodes][:, nodes]
    if subgraph.nnz == 0 or nodes.size < 2:
        return nodes
    d0 = dijkstra(subgraph, indices=0)
    finite0 = np.isfinite(d0)
    if np.count_nonzero(finite0) < 2:
        return nodes[:1]
    a = int(np.argmax(np.where(finite0, d0, -np.inf)))
    d1, pred = dijkstra(subgraph, indices=a, return_predecessors=True)
    finite1 = np.isfinite(d1)
    if np.count_nonzero(finite1) < 2:
        return nodes[[a]]
    b = int(np.argmax(np.where(finite1, d1, -np.inf)))
    path = []
    cur = b
    while cur != -9999 and cur != a:
        path.append(cur)
        cur = int(pred[cur])
    path.append(a)
    return nodes[np.array(path[::-1], dtype=int)]


def tangent_aligned_front_paths(
    rows: np.ndarray,
    cols: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    tx: np.ndarray,
    ty: np.ndarray,
    max_tangent_angle_deg: float = 45.0,
) -> list[tuple[np.ndarray, np.ndarray]]:
    n = rows.size
    ncols = x.shape[1]
    node = {int(r) * ncols + int(c): i for i, (r, c) in enumerate(zip(rows, cols))}
    src: list[int] = []
    dst: list[int] = []
    weights: list[float] = []
    min_alignment = float(np.cos(np.deg2rad(max_tangent_angle_deg)))
    for i, (r, c) in enumerate(zip(rows, cols)):
        ti = np.array([tx[r, c], ty[r, c]], dtype=float)
        nti = np.linalg.norm(ti)
        if not np.isfinite(nti) or nti <= 0:
            continue
        ti /= nti
        for dr in (-1, 0, 1):
            for dc in (-1, 0, 1):
                if dr == 0 and dc == 0:
                    continue
                j = node.get(int(r + dr) * ncols + int(c + dc))
                if j is None or j <= i:
                    continue
                tj = np.array([tx[rows[j], cols[j]], ty[rows[j], cols[j]]], dtype=float)
                ntj = np.linalg.norm(tj)
                if not np.isfinite(ntj) or ntj <= 0:
                    continue
                tj /= ntj
                step = np.array([x[rows[j], cols[j]] - x[r, c], y[rows[j], cols[j]] - y[r, c]], dtype=float)
                w = float(np.linalg.norm(step))
                if np.isfinite(w) and w > 0:
                    direction = step / w
                    if min(abs(float(np.dot(direction, ti))), abs(float(np.dot(direction, tj)))) < min_alignment:
                        continue
                    src.extend((i, j))
                    dst.extend((j, i))
                    weights.extend((w, w))
    if not src:
        return []
    graph = coo_matrix((weights, (src, dst)), shape=(n, n)).tocsr()
    degree = np.diff(graph.indptr)
    keep = degree <= 2
    if np.count_nonzero(keep) < 2:
        return []
    kept_nodes = np.nonzero(keep)[0]
    split_graph = graph[kept_nodes][:, kept_nodes]
    n_components, labels = connected_components(split_graph, directed=False)
    paths: list[tuple[np.ndarray, np.ndarray]] = []
    for component in range(n_components):
        local = np.nonzero(labels == component)[0]
        if local.size < 5:
            continue
        path_nodes = diameter_path(split_graph, local)
        original = kept_nodes[path_nodes]
        paths.append((rows[original], cols[original]))
    return paths


def extract_fronts(
    data: dict[str, np.ndarray],
    top_percent: float,
    min_major_axis_km: float,
    count: int,
    max_tangent_angle_deg: float = 45.0,
) -> list[Front]:
    m2 = data["M2"]
    valid = np.isfinite(m2)
    threshold = np.nanpercentile(m2[valid], 100.0 - top_percent)
    fa_mask = valid & (m2 >= threshold)
    structure = np.ones((3, 3), dtype=bool)
    labels, nlab = ndimage.label(fa_mask, structure=structure)

    m2_filled = fill_nearest(m2)
    m2_smooth = ndimage.gaussian_filter(m2_filled, 1.0, mode="nearest")
    dx_km = np.nanmedian(data["dxc"]) / 1000.0
    dy_km = np.nanmedian(data["dyc"]) / 1000.0
    gx, gy = np.gradient(m2_smooth, dx_km, dy_km)
    gxx, gxy1 = np.gradient(gx, dx_km, dy_km)
    gxy2, gyy = np.gradient(gy, dx_km, dy_km)
    gxy = 0.5 * (gxy1 + gxy2)
    gradb = np.hypot(data["bx"], data["by"])
    nx = np.divide(data["bx"], gradb, out=np.full_like(gradb, np.nan), where=gradb > 0)
    ny = np.divide(data["by"], gradb, out=np.full_like(gradb, np.nan), where=gradb > 0)
    tx = -ny
    ty = nx
    normal_derivative = nx * gx + ny * gy
    normal_curvature = nx * nx * gxx + 2 * nx * ny * gxy + ny * ny * gyy
    offset = -normal_derivative / normal_curvature
    ridge_mask = (
        fa_mask
        & np.isfinite(offset)
        & np.isfinite(normal_curvature)
        & (normal_curvature < 0)
        & (np.abs(offset) <= 0.5 * np.hypot(dx_km, dy_km))
    )
    ridge_labels, nridge = ndimage.label(ridge_mask, structure=structure)

    candidates: list[tuple[float, np.ndarray, np.ndarray]] = []
    for lab in range(1, nridge + 1):
        rows, cols = np.nonzero(ridge_labels == lab)
        if rows.size < 8:
            continue
        fa_ids = labels[rows, cols]
        fa_id = int(np.bincount(fa_ids[fa_ids > 0]).argmax()) if np.any(fa_ids > 0) else 0
        fa_rows, fa_cols = np.nonzero(labels == fa_id)
        if component_major_axis_km(fa_rows, fa_cols, data["x"], data["y"]) < min_major_axis_km:
            continue
        if rows.size > 3500:
            keep = np.argsort(m2[rows, cols])[-3500:]
            rows, cols = rows[keep], cols[keep]
        paths = tangent_aligned_front_paths(rows, cols, data["x"], data["y"], tx, ty, max_tangent_angle_deg)
        for path_rows, path_cols in paths:
            if path_rows.size < 2:
                continue
            length = float(np.sum(np.hypot(np.diff(data["x"][path_rows, path_cols]), np.diff(data["y"][path_rows, path_cols]))))
            if length >= 8.0:
                candidates.append((length, path_rows, path_cols))
    candidates.sort(key=lambda item: item[0], reverse=True)

    fronts: list[Front] = []
    for fl_id, (length, rows, cols) in enumerate(candidates[:count], start=1):
        x = data["x"][rows, cols]
        y = data["y"][rows, cols]
        lon = data["lon"][rows, cols]
        lat = data["lat"][rows, cols]
        if y[0] < y[-1] or (y[0] == y[-1] and x[0] > x[-1]):
            rows, cols, x, y, lon, lat = rows[::-1], cols[::-1], x[::-1], y[::-1], lon[::-1], lat[::-1]
        arc = np.r_[0.0, np.cumsum(np.hypot(np.diff(x), np.diff(y)))]
        fronts.append(Front(fl_id, rows, cols, x, y, lon, lat, arc, float(arc[-1])))
    return fronts


def sample_grid(a: np.ndarray, rows: np.ndarray, cols: np.ndarray) -> np.ndarray:
    return ndimage.map_coordinates(a, [rows, cols], order=1, mode="nearest", cval=np.nan)


def principal_axis_angle(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.rad2deg(np.arccos(np.clip(abs(float(np.dot(a, b))), 0.0, 1.0))))


def signed_axis_angle(reference: np.ndarray, axis: np.ndarray, coriolis_sign: float) -> float:
    ccw = np.rad2deg(np.arctan2(reference[0] * axis[1] - reference[1] * axis[0], float(np.dot(reference, axis))))
    ccw_acute = (ccw + 90.0) % 180.0 - 90.0
    return float(-coriolis_sign * ccw_acute)


def normalize(v: np.ndarray) -> np.ndarray:
    n = np.linalg.norm(v)
    return v / n if np.isfinite(n) and n > 0 else np.full_like(v, np.nan)


def compute_angles(front: Front, data: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    fields = {name: sample_grid(data[name], front.rows, front.cols) for name in ("ux", "uy", "vx", "vy", "bx", "by", "f")}
    ux, uy, vx, vy, bx, by, coriolis = (fields[name] for name in ("ux", "uy", "vx", "vy", "bx", "by", "f"))
    n = front.rows.size
    out = {key: np.full(n, np.nan) for key in (
        "theta", "alpha", "beta", "gamma", "epsilon", "gamma_signed", "epsilon_signed",
        "lambda_s_minus", "lambda_s_plus", "lambda_at_minus", "lambda_at_plus",
        "lambda_at_minus_over_f", "lambda_at_plus_over_f", "frontogenetic_rate_over_f",
        "sigma", "zeta", "divergence", "det_lambda",
        "identity_zeta_part", "identity_divergence_part", "identity_det_part", "identity_residual",
    )}
    out["valid"] = np.zeros(n, dtype=bool)
    sigma = np.hypot(ux - vy, uy + vx)
    zeta = vx - uy
    divergence = ux + vy
    det_lambda = ux * vy - uy * vx
    gradb2 = bx * bx + by * by
    rate = -(ux * bx * bx + (uy + vx) * bx * by + vy * by * by) / gradb2
    finite_identity = (
        np.isfinite(sigma)
        & (sigma > 0)
        & np.isfinite(zeta)
        & np.isfinite(divergence)
        & np.isfinite(det_lambda)
    )
    out["sigma"][finite_identity] = sigma[finite_identity]
    out["zeta"][finite_identity] = zeta[finite_identity]
    out["divergence"][finite_identity] = divergence[finite_identity]
    out["det_lambda"][finite_identity] = det_lambda[finite_identity]
    out["identity_zeta_part"][finite_identity] = zeta[finite_identity] ** 2 / sigma[finite_identity] ** 2
    out["identity_divergence_part"][finite_identity] = divergence[finite_identity] ** 2 / sigma[finite_identity] ** 2
    out["identity_det_part"][finite_identity] = -4.0 * det_lambda[finite_identity] / sigma[finite_identity] ** 2
    out["identity_residual"][finite_identity] = (
        out["identity_zeta_part"][finite_identity]
        + out["identity_divergence_part"][finite_identity]
        + out["identity_det_part"][finite_identity]
        - 1.0
    )
    for k in range(n):
        vals = np.array([ux[k], uy[k], vx[k], vy[k], bx[k], by[k], coriolis[k]])
        if not np.all(np.isfinite(vals)) or abs(coriolis[k]) <= 0 or sigma[k] <= 0 or gradb2[k] <= 0:
            continue
        cos_theta = np.clip(-zeta[k] / sigma[k], -1.0, 1.0)
        A = np.array([[ux[k], uy[k]], [vx[k], vy[k]]], dtype=float)
        at_eval, at_vec = np.linalg.eig(A.T)
        if np.max(np.abs(np.imag(at_eval))) > 1e-10:
            continue
        order = np.argsort(np.real(at_eval))
        at_eval = np.real(at_eval[order])
        at_vec = np.real(at_vec[:, order])
        if abs(at_eval[1] - at_eval[0]) <= 100 * np.finfo(float).eps * max(1.0, np.linalg.norm(A)):
            continue
        out["theta"][k] = np.rad2deg(np.arccos(cos_theta))
        out["lambda_at_minus"][k], out["lambda_at_plus"][k] = at_eval
        out["lambda_at_minus_over_f"][k] = at_eval[0] / coriolis[k]
        out["lambda_at_plus_over_f"][k] = at_eval[1] / coriolis[k]
        out["frontogenetic_rate_over_f"][k] = rate[k] / coriolis[k]
        S = 0.5 * (A + A.T)
        s_eval, s_vec = np.linalg.eigh(S)
        order_s = np.argsort(s_eval)
        s_eval = s_eval[order_s]
        s_vec = s_vec[:, order_s]
        if not (s_eval[0] < 0 and s_eval[1] > 0):
            continue
        blue = normalize(at_vec[:, 0])
        red = normalize(at_vec[:, 1])
        comp = normalize(s_vec[:, 0])
        ext = normalize(s_vec[:, 1])
        grad = normalize(np.array([bx[k], by[k]], dtype=float))
        out["alpha"][k] = principal_axis_angle(blue, comp)
        out["beta"][k] = principal_axis_angle(red, ext)
        out["gamma"][k] = principal_axis_angle(comp, grad)
        out["epsilon"][k] = principal_axis_angle(blue, grad)
        out["gamma_signed"][k] = signed_axis_angle(grad, comp, np.sign(coriolis[k]))
        out["epsilon_signed"][k] = signed_axis_angle(grad, blue, np.sign(coriolis[k]))
        out["lambda_s_minus"][k], out["lambda_s_plus"][k] = s_eval
        out["valid"][k] = True
    return out


def theta_lambda_drawable(angles: dict[str, np.ndarray]) -> np.ndarray:
    return (
        np.isfinite(angles["theta"])
        & np.isfinite(angles["lambda_at_minus_over_f"])
        & np.isfinite(angles["lambda_at_plus_over_f"])
        & np.isfinite(angles["frontogenetic_rate_over_f"])
    )


def contiguous_true_runs(mask: np.ndarray) -> list[tuple[int, int]]:
    mask = np.asarray(mask, dtype=bool)
    edges = np.diff(np.r_[False, mask, False].astype(int))
    starts = np.nonzero(edges == 1)[0]
    stops = np.nonzero(edges == -1)[0] - 1
    return list(zip(starts.tolist(), stops.tolist()))


def center_edges(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    if x.size == 1:
        return np.array([x[0] - 0.5, x[0] + 0.5], dtype=float)
    dx = np.diff(x)
    edges = np.empty(x.size + 1, dtype=float)
    edges[1:-1] = 0.5 * (x[:-1] + x[1:])
    edges[0] = x[0] - 0.5 * dx[0]
    edges[-1] = x[-1] + 0.5 * dx[-1]
    return edges


def subset_front(front: Front, start: int, stop: int, fl_id: int | None = None) -> Front:
    rows = front.rows[start : stop + 1].copy()
    cols = front.cols[start : stop + 1].copy()
    x = front.x[start : stop + 1].copy()
    y = front.y[start : stop + 1].copy()
    lon = front.lon[start : stop + 1].copy()
    lat = front.lat[start : stop + 1].copy()
    arc = np.r_[0.0, np.cumsum(np.hypot(np.diff(x), np.diff(y)))]
    return Front(front.fl_id if fl_id is None else fl_id, rows, cols, x, y, lon, lat, arc, float(arc[-1]))


def select_valid_segment(front: Front, angles: dict[str, np.ndarray], arc_window: tuple[float, float] | None) -> tuple[Front, tuple[float, float]]:
    drawable = theta_lambda_drawable(angles)
    if arc_window is not None:
        lo, hi = arc_window
        drawable &= (front.arc_km >= lo) & (front.arc_km <= hi)
    runs = contiguous_true_runs(drawable)
    if not runs:
        raise RuntimeError("No continuous theta/lambda drawable segment was found")
    start, stop = max(runs, key=lambda item: front.arc_km[item[1]] - front.arc_km[item[0]])
    original_arc = (float(front.arc_km[start]), float(front.arc_km[stop]))
    return subset_front(front, start, stop), original_arc


def arrow_samples(front: Front, data: dict[str, np.ndarray], n: int = 9) -> list[dict[str, np.ndarray | float]]:
    targets = np.linspace(0.0, front.length_km, n)
    rows = np.interp(targets, front.arc_km, front.rows)
    cols = np.interp(targets, front.arc_km, front.cols)
    samples = []
    for row, col, arc in zip(rows, cols, targets):
        vals = {name: float(sample_grid(data[name], np.array([row]), np.array([col]))[0]) for name in ("ux", "uy", "vx", "vy")}
        A = np.array([[vals["ux"], vals["uy"]], [vals["vx"], vals["vy"]]], dtype=float)
        if not np.all(np.isfinite(A)):
            continue
        se, sv = np.linalg.eigh(0.5 * (A + A.T))
        so = np.argsort(se)
        se, sv = se[so], sv[:, so]
        ae, av = np.linalg.eig(A.T)
        if np.max(np.abs(np.imag(ae))) <= 1e-10:
            ao = np.argsort(np.real(ae))
            ae, av = np.real(ae[ao]), np.real(av[:, ao])
        else:
            ae, av = np.array([np.nan, np.nan]), np.full((2, 2), np.nan)
        samples.append({
            "row": row,
            "col": col,
            "arc": arc,
            "x": float(sample_grid(data["x"], np.array([row]), np.array([col]))[0]),
            "y": float(sample_grid(data["y"], np.array([row]), np.array([col]))[0]),
            "S_values": se,
            "S_vectors": sv,
            "AT_values": ae,
            "AT_vectors": av,
        })
    return samples


def smoothed_display_line(front: Front, step_km: float = 0.25, sigma_km: float = 2.0) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if front.arc_km.size < 4 or front.length_km <= step_km:
        return front.arc_km, front.x, front.y
    arc = np.arange(0.0, front.length_km, step_km)
    if arc[-1] < front.length_km:
        arc = np.r_[arc, front.length_km]
    x = np.interp(arc, front.arc_km, front.x)
    y = np.interp(arc, front.arc_km, front.y)
    sigma_points = max(sigma_km / step_km, 0.0)
    xs = ndimage.gaussian_filter1d(x, sigma_points, mode="nearest")
    ys = ndimage.gaussian_filter1d(y, sigma_points, mode="nearest")
    xs[0], ys[0] = front.x[0], front.y[0]
    xs[-1], ys[-1] = front.x[-1], front.y[-1]
    return arc, xs, ys


def draw_mode(ax, x0: float, y0: float, v: np.ndarray, lam: float, color, width: float, scale: ArrowScale):
    if not np.all(np.isfinite(v)) or not np.isfinite(lam) or scale.global_max <= 0:
        return
    half = scale.half_length_km * (abs(lam) / scale.global_max) ** 0.5
    uv = half * normalize(v)
    if lam < 0:
        ax.quiver(x0 + uv[0], y0 + uv[1], -uv[0], -uv[1], angles="xy", scale_units="xy", scale=1, color=color, width=0.0028, headwidth=4.5)
        ax.quiver(x0 - uv[0], y0 - uv[1], uv[0], uv[1], angles="xy", scale_units="xy", scale=1, color=color, width=0.0028, headwidth=4.5)
    else:
        ax.quiver(x0, y0, uv[0], uv[1], angles="xy", scale_units="xy", scale=1, color=color, width=0.0028, headwidth=4.5)
        ax.quiver(x0, y0, -uv[0], -uv[1], angles="xy", scale_units="xy", scale=1, color=color, width=0.0028, headwidth=4.5)


def global_arrow_values(fronts: list[Front], data: dict[str, np.ndarray]) -> list[float]:
    values: list[float] = []
    for front in fronts:
        for sample in arrow_samples(front, data):
            values.extend(np.abs(sample["S_values"]).tolist())
            values.extend(np.abs(sample["AT_values"]).tolist())
    return [float(v) for v in values if np.isfinite(v) and v > 0]


def global_arrow_scale(fronts: list[Front], data: dict[str, np.ndarray], lmax_linear: float) -> ArrowScale:
    values = np.asarray(global_arrow_values(fronts, data), dtype=float)
    if values.size == 0:
        return ArrowScale(1.0, lmax_linear)
    return ArrowScale(float(np.nanmax(values)), lmax_linear)


def plot_preview(
    front: Front,
    data: dict[str, np.ndarray],
    angles: dict[str, np.ndarray],
    out_file: Path,
    arrow_scale: ArrowScale,
    original_arc_window: tuple[float, float] | None = None,
):
    samples = arrow_samples(front, data)
    smooth_arc, smooth_x, smooth_y = smoothed_display_line(front)
    colors = [(0.68, 0.68, 0.68), (0.42, 0.42, 0.42), (0.0, 0.25, 0.95), (0.95, 0.02, 0.02)]

    fig = plt.figure(figsize=(15.0, 9.2), dpi=220)
    gs = fig.add_gridspec(2, 2, width_ratios=[1.0, 1.05], height_ratios=[1.08, 0.82], wspace=0.16, hspace=0.48)
    ax0 = fig.add_subplot(gs[:, 0])
    ax1 = fig.add_subplot(gs[0, 1])
    ax2 = fig.add_subplot(gs[1, 1], sharex=ax1)

    margin = 14.0
    q = (
        (data["x"] >= np.nanmin(front.x) - margin)
        & (data["x"] <= np.nanmax(front.x) + margin)
        & (data["y"] >= np.nanmin(front.y) - margin)
        & (data["y"] <= np.nanmax(front.y) + margin)
    )
    m2_plot = np.where(q, data["M2"], np.nan)
    im = ax0.pcolormesh(data["x"], data["y"], m2_plot, shading="nearest", cmap="summer", rasterized=True)
    fig.colorbar(im, ax=ax0, pad=0.018, label="M2")
    ax0.plot(smooth_x, smooth_y, "k-", lw=2.4)
    for s in samples:
        draw_x = float(np.interp(s["arc"], smooth_arc, smooth_x))
        draw_y = float(np.interp(s["arc"], smooth_arc, smooth_y))
        ax0.plot(draw_x, draw_y, "ko", mfc="white", ms=3.5)
        se = np.asarray(s["S_values"])
        sv = np.asarray(s["S_vectors"])
        if np.all(np.isfinite(se)) and se[0] < 0 and se[1] > 0:
            draw_mode(ax0, draw_x, draw_y, sv[:, 0], se[0], colors[0], 1.8, arrow_scale)
            draw_mode(ax0, draw_x, draw_y, sv[:, 1], se[1], colors[1], 1.8, arrow_scale)
        ae = np.asarray(s["AT_values"])
        av = np.asarray(s["AT_vectors"])
        if np.all(np.isfinite(ae)):
            draw_mode(ax0, draw_x, draw_y, av[:, 0], ae[0], colors[2], 0.9, arrow_scale)
            draw_mode(ax0, draw_x, draw_y, av[:, 1], ae[1], colors[3], 1.2, arrow_scale)
    for color, name in zip(colors, ["Vs-", "Vs+", "V-", "V+"]):
        ax0.plot([], [], color=color, lw=2, label=name)
    ax0.legend(loc="lower center", bbox_to_anchor=(0.5, -0.18), ncol=4, frameon=False)
    ax0.set_aspect("equal", adjustable="box")
    ax0.set_xlim(np.nanmin(front.x) - margin, np.nanmax(front.x) + margin)
    ax0.set_ylim(np.nanmin(front.y) - margin, np.nanmax(front.y) + margin)
    ax0.set_xlabel("x (km)")
    ax0.set_ylabel("y (km)")
    ax0.set_title("FL and eigendirections")
    ax0.grid(True, alpha=0.25)

    x = front.arc_km
    edges = center_edges(x)
    widths = np.diff(edges)
    ratios = np.column_stack((
        -angles["lambda_at_plus_over_f"],
        -angles["lambda_at_minus_over_f"],
        angles["frontogenetic_rate_over_f"],
    ))
    max_abs = np.nanmax(np.abs(ratios))
    if not np.isfinite(max_abs) or max_abs <= 0:
        max_abs = 1.0
    mapped = 90.0 + 30.0 * ratios / max_abs
    mapped_lambda = np.clip(mapped[:, :2], 60.0, 120.0)
    ax1.bar(edges[:-1], mapped_lambda[:, 0] - 90.0, width=widths, bottom=90.0, align="edge", color=(1.0, 0.62, 0.62), edgecolor="white", linewidth=0.25, alpha=0.30, label=r"$-\lambda(V+)/f$")
    ax1.bar(edges[:-1], mapped_lambda[:, 1] - 90.0, width=widths, bottom=90.0, align="edge", color=(0.60, 0.72, 1.0), edgecolor="white", linewidth=0.25, alpha=0.30, label=r"$-\lambda(V-)/f$")
    ax1.plot(x, mapped[:, 2], "--", color=(0.25, 0.25, 0.25), lw=1.2, label="frontogenetic rate/f")
    ax1.plot(x, angles["theta"], "-", color=(0.85, 0.10, 0.10), lw=2.0, label=r"$\theta$")
    ax1.plot(x, angles["gamma_signed"], "-", color=(0.10, 0.45, 0.12), lw=2.0, label=r"signed $\gamma$")
    ax1.plot(x, angles["epsilon_signed"], "-", color=(0.05, 0.30, 0.90), lw=1.8, label=r"signed $\epsilon$")
    for yv, ls, col in [(-90, ":", ".72"), (0, "-", ".35"), (60, ":", ".72"), (90, "-", ".58"), (120, ":", ".72")]:
        ax1.axhline(yv, ls=ls, color=col, lw=1.0)
    ax1.set_xlim(0, max(float(front.length_km), 1.0))
    ax1.set_ylim(-90, 180)
    ax1.set_yticks(np.arange(-90, 181, 30))
    ax1.tick_params(labelbottom=False)
    ax1.set_ylabel("Angle (deg): + anticyclonic, - cyclonic")
    ax1.set_title(r"$\theta$, signed $\gamma/\epsilon$, $-\lambda/f$ and frontogenetic rate/f")
    ax1.grid(True, alpha=0.25)
    ax1.legend(loc="upper center", bbox_to_anchor=(0.5, -0.10), ncol=3, frameon=False, fontsize=8, handlelength=2.3)
    ax1r = ax1.twinx()
    ax1r.set_ylim(-90, 180)
    ax1r.set_yticks([60, 90, 120])
    ax1r.set_yticklabels([f"-{max_abs:.3g}", "0", f"+{max_abs:.3g}"])
    ax1r.set_ylabel(r"$-\lambda/f$ and frontogenetic rate/f")
    ax1r.tick_params(colors=".25")
    ax1r.yaxis.label.set_color(".25")

    det_part = angles["identity_det_part"]
    zeta_part = angles["identity_zeta_part"]
    div_part = angles["identity_divergence_part"]
    residual = angles["identity_residual"]
    finite_identity = np.isfinite(det_part) & np.isfinite(zeta_part) & np.isfinite(div_part) & np.isfinite(residual)
    max_residual = float(np.nanmax(np.abs(residual[finite_identity]))) if np.any(finite_identity) else np.nan
    if not np.isfinite(max_residual) or max_residual > 1e-10:
        raise RuntimeError(f"Strain identity residual check failed for FL{front.fl_id:02d}: max_abs={max_residual}")
    det_pos = np.where(finite_identity & (det_part > 0), det_part, 0.0)
    det_neg = np.where(finite_identity & (det_part < 0), det_part, 0.0)
    zeta_base = det_pos
    div_base = det_pos + np.where(finite_identity, zeta_part, 0.0)
    ax2.bar(edges[:-1], det_pos, width=widths, bottom=0.0, align="edge", color=(0.58, 0.65, 0.76), edgecolor="white", linewidth=0.25, alpha=0.78, label=r"$-4\det(\Lambda)/\sigma^2$")
    ax2.bar(edges[:-1], det_neg, width=widths, bottom=0.0, align="edge", color=(0.58, 0.65, 0.76), edgecolor="white", linewidth=0.25, alpha=0.78)
    ax2.bar(edges[:-1], zeta_part, width=widths, bottom=zeta_base, align="edge", color=(0.88, 0.50, 0.38), edgecolor="white", linewidth=0.25, alpha=0.78, label=r"$\zeta^2/\sigma^2$")
    ax2.bar(edges[:-1], div_part, width=widths, bottom=div_base, align="edge", color=(0.42, 0.70, 0.48), edgecolor="white", linewidth=0.25, alpha=0.78, label=r"$\delta^2/\sigma^2$")
    pos_top = np.where(finite_identity, div_base + div_part, np.nan)
    ymin = min(0.0, float(np.nanmin(det_neg[finite_identity])) if np.any(finite_identity) else 0.0)
    ymax = max(1.0, float(np.nanmax(pos_top[finite_identity])) if np.any(finite_identity) else 1.0)
    pad = 0.06 * max(ymax - ymin, 1.0)
    ax2.axhline(0, color=".25", lw=0.9)
    ax2.axhline(1, color=".25", lw=0.9, ls="--")
    ax2.set_xlim(0, max(float(front.length_km), 1.0))
    ax2.set_ylim(ymin - pad, ymax + pad)
    ax2.set_xlabel("Distance from upper FL endpoint (km)")
    ax2.set_ylabel("Identity terms")
    ax2.set_title(r"$\zeta^2/\sigma^2 + \delta^2/\sigma^2 - 4\det(\Lambda)/\sigma^2 = 1$")
    ax2.grid(True, alpha=0.25)
    ax2.legend(loc="lower center", bbox_to_anchor=(0.5, -0.36), ncol=3, frameon=False)

    subtitle = f"Kuroshio Extension {STAMP} FL{front.fl_id:02d} ({front.length_km:.2f} km)"
    if original_arc_window is not None:
        subtitle += f" | original arc {original_arc_window[0]:.1f}-{original_arc_window[1]:.1f} km"
    fig.suptitle(subtitle, fontweight="bold")
    out_file.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_file, dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def write_points_csv(front: Front, angles: dict[str, np.ndarray], out_file: Path):
    out_file.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "fl_id", "grid_index_top_to_bottom", "arc_from_top_km", "x_km", "y_km", "lon", "lat",
        "theta_deg", "alpha_deg", "beta_deg", "gamma_signed_deg", "epsilon_signed_deg",
        "lambda_AT_minus_over_f", "lambda_AT_plus_over_f", "frontogenetic_rate_over_f", "valid",
    ]
    with out_file.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for i in range(front.rows.size):
            writer.writerow({
                "fl_id": front.fl_id,
                "grid_index_top_to_bottom": i + 1,
                "arc_from_top_km": front.arc_km[i],
                "x_km": front.x[i],
                "y_km": front.y[i],
                "lon": front.lon[i],
                "lat": front.lat[i],
                "theta_deg": angles["theta"][i],
                "alpha_deg": angles["alpha"][i],
                "beta_deg": angles["beta"][i],
                "gamma_signed_deg": angles["gamma_signed"][i],
                "epsilon_signed_deg": angles["epsilon_signed"][i],
                "lambda_AT_minus_over_f": angles["lambda_at_minus_over_f"][i],
                "lambda_AT_plus_over_f": angles["lambda_at_plus_over_f"][i],
                "frontogenetic_rate_over_f": angles["frontogenetic_rate_over_f"][i],
                "valid": bool(angles["valid"][i]),
            })


def point_rows(front: Front, angles: dict[str, np.ndarray]) -> list[dict[str, float | int | bool]]:
    rows = []
    for i in range(front.rows.size):
        rows.append({
            "fl_id": int(front.fl_id),
            "grid_index_top_to_bottom": int(i + 1),
            "arc_from_top_km": float(front.arc_km[i]),
            "x_km": float(front.x[i]),
            "y_km": float(front.y[i]),
            "lon": float(front.lon[i]),
            "lat": float(front.lat[i]),
            "theta_deg": float(angles["theta"][i]),
            "alpha_deg": float(angles["alpha"][i]),
            "beta_deg": float(angles["beta"][i]),
            "gamma_deg": float(angles["gamma"][i]),
            "epsilon_deg": float(angles["epsilon"][i]),
            "gamma_signed_deg": float(angles["gamma_signed"][i]),
            "epsilon_signed_deg": float(angles["epsilon_signed"][i]),
            "lambda_AT_minus_over_f": float(angles["lambda_at_minus_over_f"][i]),
            "lambda_AT_plus_over_f": float(angles["lambda_at_plus_over_f"][i]),
            "frontogenetic_rate_over_f": float(angles["frontogenetic_rate_over_f"][i]),
            "valid": bool(angles["valid"][i]),
            "theta_lambda_drawable": bool(theta_lambda_drawable(angles)[i]),
        })
    return rows


def write_rows_csv(rows: list[dict[str, float | int | bool]], out_file: Path) -> None:
    out_file.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        return
    with out_file.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def finite_stats(values: np.ndarray) -> tuple[float, float, float, float, float]:
    values = values[np.isfinite(values)]
    if values.size == 0:
        return (np.nan, np.nan, np.nan, np.nan, np.nan)
    return (
        float(np.mean(values)),
        float(np.std(values, ddof=1)) if values.size > 1 else 0.0,
        float(np.median(values)),
        float(np.percentile(values, 25)),
        float(np.percentile(values, 75)),
    )


def write_summary_csv(rows: list[dict[str, float | int | bool]], out_file: Path) -> None:
    valid_rows = [r for r in rows if bool(r["valid"])]
    ids = sorted({int(r["fl_id"]) for r in rows})
    summary = []
    for label, subset in [(f"FL{i:02d}", [r for r in valid_rows if int(r["fl_id"]) == i]) for i in ids] + [("ALL", valid_rows)]:
        entry: dict[str, float | int | str] = {"scope": label, "count": len(subset)}
        for name in ("theta_deg", "gamma_signed_deg", "epsilon_signed_deg", "gamma_deg", "epsilon_deg"):
            vals = np.asarray([float(r[name]) for r in subset], dtype=float)
            mean, std, med, p25, p75 = finite_stats(vals)
            stem = name.replace("_deg", "")
            entry[f"{stem}_mean_deg"] = mean
            entry[f"{stem}_std_deg"] = std
            entry[f"{stem}_median_deg"] = med
            entry[f"{stem}_p25_deg"] = p25
            entry[f"{stem}_p75_deg"] = p75
        summary.append(entry)
    write_rows_csv(summary, out_file)


def write_theta_bin_stats(rows: list[dict[str, float | int | bool]], out_file: Path, fig_file: Path) -> None:
    valid = [r for r in rows if bool(r["valid"])]
    theta = np.asarray([float(r["theta_deg"]) for r in valid], dtype=float)
    gamma = np.asarray([float(r["gamma_signed_deg"]) for r in valid], dtype=float)
    epsilon = np.asarray([float(r["epsilon_signed_deg"]) for r in valid], dtype=float)
    edges = np.arange(0.0, 180.0 + 10.0, 10.0)
    out_rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        q = (theta >= lo) & ((theta < hi) | ((hi == edges[-1]) & (theta <= hi)))
        gm, gs, gmed, gp25, gp75 = finite_stats(gamma[q])
        em, es, emed, ep25, ep75 = finite_stats(epsilon[q])
        out_rows.append({
            "theta_bin_left_deg": lo,
            "theta_bin_right_deg": hi,
            "theta_bin_center_deg": 0.5 * (lo + hi),
            "count": int(np.count_nonzero(q)),
            "gamma_signed_mean_deg": gm,
            "gamma_signed_std_deg": gs,
            "gamma_signed_median_deg": gmed,
            "gamma_signed_p25_deg": gp25,
            "gamma_signed_p75_deg": gp75,
            "epsilon_signed_mean_deg": em,
            "epsilon_signed_std_deg": es,
            "epsilon_signed_median_deg": emed,
            "epsilon_signed_p25_deg": ep25,
            "epsilon_signed_p75_deg": ep75,
        })
    write_rows_csv(out_rows, out_file)
    centers = np.asarray([r["theta_bin_center_deg"] for r in out_rows], dtype=float)
    fig, axes = plt.subplots(2, 1, figsize=(8.5, 6.8), dpi=180, sharex=True)
    axes[0].errorbar(centers, [r["gamma_signed_mean_deg"] for r in out_rows], [r["gamma_signed_std_deg"] for r in out_rows], fmt="-o", color=(0.10, 0.45, 0.12), label="signed gamma")
    axes[0].errorbar(centers, [r["epsilon_signed_mean_deg"] for r in out_rows], [r["epsilon_signed_std_deg"] for r in out_rows], fmt="-s", color=(0.05, 0.30, 0.90), label="signed epsilon")
    axes[0].set_ylabel("Mean +/- SD (deg)")
    axes[0].set_ylim(-90, 90)
    axes[0].grid(True, alpha=0.25)
    axes[0].legend(frameon=False)
    axes[1].plot(centers, [r["gamma_signed_median_deg"] for r in out_rows], "-o", color=(0.10, 0.45, 0.12), label="signed gamma median")
    axes[1].plot(centers, [r["epsilon_signed_median_deg"] for r in out_rows], "-s", color=(0.05, 0.30, 0.90), label="signed epsilon median")
    axes[1].set_xlabel("theta (deg)")
    axes[1].set_ylabel("Median (deg)")
    axes[1].set_xlim(0, 180)
    axes[1].set_ylim(-90, 90)
    axes[1].grid(True, alpha=0.25)
    axes[1].legend(frameon=False)
    fig.suptitle("All >20 km NewFL points: signed epsilon/gamma conditioned on theta", fontweight="bold")
    fig.tight_layout()
    fig_file.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(fig_file, dpi=300, facecolor="white")
    plt.close(fig)


def write_pairwise_jpdf(
    centers_a: np.ndarray,
    centers_b: np.ndarray,
    counts: np.ndarray,
    bin_deg: float,
    x_name: str,
    y_name: str,
    x_label: str,
    y_label: str,
    x_lim: tuple[float, float],
    y_lim: tuple[float, float],
    out_csv: Path,
    out_fig: Path,
) -> None:
    keys = np.column_stack((centers_a, centers_b))
    unique, inverse = np.unique(keys, axis=0, return_inverse=True)
    pair_counts = np.zeros(unique.shape[0], dtype=np.int64)
    np.add.at(pair_counts, inverse, counts.astype(np.int64, copy=False))
    total = int(pair_counts.sum())
    out_rows = []
    for (a_center, b_center), count in zip(unique, pair_counts):
        prob = float(count / total) if total else np.nan
        out_rows.append({
            f"{x_name}_center_deg": float(a_center),
            f"{y_name}_center_deg": float(b_center),
            "count": int(count),
            "probability": prob,
            "density_per_deg2": prob / (bin_deg ** 2) if total else np.nan,
        })
    write_rows_csv(out_rows, out_csv)
    if not out_rows:
        return

    points = unique.astype(float, copy=False)
    weights = pair_counts.astype(float)
    pca_rows = []
    pca = None
    if points.shape[0] >= 2 and np.sum(weights) > 0:
        mean = np.average(points, axis=0, weights=weights)
        centered = points - mean
        cov = (centered * weights[:, None]).T @ centered / np.sum(weights)
        evals, evecs = np.linalg.eigh(cov)
        order = np.argsort(evals)[::-1]
        evals = evals[order]
        evecs = evecs[:, order]
        total_var = float(np.sum(evals))
        pca = (mean, evals, evecs, total_var)
        for idx in range(2):
            ratio = float(evals[idx] / total_var) if total_var > 0 else np.nan
            pca_rows.append({
                "component": f"PC{idx + 1}",
                f"{x_name}_mean_deg": float(mean[0]),
                f"{y_name}_mean_deg": float(mean[1]),
                "eigenvalue_deg2": float(evals[idx]),
                "std_deg": float(np.sqrt(max(evals[idx], 0.0))),
                "explained_variance_ratio": ratio,
                f"{x_name}_vector": float(evecs[0, idx]),
                f"{y_name}_vector": float(evecs[1, idx]),
            })
    if pca_rows:
        write_rows_csv(pca_rows, out_csv.with_name(f"{out_csv.stem}_pca.csv"))

    x_edges = np.arange(x_lim[0], x_lim[1] + bin_deg, bin_deg)
    y_edges = np.arange(y_lim[0], y_lim[1] + bin_deg, bin_deg)
    grid = np.full((y_edges.size - 1, x_edges.size - 1), np.nan, dtype=float)
    xi = np.floor((unique[:, 0] - x_lim[0]) / bin_deg).astype(int)
    yi = np.floor((unique[:, 1] - y_lim[0]) / bin_deg).astype(int)
    keep = (xi >= 0) & (xi < grid.shape[1]) & (yi >= 0) & (yi < grid.shape[0])
    probability = pair_counts.astype(float) / total if total else np.full_like(pair_counts, np.nan, dtype=float)
    grid[yi[keep], xi[keep]] = np.log10(np.maximum(probability[keep], np.finfo(float).tiny))

    cmap = mpl.colormaps["viridis"].copy()
    cmap.set_bad("white")
    fig, ax = plt.subplots(figsize=(8.0, 6.2), dpi=180)
    im = ax.pcolormesh(x_edges, y_edges, grid, shading="flat", cmap=cmap)
    ax.set_xlim(*x_lim)
    ax.set_ylim(*y_lim)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.grid(True, color="0.85", lw=0.5, alpha=0.55)
    x_ref = 0.0 if x_lim[0] < 0.0 < x_lim[1] else 90.0
    ax.axvline(x_ref, color="black", lw=1.15, alpha=0.9, zorder=4)
    ax.axhline(0.0, color="black", lw=1.15, alpha=0.9, zorder=4)
    if pca is not None:
        mean, evals, evecs, total_var = pca
        ax.plot(mean[0], mean[1], "ro", ms=5.2, mec="black", mew=0.8, zorder=5, label="PCA mean")
        for idx, (style, label) in enumerate((("-", "PC1"), ("--", "PC2"))):
            std = float(np.sqrt(max(evals[idx], 0.0)))
            if not np.isfinite(std) or std <= 0:
                continue
            vec = evecs[:, idx]
            t = 2.0 * std
            if abs(vec[0]) > 1e-12:
                t = min(t, abs((x_lim[0] - mean[0]) / vec[0]), abs((x_lim[1] - mean[0]) / vec[0]))
            if abs(vec[1]) > 1e-12:
                t = min(t, abs((y_lim[0] - mean[1]) / vec[1]), abs((y_lim[1] - mean[1]) / vec[1]))
            start = mean - t * vec
            stop = mean + t * vec
            ratio = float(evals[idx] / total_var) if total_var > 0 else np.nan
            ax.plot([start[0], stop[0]], [start[1], stop[1]], style, color="red", lw=1.9, zorder=6, label=f"{label} ({100 * ratio:.1f}%)")
        ax.legend(loc="upper right", frameon=True, framealpha=0.86, borderpad=0.4, fontsize=8)
    ax.set_title(f"{x_label} vs {y_label} JPDF, bin={bin_deg:g} deg")
    fig.colorbar(im, ax=ax, pad=0.018, label="log10 probability")
    fig.tight_layout()
    out_fig.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_fig, dpi=300, facecolor="white")
    plt.close(fig)


def write_teg_jpdf(rows: list[dict[str, float | int | bool]], bin_deg: float, out_csv: Path, out_dir: Path) -> None:
    valid = [r for r in rows if bool(r["valid"])]
    theta = np.asarray([float(r["theta_deg"]) for r in valid], dtype=float)
    epsilon = np.asarray([float(r["epsilon_signed_deg"]) for r in valid], dtype=float)
    gamma = np.asarray([float(r["gamma_signed_deg"]) for r in valid], dtype=float)
    q = np.isfinite(theta) & np.isfinite(epsilon) & np.isfinite(gamma)
    theta, epsilon, gamma = theta[q], epsilon[q], gamma[q]
    t_bin = np.floor(theta / bin_deg).astype(int)
    e_bin = np.floor((epsilon + 90.0) / bin_deg).astype(int)
    g_bin = np.floor((gamma + 90.0) / bin_deg).astype(int)
    keep = (t_bin >= 0) & (t_bin < int(round(180 / bin_deg))) & (e_bin >= 0) & (e_bin < int(round(180 / bin_deg))) & (g_bin >= 0) & (g_bin < int(round(180 / bin_deg)))
    keys = np.column_stack((t_bin[keep], e_bin[keep], g_bin[keep]))
    unique, counts = np.unique(keys, axis=0, return_counts=True)
    total = int(np.sum(counts))
    out_rows = []
    for (tb, eb, gb), count in zip(unique, counts):
        prob = float(count / total) if total else np.nan
        out_rows.append({
            "theta_center_deg": (int(tb) + 0.5) * bin_deg,
            "epsilon_signed_center_deg": -90.0 + (int(eb) + 0.5) * bin_deg,
            "gamma_signed_center_deg": -90.0 + (int(gb) + 0.5) * bin_deg,
            "count": int(count),
            "probability": prob,
            "density_per_deg3": prob / (bin_deg ** 3) if total else np.nan,
        })
    write_rows_csv(out_rows, out_csv)
    if not out_rows:
        return
    theta_centers = np.asarray([r["theta_center_deg"] for r in out_rows], dtype=float)
    epsilon_centers = np.asarray([r["epsilon_signed_center_deg"] for r in out_rows], dtype=float)
    gamma_centers = np.asarray([r["gamma_signed_center_deg"] for r in out_rows], dtype=float)
    counts = np.asarray([r["count"] for r in out_rows], dtype=np.int64)
    suffix = f"bin{bin_deg:g}deg"
    write_pairwise_jpdf(
        theta_centers,
        epsilon_centers,
        counts,
        bin_deg,
        "theta",
        "epsilon_signed",
        "theta (deg)",
        "signed epsilon (deg)",
        (0.0, 180.0),
        (-90.0, 90.0),
        out_csv.parent / f"theta_epsilon_jpdf_{suffix}.csv",
        out_dir / f"theta_epsilon_jpdf_{suffix}.png",
    )
    write_pairwise_jpdf(
        theta_centers,
        gamma_centers,
        counts,
        bin_deg,
        "theta",
        "gamma_signed",
        "theta (deg)",
        "signed gamma (deg)",
        (0.0, 180.0),
        (-90.0, 90.0),
        out_csv.parent / f"theta_gamma_jpdf_{suffix}.csv",
        out_dir / f"theta_gamma_jpdf_{suffix}.png",
    )
    write_pairwise_jpdf(
        epsilon_centers,
        gamma_centers,
        counts,
        bin_deg,
        "epsilon_signed",
        "gamma_signed",
        "signed epsilon (deg)",
        "signed gamma (deg)",
        (-90.0, 90.0),
        (-90.0, 90.0),
        out_csv.parent / f"epsilon_gamma_jpdf_{suffix}.csv",
        out_dir / f"epsilon_gamma_jpdf_{suffix}.png",
    )


def run_batch(args: argparse.Namespace) -> None:
    root = args.root.resolve()
    out_root = args.out or (root / f"FL_lambdaJPDF_{STAMP}" / "batch_NewFL")
    data = load_surface_fields(root)
    fronts = extract_fronts(data, args.top_percent, args.min_major_axis_km, 10000, args.max_tangent_angle_deg)
    selected_raw = [front for front in fronts if front.length_km > args.stats_min_length_km]
    selected_raw = selected_raw[: args.stats_count] if args.stats_count > 0 else selected_raw
    if not selected_raw:
        raise RuntimeError("No FL candidates passed the batch length gate")
    prepared: list[tuple[Front, dict[str, np.ndarray], tuple[float, float]]] = []
    all_rows: list[dict[str, float | int | bool]] = []
    inventory = []
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
                "plotted": False,
                "status": "no_continuous_theta_lambda_segment",
            })
            continue
        segment.fl_id = new_id
        angles = compute_angles(segment, data)
        prepared.append((segment, angles, original_arc))
        all_rows.extend(point_rows(segment, angles))
        inventory.append({
            "fl_id": new_id,
            "raw_length_km": front.length_km,
            "length_km": segment.length_km,
            "raw_original_arc_start_km": original_arc[0],
            "raw_original_arc_end_km": original_arc[1],
            "n_points": int(segment.rows.size),
            "n_valid_angles": int(np.sum(angles["valid"])),
            "n_theta_lambda_drawable": int(np.sum(theta_lambda_drawable(angles))),
            "plotted": False,
            "status": "drawable_segment_selected",
        })
    if not prepared:
        raise RuntimeError("No drawable FL segment was found in the batch candidates")
    arrow_scale = global_arrow_scale([item[0] for item in prepared[: args.plot_count]], data, args.arrow_half_length_km)
    plot_ids = set(args.plot_ids or [])
    figure_dir = out_root / "combined_figures"
    written_figures = 0
    for front, angles, original_arc in prepared[: args.plot_count]:
        if plot_ids and front.fl_id not in plot_ids:
            continue
        png = figure_dir / f"FL{front.fl_id:02d}_arrows_theta_gamma_epsilon.png"
        plot_preview(front, data, angles, png, arrow_scale, original_arc)
        written_figures += 1
        for item in inventory:
            if item["fl_id"] == front.fl_id:
                item["plotted"] = True
                break
        print(f"figure FL{front.fl_id:02d}: {front.length_km:.3f} km -> {png}")
    table_dir = out_root / "tables"
    write_rows_csv(all_rows, table_dir / "all_gt20_FL_theta_gamma_epsilon_points.csv")
    write_summary_csv(all_rows, table_dir / "all_gt20_FL_theta_gamma_epsilon_summary.csv")
    write_theta_bin_stats(
        all_rows,
        table_dir / "all_gt20_FL_theta_binned_signed_gamma_epsilon.csv",
        out_root / "all_gt20_FL_theta_gamma_epsilon_statistics.png",
    )
    write_teg_jpdf(
        all_rows,
        args.jpdf_bin_deg,
        table_dir / f"theta_epsilon_gamma_jpdf_bin{args.jpdf_bin_deg:g}deg.csv",
        out_root,
    )
    write_rows_csv(inventory, table_dir / "all_gt20_FL_inventory.csv")
    with (out_root / "task_parameters.txt").open("w", encoding="utf-8") as f:
        f.write(f"stamp={STAMP}\n")
        f.write(f"top_percent={args.top_percent}\n")
        f.write(f"min_major_axis_km={args.min_major_axis_km}\n")
        f.write(f"stats_min_length_km={args.stats_min_length_km}\n")
        f.write(f"raw_candidates={len(fronts)}\n")
        f.write(f"selected_raw_gt_length={len(selected_raw)}\n")
        f.write(f"drawable_segments={len(prepared)}\n")
        f.write(f"plot_count={min(args.plot_count, len(prepared))}\n")
        f.write(f"plot_ids={','.join(str(i) for i in sorted(plot_ids)) if plot_ids else 'all'}\n")
        f.write(f"figures_written_this_run={written_figures}\n")
        f.write(f"jpdf_bin_deg={args.jpdf_bin_deg}\n")
        f.write(f"arrow_global_max={arrow_scale.global_max:.12g}\n")
        f.write(f"arrow_direct_power_half_length_km={arrow_scale.half_length_km:.12g}\n")
    print(f"Batch complete: {len(selected_raw)} FL candidates in tables, {len(prepared)} drawable segments, {written_figures} figures written, output={out_root}")


def select_target_front(fronts: list[Front], target_bbox: list[float] | None) -> list[Front]:
    if target_bbox is None:
        return fronts
    xmin, xmax, ymin, ymax = target_bbox
    hits = []
    for front in fronts:
        cx = float(np.nanmean(front.x))
        cy = float(np.nanmean(front.y))
        if xmin <= cx <= xmax and ymin <= cy <= ymax:
            hits.append(front)
    if not hits:
        raise RuntimeError(f"No FL candidate centroid falls inside target bbox {target_bbox}")
    hits.sort(key=lambda item: item.length_km, reverse=True)
    chosen = hits[0]
    chosen.fl_id = 1
    return [chosen] + [front for front in fronts if front is not chosen]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1] / "KuroshioExt")
    parser.add_argument("--out", type=Path, default=None)
    parser.add_argument("--top-percent", type=float, default=10.0)
    parser.add_argument("--min-major-axis-km", type=float, default=10.0)
    parser.add_argument("--preview-count", type=int, default=1)
    parser.add_argument("--global-arrow-count", type=int, default=30)
    parser.add_argument("--arrow-half-length-km", type=float, default=10.0)
    parser.add_argument("--max-tangent-angle-deg", type=float, default=45.0)
    parser.add_argument("--arc-window", nargs=2, type=float, metavar=("ARC_MIN_KM", "ARC_MAX_KM"))
    parser.add_argument("--target-bbox", nargs=4, type=float, metavar=("XMIN_KM", "XMAX_KM", "YMIN_KM", "YMAX_KM"))
    parser.add_argument("--batch", action="store_true")
    parser.add_argument("--plot-count", type=int, default=50)
    parser.add_argument("--plot-ids", nargs="+", type=int)
    parser.add_argument("--stats-min-length-km", type=float, default=20.0)
    parser.add_argument("--stats-count", type=int, default=0)
    parser.add_argument("--jpdf-bin-deg", type=float, default=0.5)
    args = parser.parse_args()

    if args.batch:
        run_batch(args)
        return

    root = args.root.resolve()
    out_root = args.out or (root / f"FL_lambdaJPDF_{STAMP}" / "preview")
    data = load_surface_fields(root)
    extract_count = max(args.preview_count, args.global_arrow_count, 300 if args.target_bbox is not None else 0)
    fronts = extract_fronts(data, args.top_percent, args.min_major_axis_km, extract_count, args.max_tangent_angle_deg)
    if not fronts:
        raise RuntimeError("No FL candidates were extracted")
    fronts = select_target_front(fronts, args.target_bbox)
    arrow_scale = global_arrow_scale(fronts[: args.global_arrow_count], data, args.arrow_half_length_km)
    for front in fronts[: args.preview_count]:
        full_angles = compute_angles(front, data)
        original_arc_window = None
        if args.arc_window is not None:
            front, original_arc_window = select_valid_segment(front, full_angles, tuple(args.arc_window))
        else:
            front, original_arc_window = select_valid_segment(front, full_angles, None)
        angles = compute_angles(front, data)
        png = out_root / f"FL{front.fl_id:02d}_arrows_theta_gamma_epsilon.png"
        plot_preview(front, data, angles, png, arrow_scale, original_arc_window)
        write_points_csv(front, angles, out_root / f"FL{front.fl_id:02d}_theta_gamma_epsilon_points.csv")
        drawable = int(np.sum(theta_lambda_drawable(angles)))
        print(
            f"FL{front.fl_id:02d}: length_km={front.length_km:.3f}, "
            f"theta_lambda_drawable={drawable}/{front.rows.size}, "
            f"valid_angles={int(np.sum(angles['valid']))}/{front.rows.size}, "
            f"global_arrow_max={arrow_scale.global_max:.6g}, figure={png}"
        )


if __name__ == "__main__":
    main()

