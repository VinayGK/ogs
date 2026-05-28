#!/usr/bin/env python3
"""
audit_extract_ms33.py — auditable post-processing for ANCHORS MS33 model outputs.

Reads a PVD (MFront-bridge or native dsm_micromacro outputs), extracts per-timestep
scalars and final-state field snapshots according to a per-model YAML probe spec,
and writes CSV history + summary JSON + PDF/PNG figures. All scientific anchors
(sampling locations, MaterialID zones, time-series fields, marker times) live in
probes_model<X>.yaml — reviewable as data, not buried in code.

Usage
-----
  python audit_extract_ms33.py --model {I,III,IV,VII} \\
      --pvd /path/to/<model>.pvd [--pvd /path/to/another.pvd ...] \\
      --tag {mfront|native} \\
      --probes probes_model<X>.yaml \\
      -o <outdir>

  # Model I dd-sweep (three PVDs, one per dry density):
  python audit_extract_ms33.py --model I \\
      --pvd .../ms33_modelI_dd1400.pvd --rho-dry 1400 \\
      --pvd .../ms33_modelI_dd1600.pvd --rho-dry 1600 \\
      --pvd .../ms33_modelI_dd1800.pvd --rho-dry 1800 \\
      --tag mfront --probes probes_modelI.yaml -o out/

Stress convention
-----------------
The CSV column `mean_stress_MPa` = -tr(sigma_VTU)/3·1e-6 (compression positive).
Under BishopsSaturationCutoff with cutoff_value=1 (set in all four MS33 PRJs),
χ ≡ 0 on the unsaturated branch and p_L = 0 at saturation, so the VTU `sigma`
field is numerically identical to σ_total throughout the simulation. Rename the
column if the cutoff convention changes (cf. GUARDRAIL §0.1 / §5.1 noted in the
adjacent AGENTS.md).

Implementation notes
--------------------
* Volume averages are r-weighted (axisymmetric): V_cell = 2π·r_centroid·A_cell.
* Centroid samples are taken at the cell whose 2D-centroid is closest to the
  geometric centre of the zone (or of the whole mesh, if no MaterialID).
* Line probes sample by VTK's vtkProbeFilter at n_samples equispaced points on
  the requested line; fields not at sample locations are linearly interpolated.
* Field-availability is checked per VTU; missing fields are reported in the
  summary JSON and the corresponding CSV column is filled with NaN.
* All matplotlib output uses the Agg backend (sandbox-safe).
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import subprocess
import sys
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import vtk  # noqa: E402
import yaml  # noqa: E402
from matplotlib.collections import PolyCollection  # noqa: E402
from vtk.util.numpy_support import vtk_to_numpy  # noqa: E402

# ---------------------------------------------------------------------------
# Constants — geometry and unit conversions are derived in-place; no material
# parameters or expected values are hardcoded in this file (cf. CLAUDE.md §1.1).
# ---------------------------------------------------------------------------
SECONDS_PER_DAY = 86_400.0


def unit_suffix(unit: str) -> str:
    """Sanitize a unit string into a column-name-friendly suffix.

    'MPa' -> 'MPa'; 'kg/m^3' -> 'kg_m3'; '-' -> ''; 'mm' -> 'mm'.
    """
    if unit in ("", "-", None):
        return ""
    return (
        unit.replace("/", "_")
        .replace("^", "")
        .replace(" ", "")
        .replace(".", "")
        .replace("·", "_")
    )


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------
@dataclass
class PvdEntry:
    timestep_index: int
    time_s: float
    vtu_path: Path


@dataclass
class CaseInput:
    pvd_path: Path
    label: str          # e.g. "dd1600" for Model I; "" for III/IV/VII
    rho_dry: Optional[float] = None  # kg/m^3, if known (Model I)


@dataclass
class HistoryRow:
    case_label: str
    timestep_index: int
    time_s: float
    time_days: float
    values: Dict[str, float] = field(default_factory=dict)


# ---------------------------------------------------------------------------
# Provenance helpers
# ---------------------------------------------------------------------------
def md5_of_file(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def git_describe(repo_root: Path) -> str:
    try:
        out = subprocess.run(
            ["git", "-C", str(repo_root), "rev-parse", "HEAD"],
            check=True, capture_output=True, text=True,
        )
        return out.stdout.strip()
    except Exception:
        return "unknown"


# ---------------------------------------------------------------------------
# PVD parsing
# ---------------------------------------------------------------------------
def parse_pvd(pvd_path: Path) -> List[PvdEntry]:
    """Parse an OGS PVD file → list of (timestep_index, time_s, vtu_path)."""
    tree = ET.parse(pvd_path)
    root = tree.getroot()
    entries: List[PvdEntry] = []
    base = pvd_path.parent
    for i, ds in enumerate(root.iter("DataSet")):
        t = float(ds.attrib["timestep"])
        rel = ds.attrib["file"]
        entries.append(PvdEntry(i, t, (base / rel).resolve()))
    entries.sort(key=lambda e: e.time_s)
    return entries


def read_vtu(path: Path) -> vtk.vtkUnstructuredGrid:
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(path))
    reader.Update()
    return reader.GetOutput()


# ---------------------------------------------------------------------------
# Mesh primitives
# ---------------------------------------------------------------------------
def material_ids(grid: vtk.vtkUnstructuredGrid) -> Optional[np.ndarray]:
    arr = grid.GetCellData().GetArray("MaterialIDs")
    return None if arr is None else vtk_to_numpy(arr)


def load_material_ids_from_mesh(
    pvd_path: Path, mesh_basename: Optional[str]
) -> Optional[np.ndarray]:
    """
    OGS snapshot VTUs do not carry MaterialIDs as CellData; that field lives in
    the input mesh referenced by the PRJ. Load it once per case so per-zone
    averaging respects the actual material assignment. Returns None if the
    mesh isn't found (in which case the script falls back to a single zone).
    """
    if mesh_basename is None:
        return None
    mesh_path = pvd_path.parent / mesh_basename
    if not mesh_path.exists():
        return None
    mesh_grid = read_vtu(mesh_path)
    return material_ids(mesh_grid)


def cell_centroids_rz(grid: vtk.vtkUnstructuredGrid) -> np.ndarray:
    """Return (n_cells, 2) array of (r, z) centroids."""
    pts = vtk_to_numpy(grid.GetPoints().GetData())[:, :2]
    n = grid.GetNumberOfCells()
    cents = np.zeros((n, 2))
    for c in range(n):
        cell = grid.GetCell(c)
        pids = [cell.GetPointId(i) for i in range(cell.GetNumberOfPoints())]
        cents[c] = pts[pids].mean(axis=0)
    return cents


def cell_axisym_volumes(grid: vtk.vtkUnstructuredGrid) -> np.ndarray:
    """Cell volumes for an axisymmetric (r, z) mesh: V = 2π·r_centroid·A_rz."""
    pts = vtk_to_numpy(grid.GetPoints().GetData())[:, :2]
    n = grid.GetNumberOfCells()
    vols = np.zeros(n)
    for c in range(n):
        cell = grid.GetCell(c)
        pids = [cell.GetPointId(i) for i in range(cell.GetNumberOfPoints())]
        ring = pts[pids]
        # Shoelace area of polygon in (r,z):
        r = ring[:, 0]
        z = ring[:, 1]
        A = 0.5 * abs(np.sum(r * np.roll(z, -1) - np.roll(r, -1) * z))
        r_c = r.mean()
        vols[c] = 2.0 * np.pi * r_c * A
    return vols


def get_point_array(grid: vtk.vtkUnstructuredGrid, name: str) -> Optional[np.ndarray]:
    arr = grid.GetPointData().GetArray(name)
    return None if arr is None else vtk_to_numpy(arr)


def cell_to_nodal_mask(
    grid: vtk.vtkUnstructuredGrid, cell_mask: np.ndarray
) -> np.ndarray:
    """For each node, True if it belongs to at least one selected cell."""
    n_pts = grid.GetNumberOfPoints()
    in_zone = np.zeros(n_pts, dtype=bool)
    for c in np.flatnonzero(cell_mask):
        cell = grid.GetCell(int(c))
        for i in range(cell.GetNumberOfPoints()):
            in_zone[cell.GetPointId(i)] = True
    return in_zone


# ---------------------------------------------------------------------------
# Reductions
# ---------------------------------------------------------------------------
def mean_compressive(stress_field: np.ndarray) -> np.ndarray:
    """
    Return -tr(σ)/3 per node. OGS Kelvin layout for 2D axisym has 4 components
    [σ_rr, σ_zz, σ_θθ, σ_rz·√2]; mean compressive = -(c0+c1+c2)/3.
    """
    if stress_field.ndim == 1 or stress_field.shape[1] < 3:
        raise ValueError(
            f"mean_compressive expects a tensor field with ≥3 diagonal "
            f"components; got shape {stress_field.shape}."
        )
    return -(stress_field[:, 0] + stress_field[:, 1] + stress_field[:, 2]) / 3.0


def vector_magnitude(vec_field: np.ndarray) -> np.ndarray:
    if vec_field.ndim != 2:
        raise ValueError(f"vector_magnitude expects (n, k); got {vec_field.shape}")
    return np.linalg.norm(vec_field, axis=1)


def void_ratio(porosity: np.ndarray) -> np.ndarray:
    """e = φ/(1-φ). Guard against φ → 1."""
    eps = 1e-12
    return porosity / np.maximum(1.0 - porosity, eps)


# ---------------------------------------------------------------------------
# Volume averaging (axisymmetric, r-weighted)
# ---------------------------------------------------------------------------
def volume_avg(
    grid: vtk.vtkUnstructuredGrid,
    nodal_scalar: np.ndarray,
    cell_mask: Optional[np.ndarray] = None,
) -> float:
    """Volume-weighted (axisymmetric) average of a nodal scalar field."""
    vols = cell_axisym_volumes(grid)
    n_cells = grid.GetNumberOfCells()
    total_v = 0.0
    total_fv = 0.0
    for c in range(n_cells):
        if cell_mask is not None and not cell_mask[c]:
            continue
        cell = grid.GetCell(c)
        pids = [cell.GetPointId(i) for i in range(cell.GetNumberOfPoints())]
        f_c = float(nodal_scalar[pids].mean())
        total_v += vols[c]
        total_fv += vols[c] * f_c
    return float("nan") if total_v == 0.0 else total_fv / total_v


def centroid_value(
    grid: vtk.vtkUnstructuredGrid,
    nodal_scalar: np.ndarray,
    cell_mask: Optional[np.ndarray] = None,
) -> float:
    """Nodal scalar averaged on the single cell closest to the zone centroid."""
    cents = cell_centroids_rz(grid)
    if cell_mask is None:
        mask = np.ones(grid.GetNumberOfCells(), dtype=bool)
    else:
        mask = cell_mask.astype(bool)
    idx = np.flatnonzero(mask)
    if idx.size == 0:
        return float("nan")
    centre = cents[idx].mean(axis=0)
    d2 = ((cents[idx] - centre) ** 2).sum(axis=1)
    best = idx[int(np.argmin(d2))]
    cell = grid.GetCell(int(best))
    pids = [cell.GetPointId(i) for i in range(cell.GetNumberOfPoints())]
    return float(nodal_scalar[pids].mean())


# ---------------------------------------------------------------------------
# Field extraction with reduction
# ---------------------------------------------------------------------------
def extract_nodal_scalar(
    grid: vtk.vtkUnstructuredGrid, field_name: str, reduction: str
) -> Optional[np.ndarray]:
    """
    Read `field_name` from PointData and reduce to a nodal scalar according to
    `reduction`. Returns None if the field is not present.
    """
    raw = get_point_array(grid, field_name)
    if raw is None:
        return None
    if reduction in ("identity", "mean"):
        # `identity` is used for scalar fields (pressure, saturation, …); `mean`
        # is the reduction *over the zone*, applied later. Both consume a scalar.
        if raw.ndim > 1 and raw.shape[1] == 1:
            raw = raw[:, 0]
        if raw.ndim != 1:
            raise ValueError(
                f"{field_name}: expected scalar field for reduction '{reduction}'; "
                f"got shape {raw.shape}"
            )
        return raw.astype(float)
    if reduction == "mean_compressive":
        return mean_compressive(raw).astype(float)
    if reduction == "vector_magnitude":
        return vector_magnitude(raw).astype(float)
    if reduction == "void_ratio_mean":
        # Treat the `porosity` field nodally → e = φ/(1-φ), volume-averaged later.
        if raw.ndim > 1 and raw.shape[1] == 1:
            raw = raw[:, 0]
        return void_ratio(raw).astype(float)
    raise ValueError(f"Unknown reduction '{reduction}' for field '{field_name}'")


# ---------------------------------------------------------------------------
# Line probe (for Model IV interface stress)
# ---------------------------------------------------------------------------
def line_probe_mean_compressive(
    grid: vtk.vtkUnstructuredGrid,
    field_name: str,
    z_m: float,
    r_min_m: float,
    r_max_m: float,
    n_samples: int,
) -> float:
    """
    Sample a tensor field along the line z=const, r ∈ [r_min, r_max] at n_samples
    equispaced points using vtkProbeFilter; reduce each sample to -tr/3; return
    line-mean. Returns NaN if the field is absent or all samples lie outside the
    mesh.
    """
    if get_point_array(grid, field_name) is None:
        return float("nan")
    line = vtk.vtkLineSource()
    line.SetPoint1(r_min_m, z_m, 0.0)
    line.SetPoint2(r_max_m, z_m, 0.0)
    line.SetResolution(n_samples - 1)
    probe = vtk.vtkProbeFilter()
    probe.SetInputConnection(line.GetOutputPort())
    probe.SetSourceData(grid)
    probe.Update()
    sampled = probe.GetOutput()
    arr = sampled.GetPointData().GetArray(field_name)
    if arr is None:
        return float("nan")
    raw = vtk_to_numpy(arr)
    valid_arr = sampled.GetPointData().GetArray("vtkValidPointMask")
    valid = vtk_to_numpy(valid_arr).astype(bool) if valid_arr is not None else np.ones(len(raw), bool)
    if not valid.any():
        return float("nan")
    p = mean_compressive(raw[valid])
    return float(p.mean())


# ---------------------------------------------------------------------------
# Gap-aperture probe (Model III): aperture = gap_nominal − max u_r at interface
# ---------------------------------------------------------------------------
def gap_aperture(
    grid: vtk.vtkUnstructuredGrid,
    gap_nominal_m: float,
    interface_r_m: float,
    z_range_m: Tuple[float, float],
    tol_m: float = 1.0e-9,
) -> float:
    """
    Return the remaining gap (m) = gap_nominal − max(u_r) over the nodes whose
    r-coordinate equals `interface_r_m` within `tol_m`. If no such nodes exist,
    return NaN with a warning.
    """
    pts = vtk_to_numpy(grid.GetPoints().GetData())[:, :2]
    disp = get_point_array(grid, "displacement")
    if disp is None or disp.shape[1] < 2:
        return float("nan")
    mask = (
        (np.abs(pts[:, 0] - interface_r_m) < tol_m)
        & (pts[:, 1] >= z_range_m[0] - tol_m)
        & (pts[:, 1] <= z_range_m[1] + tol_m)
    )
    if not mask.any():
        return float("nan")
    u_r = disp[mask, 0]
    return float(gap_nominal_m - u_r.max())


# ---------------------------------------------------------------------------
# Per-VTU reduction → history row
# ---------------------------------------------------------------------------
def reduce_vtu(
    grid: vtk.vtkUnstructuredGrid,
    probes_spec: Dict[str, Any],
    case_label: str,
    entry: PvdEntry,
    missing_fields: set,
    material_ids_override: Optional[np.ndarray] = None,
) -> HistoryRow:
    row = HistoryRow(
        case_label=case_label,
        timestep_index=entry.timestep_index,
        time_s=entry.time_s,
        time_days=entry.time_s / SECONDS_PER_DAY,
    )

    # Prefer MaterialIDs from the input mesh (OGS snapshots omit CellData);
    # fall back to whatever the snapshot carries if no override was supplied.
    mat = material_ids_override if material_ids_override is not None else material_ids(grid)
    zones = probes_spec.get("zones", [])

    # Build masks per zone:
    zone_cell_mask: Dict[str, np.ndarray] = {}
    n_cells = grid.GetNumberOfCells()
    if mat is None:
        # Single-material mesh (e.g. Model I): one zone covers everything.
        for z in zones:
            zone_cell_mask[z["id"]] = np.ones(n_cells, dtype=bool)
    else:
        for z in zones:
            zone_cell_mask[z["id"]] = (mat == int(z["material_id"]))

    # Time-series quantities:
    for ts in probes_spec.get("time_series", []):
        key = ts["key"]
        # Probe-style entries route to dedicated samplers below.
        if "probe" in ts:
            continue
        field_name = ts["field"]
        reduction = ts["reduction"]
        scale = float(ts.get("scale", 1.0))
        per_zone = bool(ts.get("per_zone", False))
        suffix = unit_suffix(ts.get("unit", ""))
        suf = f"_{suffix}" if suffix else ""

        nodal = extract_nodal_scalar(grid, field_name, reduction)
        if nodal is None:
            missing_fields.add(field_name)
            if per_zone:
                for z in zones:
                    row.values[f"{key}_{z['id']}_volavg{suf}"] = float("nan")
                    row.values[f"{key}_{z['id']}_centroid{suf}"] = float("nan")
            else:
                row.values[f"{key}{suf}"] = float("nan")
            continue

        if per_zone:
            for z in zones:
                mask = zone_cell_mask[z["id"]]
                row.values[f"{key}_{z['id']}_volavg{suf}"] = volume_avg(grid, nodal, mask) * scale
                row.values[f"{key}_{z['id']}_centroid{suf}"] = centroid_value(grid, nodal, mask) * scale
        else:
            row.values[f"{key}{suf}"] = volume_avg(grid, nodal, None) * scale

    # Probe-style time-series (line probes, gap apertures, etc.):
    probes_block = probes_spec.get("probes", {})
    for ts in probes_spec.get("time_series", []):
        if "probe" not in ts:
            continue
        key = ts["key"]
        probe_name = ts["probe"]
        spec = probes_block.get(probe_name)
        if spec is None:
            row.values[f"{key}"] = float("nan")
            continue
        ptype = spec["type"]
        if ptype == "line_probe":
            g = spec["geometry"]
            row.values[f"{key}_line_mean_MPa"] = (
                line_probe_mean_compressive(
                    grid,
                    spec["field"],
                    float(g["z_m"]),
                    float(g["r_min_m"]),
                    float(g["r_max_m"]),
                    int(g["n_samples"]),
                )
                * 1.0e-6
            )
        elif ptype == "radial_extension_at_interface":
            row.values[f"{key}_mm"] = (
                gap_aperture(
                    grid,
                    float(spec["gap_nominal_m"]),
                    float(spec["interface_r_m"]),
                    tuple(spec["z_range_m"]),
                )
                * 1.0e3
            )
        else:
            row.values[f"{key}"] = float("nan")

    return row


# ---------------------------------------------------------------------------
# Triptych rendering
# ---------------------------------------------------------------------------
def render_field_panel(
    ax: plt.Axes,
    grid: vtk.vtkUnstructuredGrid,
    panel_spec: Dict[str, Any],
    title: str,
) -> None:
    """Render one field panel of the final-state triptych."""
    nodal = extract_nodal_scalar(grid, panel_spec["source"], panel_spec["reduction"])
    pts = vtk_to_numpy(grid.GetPoints().GetData())[:, :2]
    scale = float(panel_spec.get("scale", 1.0))
    unit = panel_spec.get("unit", "")
    if nodal is None:
        ax.set_title(f"{title}\n(field '{panel_spec['source']}' missing)")
        ax.set_xlabel("r [m]")
        ax.set_ylabel("z [m]")
        return
    values = nodal * scale

    polys: List[np.ndarray] = []
    cell_values: List[float] = []
    for c in range(grid.GetNumberOfCells()):
        cell = grid.GetCell(c)
        pids = [cell.GetPointId(i) for i in range(cell.GetNumberOfPoints())]
        polys.append(pts[pids])
        cell_values.append(float(values[pids].mean()))

    coll = PolyCollection(polys, array=np.array(cell_values), edgecolors="0.6", linewidths=0.15)
    ax.add_collection(coll)
    ax.set_xlim(pts[:, 0].min(), pts[:, 0].max())
    ax.set_ylim(pts[:, 1].min(), pts[:, 1].max())
    ax.set_aspect("equal")
    ax.set_xlabel("r [m]")
    ax.set_ylabel("z [m]")
    ax.set_title(title)
    cbar = plt.colorbar(coll, ax=ax, fraction=0.046, pad=0.04)
    if unit:
        cbar.set_label(unit)


def render_triptych(
    grid: vtk.vtkUnstructuredGrid,
    probes_spec: Dict[str, Any],
    out_path: Path,
    suptitle: str,
) -> None:
    panels = probes_spec["triptych"]["fields"]
    n = len(panels)
    fig, axes = plt.subplots(1, n, figsize=(4.5 * n, 5.0))
    if n == 1:
        axes = [axes]
    for ax, p in zip(axes, panels):
        title = p["panel"]
        render_field_panel(ax, grid, p, title)
    fig.suptitle(suptitle)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.96))
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Time-series figures
# ---------------------------------------------------------------------------
def add_marker_lines(ax: plt.Axes, marker_days: Sequence[float]) -> None:
    for d in marker_days:
        ax.axvline(d, color="0.7", linestyle="--", linewidth=0.7, zorder=0)


def _pretty_label(column_name: str) -> str:
    """Strip a unit suffix off a column name for legend use; humanize underscores."""
    # Try to peel a known unit suffix:
    for suf in ("_MPa", "_kg_m3", "_mm", "_m"):
        if column_name.endswith(suf):
            return column_name[: -len(suf)].replace("_", " ")
    return column_name.replace("_", " ")


def _axis_label(column_names: Sequence[str]) -> str:
    """Build a clean axis label from one or more column names.

    Uses the shared unit suffix if all columns end with the same one; otherwise
    falls back to a comma-joined humanized list.
    """
    units = []
    for c in column_names:
        for suf in ("_MPa", "_kg_m3", "_mm", "_m"):
            if c.endswith(suf):
                units.append(suf.lstrip("_"))
                break
        else:
            units.append("")
    if len(set(units)) == 1 and units[0]:
        labels = [_pretty_label(c) for c in column_names]
        if len({lbl.split()[0] for lbl in labels}) == 1:
            common = labels[0].split()[0]
            return f"{common} [{units[0]}]"
        return f"{', '.join(labels)} [{units[0]}]"
    return ", ".join(_pretty_label(c) for c in column_names)


def render_time_series_figure(
    history: List[HistoryRow],
    xaxis: str,
    yaxis,
    title: str,
    out_path: Path,
    marker_days: Sequence[float],
    cases: Sequence[str],
) -> None:
    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    ys = [yaxis] if isinstance(yaxis, str) else list(yaxis)
    multi_case = len(cases) > 1
    multi_y = len(ys) > 1
    for case in cases:
        rows_case = [r for r in history if r.case_label == case]
        if not rows_case:
            continue
        x = np.array([_lookup(r, xaxis) for r in rows_case])
        for y_key in ys:
            y = np.array([_lookup(r, y_key) for r in rows_case])
            if multi_case and multi_y:
                label = f"{case} / {_pretty_label(y_key)}"
            elif multi_case:
                label = case
            elif multi_y:
                label = _pretty_label(y_key)
            else:
                label = _pretty_label(y_key)
            ax.plot(x, y, label=label)
            if xaxis == "time_days" and marker_days:
                xm = np.array(marker_days, dtype=float)
                ym = np.interp(xm, x, y, left=np.nan, right=np.nan)
                ax.plot(xm, ym, "o", color=ax.lines[-1].get_color(), zorder=5)
    if xaxis == "time_days":
        add_marker_lines(ax, marker_days)
    ax.set_xlabel("time [days]" if xaxis == "time_days" else _axis_label([xaxis]))
    ax.set_ylabel(_axis_label(ys))
    ax.set_title(title)
    ax.legend(fontsize=9, loc="best")
    ax.grid(True, linewidth=0.3, alpha=0.4)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def _lookup(row: HistoryRow, key: str) -> float:
    if key == "time_days":
        return row.time_days
    if key == "time_s":
        return row.time_s
    return row.values.get(key, float("nan"))


# ---------------------------------------------------------------------------
# Model-I Villar benchmark figure (dd sweep)
# ---------------------------------------------------------------------------
def villar_target_MPa(rho_d_kg_m3: float) -> float:
    """Villar / Lloret correlation Ps[MPa] = exp(6.77·ρ_d[g/cm^3] − 9.07).

    Source: documented in the MS33 skill anchors-ms33-workflow and in the
    audit_extract_ms33.py docstring; this is the EURAD-2 MS33 reference
    correlation, NOT a fitted result of these simulations.
    """
    return float(np.exp(6.77 * rho_d_kg_m3 / 1000.0 - 9.07))


def render_villar_benchmark(
    history: List[HistoryRow],
    cases: Sequence[CaseInput],
    out_path: Path,
) -> None:
    rho_grid = np.linspace(1300.0, 1900.0, 200) / 1000.0
    ps_curve = np.array([villar_target_MPa(r * 1000.0) for r in rho_grid])

    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    ax.plot(rho_grid * 1000.0, ps_curve, "k--", label="Villar correlation")
    # OGS simulation final-state points: pick the last row per case.
    for c in cases:
        if c.rho_dry is None:
            continue
        rows_case = [r for r in history if r.case_label == c.label]
        if not rows_case:
            continue
        ps_sim = rows_case[-1].values.get("mean_stress_bentonite_volavg_MPa", float("nan"))
        ax.plot(c.rho_dry, ps_sim, "o", markersize=8, label=f"OGS {c.label}: {ps_sim:.3f} MPa")
    ax.set_xlabel("dry density ρ_d [kg/m^3]")
    ax.set_ylabel("swelling pressure Ps [MPa]")
    ax.set_yscale("log")
    ax.set_title("Model I — Villar benchmark")
    ax.legend(fontsize=8, loc="best")
    ax.grid(True, which="both", linewidth=0.3, alpha=0.4)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


# ---------------------------------------------------------------------------
# CSV / JSON writers
# ---------------------------------------------------------------------------
def write_history_csv(history: List[HistoryRow], out_path: Path) -> List[str]:
    if not history:
        out_path.write_text("")
        return []
    columns: List[str] = ["case_label", "timestep_index", "time_s", "time_days"]
    seen = set(columns)
    for r in history:
        for k in r.values.keys():
            if k not in seen:
                columns.append(k)
                seen.add(k)
    with out_path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(columns)
        for r in history:
            row = [r.case_label, r.timestep_index, r.time_s, r.time_days]
            row += [r.values.get(k, "") for k in columns[4:]]
            w.writerow(row)
    return columns


def write_summary_json(
    out_path: Path,
    *,
    model: str,
    tag: str,
    cases: Sequence[CaseInput],
    probes_path: Path,
    probes_spec: Dict[str, Any],
    csv_columns: Sequence[str],
    missing_fields: set,
    figures: Sequence[Path],
    ogs_repo: Path,
) -> None:
    payload = {
        "tool": "audit_extract_ms33.py",
        "model": model,
        "tag": tag,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "ogs_git_head": git_describe(ogs_repo),
        "probes_yaml": {
            "path": str(probes_path),
            "md5": md5_of_file(probes_path),
        },
        "cases": [
            {
                "label": c.label,
                "pvd_path": str(c.pvd_path),
                "pvd_md5": md5_of_file(c.pvd_path),
                "rho_dry_kg_m3": c.rho_dry,
            }
            for c in cases
        ],
        "csv_columns": list(csv_columns),
        "missing_fields": sorted(missing_fields),
        "figures": [str(p) for p in figures],
        "probes_spec_recap": {
            "zones": probes_spec.get("zones", []),
            "markers_days": probes_spec.get("markers_days", []),
            "stress_column_label": "mean_stress_MPa",
            "stress_column_doc": (
                "Under BishopsSaturationCutoff cutoff_value=1, sigma_eff is "
                "identically sigma_total throughout these simulations; the "
                "VTU 'sigma' field maps to both interpretations."
            ),
        },
    }
    out_path.write_text(json.dumps(payload, indent=2))


# ---------------------------------------------------------------------------
# Main driver
# ---------------------------------------------------------------------------
def main(argv: Optional[Sequence[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument("--model", required=True, choices=["I", "III", "IV", "VII"])
    p.add_argument(
        "--pvd",
        required=True,
        action="append",
        type=Path,
        help="path to a .pvd file; may be repeated (Model I dd-sweep)",
    )
    p.add_argument(
        "--rho-dry",
        action="append",
        type=float,
        help="dry density [kg/m^3] for the corresponding --pvd (Model I); "
             "repeat once per --pvd",
    )
    p.add_argument("--tag", required=True, choices=["mfront", "native"])
    p.add_argument("--probes", required=True, type=Path)
    p.add_argument("-o", "--out", required=True, type=Path)
    args = p.parse_args(argv)

    out_dir: Path = args.out
    out_dir.mkdir(parents=True, exist_ok=True)
    probes_spec = yaml.safe_load(args.probes.read_text())

    # Build CaseInput list:
    pvds: List[Path] = args.pvd
    rhos = args.rho_dry or [None] * len(pvds)
    if len(rhos) != len(pvds):
        print(
            f"--rho-dry was given {len(rhos)} times but --pvd was given {len(pvds)} "
            f"times; they must match (one --rho-dry per --pvd) or --rho-dry omitted.",
            file=sys.stderr,
        )
        return 2
    cases: List[CaseInput] = []
    for path, rho in zip(pvds, rhos):
        label = ""
        if rho is not None:
            label = f"dd{int(round(rho))}"
        else:
            label = path.stem
        cases.append(CaseInput(pvd_path=path.resolve(), label=label, rho_dry=rho))

    # Reduce all timesteps of all cases:
    history: List[HistoryRow] = []
    missing_fields: set = set()
    final_grids: Dict[str, vtk.vtkUnstructuredGrid] = {}
    mesh_basename = probes_spec.get("mesh")  # YAML-declared input mesh
    expects_multi_zone = len(probes_spec.get("zones", [])) > 1
    for case in cases:
        entries = parse_pvd(case.pvd_path)
        mat_ids = load_material_ids_from_mesh(case.pvd_path, mesh_basename)
        # Only warn if the probe spec defines multiple zones — single-material
        # meshes legitimately omit MaterialIDs.
        if mat_ids is None and mesh_basename is not None and expects_multi_zone:
            print(
                f"  (warn) MaterialIDs not loadable from {mesh_basename}; "
                f"per-zone averaging will fall back to single-zone.",
                file=sys.stderr,
            )
        for e in entries:
            grid = read_vtu(e.vtu_path)
            row = reduce_vtu(grid, probes_spec, case.label, e, missing_fields, mat_ids)
            history.append(row)
        if entries:
            final_grids[case.label] = read_vtu(entries[-1].vtu_path)

    # Stamp out file names with model + tag:
    prefix = f"model{args.model}_{args.tag}"
    history_csv = out_dir / f"{prefix}_history.csv"
    summary_json = out_dir / f"{prefix}_summary.json"
    figures_out: List[Path] = []

    columns = write_history_csv(history, history_csv)

    # Triptych(s): one per case (Model I has up to three; III/IV/VII have one).
    # If there's a single case and no explicit rho-dry, drop the redundant
    # case suffix so the filename stays clean (e.g. modelIV_mfront_triptych.png
    # rather than modelIV_mfront_ms33_modelIV_pellets_triptych.png).
    single_case_no_rho = (len(cases) == 1) and (cases[0].rho_dry is None)
    for case in cases:
        grid = final_grids.get(case.label)
        if grid is None:
            continue
        if single_case_no_rho:
            triptych_path = out_dir / f"{prefix}_triptych.png"
            suptitle = f"Model {args.model} ({args.tag}) — final state"
        else:
            triptych_path = out_dir / f"{prefix}_{case.label}_triptych.png"
            suptitle = f"Model {args.model} ({args.tag}, {case.label}) — final state"
        render_triptych(grid, probes_spec, triptych_path, suptitle=suptitle)
        figures_out.append(triptych_path)

    # Time-series / p–s / villar figures from YAML 'figures' block:
    marker_days = probes_spec.get("markers_days", [])
    case_labels = [c.label for c in cases]
    for fig_spec in probes_spec.get("figures", []):
        key = fig_spec["key"]
        if fig_spec.get("kind") == "dd_sweep" or key == "villar_benchmark":
            if args.model == "I":
                fp = out_dir / f"{prefix}_villar_benchmark.pdf"
                render_villar_benchmark(history, cases, fp)
                figures_out.append(fp)
            continue
        fp = out_dir / f"{prefix}_{key}.pdf"
        render_time_series_figure(
            history,
            fig_spec["xaxis"],
            fig_spec["yaxis"],
            fig_spec.get("title", key),
            fp,
            marker_days,
            case_labels,
        )
        figures_out.append(fp)

    ogs_repo = Path(__file__).resolve().parents[3]
    write_summary_json(
        summary_json,
        model=args.model,
        tag=args.tag,
        cases=cases,
        probes_path=args.probes.resolve(),
        probes_spec=probes_spec,
        csv_columns=columns,
        missing_fields=missing_fields,
        figures=figures_out,
        ogs_repo=ogs_repo,
    )

    print(f"wrote {history_csv}")
    print(f"wrote {summary_json}")
    for f in figures_out:
        print(f"wrote {f}")
    if missing_fields:
        print(f"  (note) fields missing from VTUs: {sorted(missing_fields)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
