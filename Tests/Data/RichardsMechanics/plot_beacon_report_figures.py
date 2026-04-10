#!/usr/bin/env python3
"""Generate BEACON-style comparison figures for native and MFront runs.

This helper script re-runs selected BEACON cases (1a01, 1b, 1c), extracts
stress and density profile indicators from VTK outputs, and writes report-like
plots plus CSV/JSON summary artifacts under ``Tests/Data/RichardsMechanics``.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import tempfile
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", tempfile.mkdtemp(prefix="mplconfig-"))
os.environ.setdefault("XDG_CACHE_HOME", tempfile.mkdtemp(prefix="xdg-cache-"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy


HERE = Path(__file__).resolve().parent
NATIVE_DEFAULT = Path("/Users/vinaykumar/git/build/release-native-beacon/bin/ogs")
MFRONT_DEFAULT = Path("/Users/vinaykumar/git/build/release-mfront-tpm/bin/ogs")


CASE_CONFIGS = {
    "1a01": {
        "native_project": HERE / "beacon_1a01_inflow_unstructured_batch.prj",
        "bridge_project": HERE / "beacon_1a01_dsm_micromacro_mcc_inflow_unstructured_batch.prj",
        "native_prefix": "beacon_1a01_native_reportplot",
        "bridge_prefix": "beacon_1a01_bridge_reportplot",
        "t_end_s": 1.0e5,
        "output_times_s": [
            10.0,
            30.0,
            100.0,
            300.0,
            1.0e3,
            3.0e3,
            1.0e4,
            3.0e4,
            5.0e4,
            7.5e4,
            1.0e5,
        ],
        "profile_bins": 4,
        "profile_height_mm": [2.5, 7.5, 12.5, 17.5],
        "report_profile_kg_m3": [1466.0, 1454.0, 1427.0, 1353.0],
        "report_axial_kpa": 604.0,
        "report_radial_kpa": 994.0,
    },
    "1b": {
        "native_project": HERE / "beacon_1b_unstructured_batch.prj",
        "bridge_project": HERE / "beacon_1b_dsm_micromacro_mcc_unstructured_batch.prj",
        "native_prefix": "beacon_1b_native_reportplot",
        "bridge_prefix": "beacon_1b_bridge_reportplot",
        "t_end_s": 4.32e7,
        "output_times_s": [
            1.08e5,
            5.40e5,
            1.08e6,
            2.16e6,
            5.40e6,
            1.08e7,
            1.62e7,
            2.16e7,
            3.24e7,
            4.32e7,
        ],
        "profile_bins": 10,
        "report_density_mean_kg_m3": 1520.0,
        "report_density_note": "nonzero and stabilized after about 500 days",
    },
    "1c": {
        "bridge_project": HERE / "beacon_1c_dsm_micromacro_mcc_bridge.prj",
        "bridge_prefix": "beacon_1c_bridge_reportplot",
        "t_end_s": 1.0e3,
        "output_times_s": [50.0, 100.0, 150.0, 250.0, 400.0, 550.0, 700.0, 850.0, 1.0e3],
        "profile_bins": 20,
        "block_height_m": 0.0485,
    },
}


@dataclass
class BoundarySeries:
    label: str
    times_s: np.ndarray
    axial_kpa: np.ndarray
    radial_kpa: np.ndarray


@dataclass
class ZoneStressSeries:
    label: str
    times_s: np.ndarray
    block_axial_kpa: np.ndarray
    block_radial_kpa: np.ndarray
    pellet_axial_kpa: np.ndarray
    pellet_radial_kpa: np.ndarray


def write_project_copy(source: Path, target: Path, prefix: str, output_times_s: list[float]) -> None:
    """Write a temporary project with absolute mesh/geometry and custom output prefix/times."""
    target.parent.mkdir(parents=True, exist_ok=True)
    tree = ET.parse(source)
    root = tree.getroot()

    mesh = root.find("./mesh")
    geometry = root.find("./geometry")
    if mesh is None or geometry is None or mesh.text is None or geometry.text is None:
        raise RuntimeError(f"Could not resolve mesh/geometry in {source}")

    mesh.text = str((source.parent / mesh.text).resolve())
    geometry.text = str((source.parent / geometry.text).resolve())

    output = root.find("./time_loop/output")
    if output is None:
        raise RuntimeError(f"No output block found in {source}")

    prefix_node = output.find("prefix")
    if prefix_node is None:
        prefix_node = ET.SubElement(output, "prefix")
    prefix_node.text = prefix

    fixed_times = output.find("fixed_output_times")
    if fixed_times is None:
        fixed_times = ET.SubElement(output, "fixed_output_times")
    fixed_times.text = " ".join(f"{t:.6f}" for t in output_times_s)

    tree.write(target, encoding="ISO-8859-1", xml_declaration=True)


def run_ogs(ogs: Path, project: Path, output_dir: Path) -> None:
    """Execute one OGS project and raise with captured logs on failure."""
    output_dir.mkdir(parents=True, exist_ok=True)
    result = subprocess.run(
        [str(ogs), "-o", str(output_dir), str(project)],
        cwd=HERE,
        text=True,
        capture_output=True,
        check=False,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"ogs failed for '{project.name}' with '{ogs}'.\n"
            f"stdout:\n{result.stdout}\n"
            f"stderr:\n{result.stderr}"
        )


def time_from_filename(path: Path) -> float:
    marker = "_t_"
    if marker not in path.stem:
        raise RuntimeError(f"Could not parse time from {path.name}")
    return float(path.stem.split(marker, 1)[1])


def load_vtu(path: Path) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(path))
    reader.Update()
    grid = reader.GetOutput()
    points = np.array([grid.GetPoint(i) for i in range(grid.GetNumberOfPoints())])
    point_data = {
        grid.GetPointData().GetArrayName(i): vtk_to_numpy(grid.GetPointData().GetArray(i))
        for i in range(grid.GetPointData().GetNumberOfArrays())
    }
    return points, point_data


def boundary_mean(points: np.ndarray, values: np.ndarray, axis: int, side: str) -> np.ndarray:
    target = points[:, axis].max() if side == "max" else points[:, axis].min()
    mask = np.isclose(points[:, axis], target, atol=1e-10)
    return np.asarray(values[mask]).mean(axis=0)


def binned_profile(
    points: np.ndarray,
    values: np.ndarray,
    *,
    n_bins: int,
    y_max_override: float | None = None,
) -> list[dict[str, float]]:
    y = points[:, 1]
    ymin = float(y.min())
    ymax = float(y.max() if y_max_override is None else y_max_override)
    edges = np.linspace(ymin, ymax, n_bins + 1)
    out: list[dict[str, float]] = []
    for i in range(n_bins):
        lo = edges[i]
        hi = edges[i + 1]
        if i == n_bins - 1:
            mask = (y >= lo) & (y <= hi)
        else:
            mask = (y >= lo) & (y < hi)
        out.append(
            {
                "center_mm": 1000.0 * 0.5 * (lo + hi),
                "mean": float(np.asarray(values[mask]).mean()),
                "count": int(mask.sum()),
            }
        )
    return out


def collect_boundary_series(run_dir: Path, prefix: str, label: str) -> BoundarySeries:
    """Read all VTU files for one prefix and return top/outer boundary stress time series."""
    paths = sorted(run_dir.glob(f"{prefix}_t_*.vtu"), key=time_from_filename)
    if not paths:
        raise RuntimeError(f"No VTU outputs found for {prefix} in {run_dir}")

    times_s = []
    axial_kpa = []
    radial_kpa = []
    for path in paths:
        points, data = load_vtu(path)
        sigma_top = boundary_mean(points, data["sigma"], axis=1, side="max")
        sigma_outer = boundary_mean(points, data["sigma"], axis=0, side="max")
        times_s.append(time_from_filename(path))
        axial_kpa.append(abs(float(sigma_top[1])) / 1000.0)
        radial_kpa.append(abs(float(sigma_outer[0])) / 1000.0)

    order = np.argsort(np.asarray(times_s))
    return BoundarySeries(
        label=label,
        times_s=np.asarray(times_s)[order],
        axial_kpa=np.asarray(axial_kpa)[order],
        radial_kpa=np.asarray(radial_kpa)[order],
    )


def collect_density_profile(run_dir: Path, prefix: str, n_bins: int) -> list[dict[str, float]]:
    paths = sorted(run_dir.glob(f"{prefix}_t_*.vtu"), key=time_from_filename)
    if not paths:
        raise RuntimeError(f"No VTU outputs found for {prefix} in {run_dir}")
    points, data = load_vtu(paths[-1])
    return binned_profile(points, data["dry_density_solid"], n_bins=n_bins)


def collect_zone_stress_series(run_dir: Path, prefix: str, label: str, block_height_m: float) -> ZoneStressSeries:
    paths = sorted(run_dir.glob(f"{prefix}_t_*.vtu"), key=time_from_filename)
    if not paths:
        raise RuntimeError(f"No VTU outputs found for {prefix} in {run_dir}")

    times_s = []
    block_axial_kpa = []
    block_radial_kpa = []
    pellet_axial_kpa = []
    pellet_radial_kpa = []

    for path in paths:
        points, data = load_vtu(path)
        sigma = np.asarray(data["sigma"])
        y = points[:, 1]
        block_mask = y <= (block_height_m + 1e-10)
        pellet_mask = y >= (block_height_m - 1e-10)

        times_s.append(time_from_filename(path))
        block_axial_kpa.append(abs(float(np.asarray(sigma[block_mask, 1]).mean())) / 1000.0)
        block_radial_kpa.append(abs(float(np.asarray(sigma[block_mask, 0]).mean())) / 1000.0)
        pellet_axial_kpa.append(abs(float(np.asarray(sigma[pellet_mask, 1]).mean())) / 1000.0)
        pellet_radial_kpa.append(abs(float(np.asarray(sigma[pellet_mask, 0]).mean())) / 1000.0)

    order = np.argsort(np.asarray(times_s))
    return ZoneStressSeries(
        label=label,
        times_s=np.asarray(times_s)[order],
        block_axial_kpa=np.asarray(block_axial_kpa)[order],
        block_radial_kpa=np.asarray(block_radial_kpa)[order],
        pellet_axial_kpa=np.asarray(pellet_axial_kpa)[order],
        pellet_radial_kpa=np.asarray(pellet_radial_kpa)[order],
    )


def collect_multi_profile(
    run_dir: Path,
    prefix: str,
    n_bins: int,
) -> list[dict[str, float]]:
    paths = sorted(run_dir.glob(f"{prefix}_t_*.vtu"), key=time_from_filename)
    if not paths:
        raise RuntimeError(f"No VTU outputs found for {prefix} in {run_dir}")
    points, data = load_vtu(paths[-1])
    porosity_profile = binned_profile(points, data["porosity"], n_bins=n_bins)
    transport_profile = binned_profile(points, data["transport_porosity"], n_bins=n_bins)
    saturation_profile = binned_profile(points, data["saturation"], n_bins=n_bins)

    rows = []
    for i in range(n_bins):
        rows.append(
            {
                "center_mm": porosity_profile[i]["center_mm"],
                "porosity": porosity_profile[i]["mean"],
                "transport_porosity": transport_profile[i]["mean"],
                "saturation": saturation_profile[i]["mean"],
            }
        )
    return rows


def collect_profile_by_unique_height(run_dir: Path, prefix: str) -> list[dict[str, float]]:
    paths = sorted(run_dir.glob(f"{prefix}_t_*.vtu"), key=time_from_filename)
    if not paths:
        raise RuntimeError(f"No VTU outputs found for {prefix} in {run_dir}")
    points, data = load_vtu(paths[-1])
    y = np.round(points[:, 1], decimals=10)
    rows = []
    for height in np.unique(y):
        mask = np.isclose(y, height, atol=1e-10)
        rows.append(
            {
                "center_mm": float(1000.0 * height),
                "porosity": float(np.asarray(data["porosity"][mask]).mean()),
                "transport_porosity": float(np.asarray(data["transport_porosity"][mask]).mean()),
                "saturation": float(np.asarray(data["saturation"][mask]).mean()),
            }
        )
    return rows


def write_boundary_series_csv(path: Path, series_list: list[BoundarySeries]) -> None:
    with path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["series", "time_s", "time_h", "time_d", "axial_kpa", "radial_kpa"])
        for series in series_list:
            for row in zip(
                series.times_s,
                series.times_s / 3600.0,
                series.times_s / 86400.0,
                series.axial_kpa,
                series.radial_kpa,
                strict=True,
            ):
                writer.writerow([series.label, *row])


def write_density_profile_csv(path: Path, report_rows: list[dict[str, float]] | None, native_rows: list[dict[str, float]] | None, bridge_rows: list[dict[str, float]] | None) -> None:
    with path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["height_mm", "report_density_kg_m3", "native_density_kg_m3", "bridge_density_kg_m3"])
        n_rows = 0
        for rows in (report_rows, native_rows, bridge_rows):
            if rows is not None:
                n_rows = max(n_rows, len(rows))
        for i in range(n_rows):
            height = None
            if report_rows is not None and i < len(report_rows):
                height = report_rows[i]["center_mm"]
            elif native_rows is not None and i < len(native_rows):
                height = native_rows[i]["center_mm"]
            elif bridge_rows is not None and i < len(bridge_rows):
                height = bridge_rows[i]["center_mm"]
            writer.writerow(
                [
                    height,
                    None if report_rows is None else report_rows[i]["mean"],
                    None if native_rows is None else native_rows[i]["mean"],
                    None if bridge_rows is None else bridge_rows[i]["mean"],
                ]
            )


def write_zone_series_csv(path: Path, series: ZoneStressSeries) -> None:
    with path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "series",
                "time_s",
                "time_min",
                "block_axial_kpa",
                "block_radial_kpa",
                "pellet_axial_kpa",
                "pellet_radial_kpa",
            ]
        )
        for row in zip(
            series.times_s,
            series.times_s / 60.0,
            series.block_axial_kpa,
            series.block_radial_kpa,
            series.pellet_axial_kpa,
            series.pellet_radial_kpa,
            strict=True,
        ):
            writer.writerow([series.label, *row])


def write_zone_profile_csv(path: Path, rows: list[dict[str, float]]) -> None:
    with path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["height_mm", "porosity", "transport_porosity", "saturation"])
        for row in rows:
            writer.writerow([row["center_mm"], row["porosity"], row["transport_porosity"], row["saturation"]])


def plot_1a01_panel(
    path: Path,
    series_list: list[BoundarySeries],
    native_profile: list[dict[str, float]],
    bridge_profile: list[dict[str, float]],
    config: dict,
) -> None:
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(11.8, 4.6))

    colors = {"Native OGS": "#d62728", "MFront OGS": "#1f77b4"}
    for series in series_list:
        color = colors.get(series.label, "#333333")
        ax0.plot(series.times_s / 3600.0, series.axial_kpa, color=color, linewidth=2.0, label=f"{series.label} axial")
        ax0.plot(
            series.times_s / 3600.0,
            series.radial_kpa,
            color=color,
            linewidth=1.8,
            linestyle="--",
            label=f"{series.label} radial",
        )

    ax0.axhline(config["report_axial_kpa"], color="#666666", linewidth=1.5, linestyle="-.", label="BEACON axial target")
    ax0.axhline(config["report_radial_kpa"], color="#999999", linewidth=1.5, linestyle=":", label="BEACON radial target")
    ax0.set_title("1a01 Stage-1 Swelling Pressure")
    ax0.set_xlabel("Time (h)")
    ax0.set_ylabel("Swelling pressure (kPa)")
    ax0.set_xlim(left=0.0)
    ax0.set_ylim(bottom=0.0)
    ax0.grid(True, alpha=0.3)
    ax0.legend(frameon=True, fontsize=8)

    report_heights = np.asarray(config["profile_height_mm"])
    report_density = np.asarray(config["report_profile_kg_m3"])
    ax1.plot(report_density, report_heights, color="#444444", linewidth=2.0, marker="s", label="BEACON post-mortem")
    ax1.plot(
        [row["mean"] for row in native_profile],
        [row["center_mm"] for row in native_profile],
        color=colors["Native OGS"],
        linewidth=1.9,
        marker="o",
        label="Native OGS stage-1",
    )
    ax1.plot(
        [row["mean"] for row in bridge_profile],
        [row["center_mm"] for row in bridge_profile],
        color=colors["MFront OGS"],
        linewidth=1.9,
        marker="^",
        label="MFront OGS stage-1",
    )
    ax1.set_title("1a01 Dry-Density Profile")
    ax1.set_xlabel(r"Dry density (kg/m$^3$)")
    ax1.set_ylabel("Height from bottom (mm)")
    ax1.set_ylim(0.0, 20.0)
    ax1.grid(True, alpha=0.3)
    ax1.legend(frameon=True, fontsize=8, loc="lower left")

    fig.tight_layout()
    fig.savefig(path, dpi=220)
    plt.close(fig)


def plot_1b_panel(
    path: Path,
    series_list: list[BoundarySeries],
    native_profile: list[dict[str, float]],
    bridge_profile: list[dict[str, float]],
    config: dict,
) -> None:
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(11.8, 4.6))

    colors = {"Native OGS": "#d62728", "MFront OGS": "#1f77b4"}
    for series in series_list:
        color = colors.get(series.label, "#333333")
        ax0.plot(series.times_s / 86400.0, series.axial_kpa, color=color, linewidth=2.0, label=f"{series.label} axial")
        ax0.plot(
            series.times_s / 86400.0,
            series.radial_kpa,
            color=color,
            linewidth=1.8,
            linestyle="--",
            label=f"{series.label} radial",
        )

    ax0.set_title("1b Swelling Pressure Evolution")
    ax0.set_xlabel("Time (days)")
    ax0.set_ylabel("Swelling pressure (kPa)")
    ax0.set_xlim(left=0.0, right=max(series_list[0].times_s.max() / 86400.0, 500.0))
    ax0.set_ylim(bottom=0.0)
    ax0.grid(True, alpha=0.3)
    ax0.text(
        0.02,
        0.96,
        "BEACON note: nonzero and stabilized after about 500 days",
        transform=ax0.transAxes,
        va="top",
        ha="left",
        fontsize=8,
        bbox={"facecolor": "white", "edgecolor": "#cccccc", "alpha": 0.85},
    )
    ax0.legend(frameon=True, fontsize=8)

    ax1.plot(
        [row["mean"] for row in native_profile],
        [row["center_mm"] for row in native_profile],
        color=colors["Native OGS"],
        linewidth=1.9,
        marker="o",
        label="Native OGS",
    )
    ax1.plot(
        [row["mean"] for row in bridge_profile],
        [row["center_mm"] for row in bridge_profile],
        color=colors["MFront OGS"],
        linewidth=1.9,
        marker="^",
        label="MFront OGS",
    )
    ax1.axvline(
        config["report_density_mean_kg_m3"],
        color="#666666",
        linewidth=1.6,
        linestyle="-.",
        label="BEACON mean dry density",
    )
    ax1.set_title("1b Dry-Density Profile")
    ax1.set_xlabel(r"Dry density (kg/m$^3$)")
    ax1.set_ylabel("Height from bottom (mm)")
    ax1.set_ylim(0.0, 105.15)
    ax1.grid(True, alpha=0.3)
    ax1.legend(frameon=True, fontsize=8, loc="lower right")

    fig.tight_layout()
    fig.savefig(path, dpi=220)
    plt.close(fig)


def plot_1c_panel(
    path: Path,
    zone_series: ZoneStressSeries,
    profile_rows: list[dict[str, float]],
    block_height_m: float,
) -> None:
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(11.8, 4.6))

    time_min = zone_series.times_s / 60.0
    ax0.plot(time_min, zone_series.block_axial_kpa, color="#2ca02c", linewidth=2.0, label="Block axial")
    ax0.plot(time_min, zone_series.block_radial_kpa, color="#2ca02c", linewidth=1.8, linestyle="--", label="Block radial")
    ax0.plot(time_min, zone_series.pellet_axial_kpa, color="#ff7f0e", linewidth=2.0, label="Pellet axial")
    ax0.plot(time_min, zone_series.pellet_radial_kpa, color="#ff7f0e", linewidth=1.8, linestyle="--", label="Pellet radial")
    ax0.set_title("1c Zone-Mean Total Stress")
    ax0.set_xlabel("Time (min)")
    ax0.set_ylabel("Mean total stress magnitude (kPa)")
    ax0.set_xlim(left=0.0)
    ax0.set_ylim(bottom=0.0)
    ax0.grid(True, alpha=0.3)
    ax0.legend(frameon=True, fontsize=8)

    heights = [row["center_mm"] for row in profile_rows]
    ax1.plot([row["porosity"] for row in profile_rows], heights, color="#1f77b4", linewidth=2.0, label="Porosity")
    ax1.plot(
        [row["transport_porosity"] for row in profile_rows],
        heights,
        color="#d62728",
        linewidth=1.8,
        linestyle="--",
        label="Transport porosity",
    )
    ax1.plot(
        [row["saturation"] for row in profile_rows],
        heights,
        color="#9467bd",
        linewidth=1.8,
        linestyle="-.",
        label="Saturation",
    )
    ax1.axhline(1000.0 * block_height_m, color="#666666", linewidth=1.4, linestyle=":", label="Block/pellet interface")
    ax1.set_title("1c Final Vertical Profile")
    ax1.set_xlabel("Field value (-)")
    ax1.set_ylabel("Height from bottom (mm)")
    ax1.set_ylim(0.0, 100.0)
    ax1.grid(True, alpha=0.3)
    ax1.legend(frameon=True, fontsize=8, loc="lower right")

    fig.tight_layout()
    fig.savefig(path, dpi=220)
    plt.close(fig)


def run_uniform_case(
    case: str,
    config: dict,
    native_ogs: Path,
    bridge_ogs: Path,
    out_dir: Path,
    run_cases: bool,
) -> dict:
    native_dir = out_dir / case / "native"
    bridge_dir = out_dir / case / "bridge"
    native_copy = out_dir / case / f"{case}_native_reportplot.prj"
    bridge_copy = out_dir / case / f"{case}_bridge_reportplot.prj"

    write_project_copy(config["native_project"], native_copy, config["native_prefix"], config["output_times_s"])
    write_project_copy(config["bridge_project"], bridge_copy, config["bridge_prefix"], config["output_times_s"])

    if run_cases:
        run_ogs(native_ogs, native_copy, native_dir)
        run_ogs(bridge_ogs, bridge_copy, bridge_dir)

    native_series = collect_boundary_series(native_dir, config["native_prefix"], "Native OGS")
    bridge_series = collect_boundary_series(bridge_dir, config["bridge_prefix"], "MFront OGS")
    native_profile = collect_density_profile(native_dir, config["native_prefix"], config["profile_bins"])
    bridge_profile = collect_density_profile(bridge_dir, config["bridge_prefix"], config["profile_bins"])

    panel_path = out_dir / f"beacon_{case}_report_style_panel.png"
    series_csv = out_dir / f"beacon_{case}_series.csv"
    profile_csv = out_dir / f"beacon_{case}_profile.csv"

    write_boundary_series_csv(series_csv, [native_series, bridge_series])
    report_rows = None
    if case == "1a01":
        report_rows = [
            {"center_mm": h, "mean": v}
            for h, v in zip(config["profile_height_mm"], config["report_profile_kg_m3"], strict=True)
        ]
        plot_1a01_panel(panel_path, [native_series, bridge_series], native_profile, bridge_profile, config)
    elif case == "1b":
        plot_1b_panel(panel_path, [native_series, bridge_series], native_profile, bridge_profile, config)
    write_density_profile_csv(profile_csv, report_rows, native_profile, bridge_profile)

    summary = {
        "case": case,
        "panel": str(panel_path),
        "series_csv": str(series_csv),
        "profile_csv": str(profile_csv),
        "native_final": {
            "time_s": float(native_series.times_s[-1]),
            "axial_kpa": float(native_series.axial_kpa[-1]),
            "radial_kpa": float(native_series.radial_kpa[-1]),
            "dry_density_mean_kg_m3": float(np.mean([row["mean"] for row in native_profile])),
        },
        "bridge_final": {
            "time_s": float(bridge_series.times_s[-1]),
            "axial_kpa": float(bridge_series.axial_kpa[-1]),
            "radial_kpa": float(bridge_series.radial_kpa[-1]),
            "dry_density_mean_kg_m3": float(np.mean([row["mean"] for row in bridge_profile])),
        },
    }
    return summary


def run_1c_case(
    config: dict,
    bridge_ogs: Path,
    out_dir: Path,
    run_cases: bool,
) -> dict:
    case = "1c"
    bridge_dir = out_dir / case / "bridge"
    bridge_copy = out_dir / case / "1c_bridge_reportplot.prj"

    write_project_copy(config["bridge_project"], bridge_copy, config["bridge_prefix"], config["output_times_s"])

    if run_cases:
        run_ogs(bridge_ogs, bridge_copy, bridge_dir)

    zone_series = collect_zone_stress_series(bridge_dir, config["bridge_prefix"], "MFront OGS", config["block_height_m"])
    profile_rows = collect_profile_by_unique_height(bridge_dir, config["bridge_prefix"])

    panel_path = out_dir / "beacon_1c_report_style_panel.png"
    series_csv = out_dir / "beacon_1c_zone_series.csv"
    profile_csv = out_dir / "beacon_1c_profile.csv"

    plot_1c_panel(panel_path, zone_series, profile_rows, config["block_height_m"])
    write_zone_series_csv(series_csv, zone_series)
    write_zone_profile_csv(profile_csv, profile_rows)

    summary = {
        "case": case,
        "panel": str(panel_path),
        "series_csv": str(series_csv),
        "profile_csv": str(profile_csv),
        "bridge_final": {
            "time_s": float(zone_series.times_s[-1]),
            "block_axial_kpa": float(zone_series.block_axial_kpa[-1]),
            "block_radial_kpa": float(zone_series.block_radial_kpa[-1]),
            "pellet_axial_kpa": float(zone_series.pellet_axial_kpa[-1]),
            "pellet_radial_kpa": float(zone_series.pellet_radial_kpa[-1]),
        },
    }
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cases", nargs="+", choices=["1a01", "1b", "1c"], default=["1a01", "1b", "1c"])
    parser.add_argument("--native-ogs", type=Path, default=NATIVE_DEFAULT)
    parser.add_argument("--bridge-ogs", type=Path, default=MFRONT_DEFAULT)
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=HERE / "BEACON_report_figures",
    )
    parser.add_argument("--skip-run", action="store_true")
    args = parser.parse_args()

    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    summaries = []
    for case in args.cases:
        if case in ("1a01", "1b"):
            summaries.append(
                run_uniform_case(
                    case,
                    CASE_CONFIGS[case],
                    args.native_ogs.resolve(),
                    args.bridge_ogs.resolve(),
                    out_dir,
                    run_cases=not args.skip_run,
                )
            )
        else:
            summaries.append(
                run_1c_case(
                    CASE_CONFIGS["1c"],
                    args.bridge_ogs.resolve(),
                    out_dir,
                    run_cases=not args.skip_run,
                )
            )

    summary_path = out_dir / "beacon_report_plot_summary.json"
    summary_path.write_text(json.dumps(summaries, indent=2) + "\n")
    print(json.dumps({"summary": str(summary_path), "cases": summaries}, indent=2))


if __name__ == "__main__":
    main()
