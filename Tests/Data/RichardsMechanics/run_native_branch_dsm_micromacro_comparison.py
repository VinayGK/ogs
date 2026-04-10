#!/usr/bin/env python3
"""Compare native dsm_micromacro branch and MFront bridge outputs on BEACON cases.

The script executes paired projects, extracts final-state stress/density
metrics from VTK outputs, computes field-wise maximum absolute differences, and
writes CSV/JSON/PNG artifacts for parity tracking.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import subprocess
import tempfile
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy

os.environ.setdefault("MPLCONFIGDIR", tempfile.mkdtemp(prefix="mplconfig-"))
os.environ.setdefault("XDG_CACHE_HOME", tempfile.mkdtemp(prefix="xdg-cache-"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


HERE = Path(__file__).resolve().parent
DEFAULT_NATIVE_SOURCE = Path("/Users/vinaykumar/git/ogs-native-dsm-transition")
DEFAULT_MFRONT_SOURCE = HERE.parents[2]
DEFAULT_NATIVE_OGS = Path("/Users/vinaykumar/git/build/release-native-beacon/bin/ogs")
DEFAULT_MFRONT_OGS = Path("/Users/vinaykumar/git/build/release-mfront-tpm/bin/ogs")


@dataclass(frozen=True)
class CaseConfig:
    native_project: str
    mfront_project: str


CASE_CONFIGS: dict[str, CaseConfig] = {
    "1a01_smoke": CaseConfig(
        native_project="beacon_1a01_dsm_micromacro_smoke.prj",
        mfront_project="beacon_1a01_dsm_micromacro_mcc_bridge.prj",
    ),
    "1b_smoke": CaseConfig(
        native_project="beacon_1b_dsm_micromacro_smoke.prj",
        mfront_project="beacon_1b_dsm_micromacro_mcc_bridge.prj",
    ),
    "1c_smoke": CaseConfig(
        native_project="beacon_1c_dsm_micromacro_smoke.prj",
        mfront_project="beacon_1c_dsm_micromacro_mcc_bridge.prj",
    ),
    "1a01_inflow": CaseConfig(
        native_project="beacon_1a01_dsm_micromacro_inflow.prj",
        mfront_project="beacon_1a01_dsm_micromacro_mcc_inflow_bridge.prj",
    ),
}


def run_ogs(executable: Path, project: Path, output_dir: Path) -> None:
    """Run OGS for one project and capture logs for reproducible failure diagnosis."""
    output_dir.mkdir(parents=True, exist_ok=True)
    result = subprocess.run(
        [str(executable), "-o", str(output_dir), str(project)],
        cwd=project.parent,
        text=True,
        capture_output=True,
        check=False,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"ogs failed for {project} with {executable}\n"
            f"stdout:\n{result.stdout}\n"
            f"stderr:\n{result.stderr}"
        )


def project_output_prefix(project: Path) -> str:
    root = ET.parse(project).getroot()
    prefix = root.findtext("./time_loop/output/prefix")
    if prefix is None:
        raise RuntimeError(f"No output prefix in {project}")
    return prefix.strip()


def parse_time_from_vtu(path: Path) -> float:
    marker = "_t_"
    stem = path.stem
    if marker not in stem:
        raise RuntimeError(f"Cannot parse output time from {path.name}")
    return float(stem.split(marker, 1)[1])


def latest_vtu(output_dir: Path, prefix: str) -> Path:
    files = sorted(output_dir.glob(f"{prefix}_t_*.vtu"), key=parse_time_from_vtu)
    if not files:
        raise FileNotFoundError(f"No VTU outputs in {output_dir} for prefix {prefix}")
    return files[-1]


def load_vtu(path: Path) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(path))
    reader.Update()
    grid = reader.GetOutput()
    points = np.array([grid.GetPoint(i) for i in range(grid.GetNumberOfPoints())])
    arrays: dict[str, np.ndarray] = {}
    for i in range(grid.GetPointData().GetNumberOfArrays()):
        array = grid.GetPointData().GetArray(i)
        if array is None:
            continue
        name = grid.GetPointData().GetArrayName(i)
        if name is None:
            continue
        arrays[name] = vtk_to_numpy(array)
    return points, arrays


def boundary_mean(points: np.ndarray, values: np.ndarray, axis: int, side: str) -> np.ndarray:
    boundary_value = points[:, axis].max() if side == "max" else points[:, axis].min()
    mask = np.isclose(points[:, axis], boundary_value, atol=1e-10)
    return np.asarray(values[mask]).mean(axis=0)


def component_or_scalar(arr: np.ndarray, preferred_index: int) -> float:
    values = np.asarray(arr).reshape(-1)
    if values.size == 0:
        return float("nan")
    idx = preferred_index if preferred_index < values.size else 0
    return float(values[idx])


def evaluate_metrics(points: np.ndarray, arrays: dict[str, np.ndarray]) -> dict[str, float]:
    """Extract comparable boundary stress and scalar-field means from final VTU data."""
    sigma = arrays.get("sigma")
    swelling = arrays.get("swelling_stress")
    dry_density = arrays.get("dry_density_solid")

    sigma_top = boundary_mean(points, sigma, axis=1, side="max") if sigma is not None else np.array([np.nan, np.nan])
    sigma_side = boundary_mean(points, sigma, axis=0, side="max") if sigma is not None else np.array([np.nan, np.nan])
    swell_top = (
        boundary_mean(points, swelling, axis=1, side="max")
        if swelling is not None
        else np.array([np.nan, np.nan])
    )
    swell_side = (
        boundary_mean(points, swelling, axis=0, side="max")
        if swelling is not None
        else np.array([np.nan, np.nan])
    )

    metrics: dict[str, float] = {
        "axial_sigma_kpa": abs(component_or_scalar(sigma_top, 1)) / 1000.0,
        "radial_sigma_kpa": abs(component_or_scalar(sigma_side, 0)) / 1000.0,
        "axial_swelling_stress_kpa": abs(component_or_scalar(swell_top, 1)) / 1000.0,
        "radial_swelling_stress_kpa": abs(component_or_scalar(swell_side, 0)) / 1000.0,
        "dry_density_mean_kg_m3": float(np.asarray(dry_density).mean()) if dry_density is not None else float("nan"),
    }

    optional_scalar_fields = [
        "pressure",
        "saturation",
        "micro_pressure",
        "micro_saturation",
        "micro_water_content",
        "micro_porosity",
        "micro_exchange_source",
        "phi_m",
        "phi_M",
        "n_l",
    ]
    for field in optional_scalar_fields:
        data = arrays.get(field)
        if data is not None:
            metrics[f"{field}_mean"] = float(np.asarray(data).mean())
    return metrics


def max_abs_diff(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.abs(np.asarray(a) - np.asarray(b)).max())


def compare_fields(
    native_points: np.ndarray,
    native_data: dict[str, np.ndarray],
    mfront_points: np.ndarray,
    mfront_data: dict[str, np.ndarray],
) -> dict[str, float]:
    """Return max-abs field differences on shared point ordering."""
    if not np.allclose(native_points, mfront_points):
        raise RuntimeError("Point coordinates differ between native and MFront outputs.")

    diffs: dict[str, float] = {}
    candidate_fields = [
        "pressure",
        "saturation",
        "sigma",
        "swelling_stress",
        "dry_density_solid",
    ]
    for field in candidate_fields:
        if field in native_data and field in mfront_data:
            diffs[field] = max_abs_diff(native_data[field], mfront_data[field])
    return diffs


def write_csv(rows: list[dict[str, object]], path: Path) -> None:
    columns = [
        "case",
        "native_axial_sigma_kpa",
        "mfront_axial_sigma_kpa",
        "native_radial_sigma_kpa",
        "mfront_radial_sigma_kpa",
        "native_dry_density_mean_kg_m3",
        "mfront_dry_density_mean_kg_m3",
        "native_pressure_mean",
        "mfront_pressure_mean",
        "max_abs_diff_pressure",
        "max_abs_diff_sigma",
        "max_abs_diff_swelling_stress",
    ]
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def make_stress_plot(rows: list[dict[str, object]], output_png: Path) -> None:
    cases = [str(r["case"]) for r in rows]
    x = np.arange(len(cases))
    width = 0.35

    native_axial = np.array([float(r["native_axial_sigma_kpa"]) for r in rows])
    mfront_axial = np.array([float(r["mfront_axial_sigma_kpa"]) for r in rows])
    native_radial = np.array([float(r["native_radial_sigma_kpa"]) for r in rows])
    mfront_radial = np.array([float(r["mfront_radial_sigma_kpa"]) for r in rows])

    fig, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=True)

    axes[0].bar(x - width / 2, native_axial, width=width, label="native branch")
    axes[0].bar(x + width / 2, mfront_axial, width=width, label="MFront bridge branch")
    axes[0].set_ylabel("axial stress [kPa]")
    axes[0].grid(axis="y", alpha=0.3)
    axes[0].legend()

    axes[1].bar(x - width / 2, native_radial, width=width, label="native branch")
    axes[1].bar(x + width / 2, mfront_radial, width=width, label="MFront bridge branch")
    axes[1].set_ylabel("radial stress [kPa]")
    axes[1].set_xticks(x, cases, rotation=15, ha="right")
    axes[1].grid(axis="y", alpha=0.3)

    fig.tight_layout()
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=200)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--native-source-dir", type=Path, default=DEFAULT_NATIVE_SOURCE)
    parser.add_argument("--mfront-source-dir", type=Path, default=DEFAULT_MFRONT_SOURCE)
    parser.add_argument("--native-ogs", type=Path, default=DEFAULT_NATIVE_OGS)
    parser.add_argument("--mfront-ogs", type=Path, default=DEFAULT_MFRONT_OGS)
    parser.add_argument(
        "--output-root",
        type=Path,
        required=True,
        help="Transient runtime output folder for VTU files.",
    )
    parser.add_argument("--json-out", type=Path, required=True)
    parser.add_argument("--csv-out", type=Path, required=True)
    parser.add_argument("--plot-out", type=Path, required=True)
    args = parser.parse_args()

    native_data_dir = args.native_source_dir / "Tests" / "Data" / "RichardsMechanics"
    mfront_data_dir = args.mfront_source_dir / "Tests" / "Data" / "RichardsMechanics"

    rows: list[dict[str, object]] = []
    report: dict[str, object] = {
        "native_source_dir": str(args.native_source_dir.resolve()),
        "mfront_source_dir": str(args.mfront_source_dir.resolve()),
        "native_ogs": str(args.native_ogs.resolve()),
        "mfront_ogs": str(args.mfront_ogs.resolve()),
        "cases": {},
    }

    for case_name, case in CASE_CONFIGS.items():
        native_project = native_data_dir / case.native_project
        mfront_project = mfront_data_dir / case.mfront_project
        native_run_dir = args.output_root / case_name / "native"
        mfront_run_dir = args.output_root / case_name / "mfront"

        run_ogs(args.native_ogs, native_project, native_run_dir)
        run_ogs(args.mfront_ogs, mfront_project, mfront_run_dir)

        native_prefix = project_output_prefix(native_project)
        mfront_prefix = project_output_prefix(mfront_project)
        native_vtu = latest_vtu(native_run_dir, native_prefix)
        mfront_vtu = latest_vtu(mfront_run_dir, mfront_prefix)

        native_points, native_arrays = load_vtu(native_vtu)
        mfront_points, mfront_arrays = load_vtu(mfront_vtu)

        native_metrics = evaluate_metrics(native_points, native_arrays)
        mfront_metrics = evaluate_metrics(mfront_points, mfront_arrays)
        diffs = compare_fields(native_points, native_arrays, mfront_points, mfront_arrays)

        report_case = {
            "native_project": str(native_project),
            "mfront_project": str(mfront_project),
            "native_vtu": str(native_vtu),
            "mfront_vtu": str(mfront_vtu),
            "native_metrics": native_metrics,
            "mfront_metrics": mfront_metrics,
            "max_abs_diff": diffs,
        }
        casted_cases = report["cases"]
        assert isinstance(casted_cases, dict)
        casted_cases[case_name] = report_case

        rows.append(
            {
                "case": case_name,
                "native_axial_sigma_kpa": native_metrics["axial_sigma_kpa"],
                "mfront_axial_sigma_kpa": mfront_metrics["axial_sigma_kpa"],
                "native_radial_sigma_kpa": native_metrics["radial_sigma_kpa"],
                "mfront_radial_sigma_kpa": mfront_metrics["radial_sigma_kpa"],
                "native_dry_density_mean_kg_m3": native_metrics["dry_density_mean_kg_m3"],
                "mfront_dry_density_mean_kg_m3": mfront_metrics["dry_density_mean_kg_m3"],
                "native_pressure_mean": native_metrics.get("pressure_mean", math.nan),
                "mfront_pressure_mean": mfront_metrics.get("pressure_mean", math.nan),
                "max_abs_diff_pressure": diffs.get("pressure", math.nan),
                "max_abs_diff_sigma": diffs.get("sigma", math.nan),
                "max_abs_diff_swelling_stress": diffs.get("swelling_stress", math.nan),
            }
        )

    args.json_out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_csv(rows, args.csv_out)
    make_stress_plot(rows, args.plot_out)

    print(json.dumps({"json_out": str(args.json_out), "csv_out": str(args.csv_out), "plot_out": str(args.plot_out)}, indent=2))


if __name__ == "__main__":
    main()
