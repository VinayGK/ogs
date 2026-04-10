#!/usr/bin/env python3
"""Probe 1b dry-density/multiplier calibration sensitivity for native vs MFront.

This script does not overwrite committed project files. It writes temporary
copies, runs OGS, extracts final stresses/density, and stores a compact
CSV+JSON report under _outputs/.
"""

from __future__ import annotations

import csv
import json
import math
import shutil
import subprocess
import tempfile
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy


ROOT = Path(__file__).resolve().parent
OUT_ROOT = ROOT / "_outputs" / "beacon_1b_ddcurve_probe"
ANCHORS_CALIB_CSV = ROOT / "ANCHORS_MS33_ModelI" / "villar_dense_dd_calibration.csv"

NATIVE_OGS = Path("/Users/vinaykumar/git/build/release-native-transition-mfront/bin/ogs")
MFRONT_OGS = Path("/Users/vinaykumar/git/build/release-mfront-tpm/bin/ogs")

NATIVE_PROJECT = ROOT / "beacon_1b_unstructured_batch.prj"
MFRONT_PROJECT = ROOT / "beacon_1b_notebook_mcc_unstructured_batch.prj"

HAMAKER_REFERENCE_J = 5.1e-21
TARGET_DD_KG_M3 = 1520.0
RHO_SOLID_KG_M3 = 2780.0


@dataclass(frozen=True)
class CaseResult:
    case_id: str
    implementation: str
    run_status: str
    dry_density_target_kg_m3: float
    multiplier_used: float
    hamaker_constant_j: float
    axial_sigma_kpa: float
    radial_sigma_kpa: float
    axial_swelling_stress_kpa: float
    radial_swelling_stress_kpa: float
    dry_density_mean_kg_m3: float
    note: str


def load_multiplier_curve(path: Path) -> tuple[np.ndarray, np.ndarray]:
    with path.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    dd = np.array([float(r["dry_density_kg_m3"]) for r in rows], dtype=float)
    mult = np.array([float(r["vdw_multiplier"]) for r in rows], dtype=float)
    order = np.argsort(dd)
    return dd[order], mult[order]


def interpolate_multiplier(dd: float, dd_curve: np.ndarray, mult_curve: np.ndarray) -> float:
    return float(np.interp(dd, dd_curve, mult_curve, left=mult_curve[0], right=mult_curve[-1]))


def set_parameter_value(root: ET.Element, name: str, value: float) -> bool:
    node = root.find(f"./parameters/parameter[name='{name}']/value")
    if node is None:
        return False
    node.text = f"{value:.16g}"
    return True


def set_native_hamaker(root: ET.Element, hamaker_j: float) -> None:
    node = root.find("./processes/process/potential_exchange/hamaker_constant")
    if node is None:
        raise RuntimeError("Missing potential_exchange/hamaker_constant in native project.")
    node.text = f"{hamaker_j:.16g}"


def set_mfront_hamaker(root: ET.Element, hamaker_j: float) -> None:
    if not set_parameter_value(root, "HamakerConstant", hamaker_j):
        raise RuntimeError("Missing HamakerConstant parameter in MFront project.")


def absolutize_mesh_and_geometry(root: ET.Element, source_dir: Path) -> None:
    mesh = root.find("./mesh")
    geometry = root.find("./geometry")
    if mesh is not None and mesh.text:
        mesh.text = str((source_dir / mesh.text.strip()).resolve())
    if geometry is not None and geometry.text:
        geometry.text = str((source_dir / geometry.text.strip()).resolve())


def set_output_prefix(root: ET.Element, prefix: str) -> None:
    node = root.find("./time_loop/output/prefix")
    if node is None:
        raise RuntimeError("Could not find output prefix.")
    node.text = prefix


def set_solver_controls(root: ET.Element) -> None:
    for node in root.findall("./nonlinear_solvers/nonlinear_solver/max_iter"):
        node.text = "60"
    for node in root.findall(".//convergence_criterion/abstols"):
        if node.text is None:
            continue
        parts = node.text.strip().split()
        if not parts:
            continue
        parts[0] = "5e-7"
        node.text = " ".join(parts)


def build_project_copy(
    *,
    source: Path,
    target: Path,
    prefix: str,
    implementation: str,
    multiplier: float,
    apply_dd_target: bool,
    apply_micro_activation: bool,
) -> None:
    root = ET.parse(source).getroot()
    absolutize_mesh_and_geometry(root, source.parent)
    set_output_prefix(root, prefix)
    set_solver_controls(root)

    if apply_dd_target:
        phi0 = 1.0 - TARGET_DD_KG_M3 / RHO_SOLID_KG_M3
        set_parameter_value(root, "phi0", phi0)
        # Keep transport porosity below total porosity.
        set_parameter_value(root, "phi_tr0", max(0.0, phi0 - 0.1))

    hamaker_j = HAMAKER_REFERENCE_J * multiplier
    if implementation == "native":
        set_native_hamaker(root, hamaker_j)
    elif implementation == "mfront":
        set_mfront_hamaker(root, hamaker_j)
        if apply_micro_activation:
            # Switch on the same micro branch used in dense dd calibration.
            set_parameter_value(root, "NotebookSwellingSlope", 0.1)
            set_parameter_value(root, "NotebookSaturationMode", 1.0)
            set_parameter_value(root, "NotebookLocalSolveMode", 0.0)
            set_parameter_value(root, "MicroPotentialConvention", 1.0)
            set_parameter_value(root, "SpecificSurface", 523.0)
            set_parameter_value(root, "n_l0", 0.00153385355)
            set_parameter_value(root, "rho_lR0", 2276.031917690513)
    else:
        raise ValueError(implementation)

    target.parent.mkdir(parents=True, exist_ok=True)
    ET.ElementTree(root).write(target, encoding="ISO-8859-1", xml_declaration=True)


def run_ogs(ogs_bin: Path, project: Path, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    result = subprocess.run(
        [str(ogs_bin), "-o", str(output_dir), str(project)],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"ogs failed for {project.name}\n"
            f"stdout:\n{result.stdout}\n"
            f"stderr:\n{result.stderr}"
        )


def parse_output_time(path: Path) -> float:
    marker = "_t_"
    if marker not in path.stem:
        return -math.inf
    return float(path.stem.split(marker, 1)[1])


def latest_vtu(output_dir: Path, prefix: str) -> Path:
    files = sorted(output_dir.glob(f"{prefix}_t_*.vtu"), key=parse_output_time)
    if not files:
        raise FileNotFoundError(f"No VTU files found for prefix {prefix}")
    return files[-1]


def load_vtu(path: Path) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(path))
    reader.Update()
    grid = reader.GetOutput()
    points = np.array([grid.GetPoint(i) for i in range(grid.GetNumberOfPoints())])
    arrays: dict[str, np.ndarray] = {}
    for i in range(grid.GetPointData().GetNumberOfArrays()):
        arr = grid.GetPointData().GetArray(i)
        name = grid.GetPointData().GetArrayName(i)
        if arr is None or name is None:
            continue
        arrays[name] = vtk_to_numpy(arr)
    return points, arrays


def boundary_mean(points: np.ndarray, values: np.ndarray, axis: int, side: str) -> np.ndarray:
    target = points[:, axis].max() if side == "max" else points[:, axis].min()
    mask = np.isclose(points[:, axis], target, atol=1e-10)
    return np.asarray(values[mask]).mean(axis=0)


def evaluate_metrics(vtu_path: Path) -> tuple[float, float, float, float, float]:
    points, arrays = load_vtu(vtu_path)
    sigma = arrays["sigma"]
    swelling = arrays.get("swelling_stress")
    dry_density = arrays.get("dry_density_solid")

    sigma_top = boundary_mean(points, sigma, axis=1, side="max")
    sigma_side = boundary_mean(points, sigma, axis=0, side="max")

    if swelling is None:
        sw_top = np.array([math.nan, math.nan])
        sw_side = np.array([math.nan, math.nan])
    else:
        sw_top = boundary_mean(points, swelling, axis=1, side="max")
        sw_side = boundary_mean(points, swelling, axis=0, side="max")

    dd_mean = float(np.asarray(dry_density).mean()) if dry_density is not None else math.nan
    return (
        abs(float(sigma_top[1])) / 1e3,
        abs(float(sigma_side[0])) / 1e3,
        abs(float(sw_top[1])) / 1e3,
        abs(float(sw_side[0])) / 1e3,
        dd_mean,
    )


def run_case(
    *,
    case_id: str,
    implementation: str,
    source_project: Path,
    ogs_bin: Path,
    multiplier: float,
    apply_dd_target: bool,
    apply_micro_activation: bool,
    note: str,
    tmpdir: Path,
) -> CaseResult:
    prefix = f"{case_id}_{implementation}"
    project_copy = tmpdir / f"{prefix}.prj"
    out_dir = tmpdir / f"out_{prefix}"
    hamaker_j = HAMAKER_REFERENCE_J * multiplier

    build_project_copy(
        source=source_project,
        target=project_copy,
        prefix=prefix,
        implementation=implementation,
        multiplier=multiplier,
        apply_dd_target=apply_dd_target,
        apply_micro_activation=apply_micro_activation,
    )

    try:
        run_ogs(ogs_bin, project_copy, out_dir)
        ax, rad, ax_sw, rad_sw, dd = evaluate_metrics(latest_vtu(out_dir, prefix))
        return CaseResult(
            case_id=case_id,
            implementation=implementation,
            run_status="success",
            dry_density_target_kg_m3=TARGET_DD_KG_M3 if apply_dd_target else math.nan,
            multiplier_used=multiplier,
            hamaker_constant_j=hamaker_j,
            axial_sigma_kpa=ax,
            radial_sigma_kpa=rad,
            axial_swelling_stress_kpa=ax_sw,
            radial_swelling_stress_kpa=rad_sw,
            dry_density_mean_kg_m3=dd,
            note=note,
        )
    except Exception as exc:
        return CaseResult(
            case_id=case_id,
            implementation=implementation,
            run_status="failed",
            dry_density_target_kg_m3=TARGET_DD_KG_M3 if apply_dd_target else math.nan,
            multiplier_used=multiplier,
            hamaker_constant_j=hamaker_j,
            axial_sigma_kpa=math.nan,
            radial_sigma_kpa=math.nan,
            axial_swelling_stress_kpa=math.nan,
            radial_swelling_stress_kpa=math.nan,
            dry_density_mean_kg_m3=math.nan,
            note=f"{note}; error={str(exc).splitlines()[-1][:220]}",
        )


def write_rows_csv(rows: list[CaseResult], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    keys = list(CaseResult.__dataclass_fields__.keys())
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        for row in rows:
            writer.writerow(row.__dict__)


def main() -> None:
    dd_curve, mult_curve = load_multiplier_curve(ANCHORS_CALIB_CSV)
    interpolated_multiplier = interpolate_multiplier(TARGET_DD_KG_M3, dd_curve, mult_curve)

    tmpdir = Path(tempfile.mkdtemp(prefix="beacon_1b_ddcurve_probe_"))
    rows: list[CaseResult] = []
    try:
        rows.append(
            run_case(
                case_id="baseline_1b",
                implementation="native",
                source_project=NATIVE_PROJECT,
                ogs_bin=NATIVE_OGS,
                multiplier=1.0,
                apply_dd_target=False,
                apply_micro_activation=False,
                note="Committed native 1b deck (qualitative report comparison).",
                tmpdir=tmpdir,
            )
        )
        rows.append(
            run_case(
                case_id="baseline_1b",
                implementation="mfront",
                source_project=MFRONT_PROJECT,
                ogs_bin=MFRONT_OGS,
                multiplier=1.0,
                apply_dd_target=False,
                apply_micro_activation=False,
                note="Committed MFront 1b deck (qualitative report comparison).",
                tmpdir=tmpdir,
            )
        )
        rows.append(
            run_case(
                case_id="dd_target_probe",
                implementation="native",
                source_project=NATIVE_PROJECT,
                ogs_bin=NATIVE_OGS,
                multiplier=interpolated_multiplier,
                apply_dd_target=True,
                apply_micro_activation=False,
                note=(
                    "Native 1b with DD target + interpolated multiplier from dense Villar curve; "
                    "native 1b swelling remains on saturation-dependent law."
                ),
                tmpdir=tmpdir,
            )
        )
        rows.append(
            run_case(
                case_id="dd_target_probe",
                implementation="mfront",
                source_project=MFRONT_PROJECT,
                ogs_bin=MFRONT_OGS,
                multiplier=interpolated_multiplier,
                apply_dd_target=True,
                apply_micro_activation=True,
                note=(
                    "MFront 1b with DD target + micro activation + interpolated multiplier "
                    "from dense Villar curve."
                ),
                tmpdir=tmpdir,
            )
        )
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)

    out_csv = OUT_ROOT / "beacon_1b_ddcurve_probe_rows.csv"
    out_json = OUT_ROOT / "beacon_1b_ddcurve_probe_summary.json"
    write_rows_csv(rows, out_csv)

    summary = {
        "target_dry_density_kg_m3": TARGET_DD_KG_M3,
        "rho_solid_kg_m3": RHO_SOLID_KG_M3,
        "interpolated_multiplier_from_dense_villar_curve": interpolated_multiplier,
        "reference_curve_csv": str(ANCHORS_CALIB_CSV),
        "rows": [r.__dict__ for r in rows],
        "notes": [
            "BEACON D5.1.1 provides a qualitative 1b stress target ('nonzero and stabilized after about 500 days') but no tabulated plateau value.",
            "The DD-target probe therefore evaluates sensitivity to DD-corrected setup plus multiplier transfer from the dense Villar calibration.",
            "Native and MFront 1b decks are not equivalent in swelling-source formulation by default.",
        ],
    }
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    print(f"Wrote: {out_csv}")
    print(f"Wrote: {out_json}")


if __name__ == "__main__":
    main()
