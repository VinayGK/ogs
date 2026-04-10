#!/usr/bin/env python3
"""Run the dry-density calibrated native-vs-MFront benchmark comparison.

The script reuses the ANCHORS calibration curves, reruns BEACON with
interpolated calibrated multipliers, and writes combined CSV/JSON summaries
plus pairwise delta tables for the transition note.
"""

from __future__ import annotations

import argparse
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
ANCHORS_DIR = ROOT / "ANCHORS_MS33_ModelI"
MFRONT_CALIBRATION_CSV = ANCHORS_DIR / "villar_dense_dd_calibration.csv"
NATIVE_CALIBRATION_CSV = ANCHORS_DIR / "villar_dense_dd_native_notebook_calibration.csv"
HAMAKER_REFERENCE_J = 5.1e-21
OUTPUT_ROOT = ROOT / "_outputs" / "calibrated_native_mfront_comparison"


@dataclass(frozen=True)
class BeaconCase:
    case_id: str
    native_project: Path
    mfront_project: Path


BEACON_CASES = [
    BeaconCase(
        case_id="1a01",
        native_project=ROOT / "beacon_1a01_inflow_unstructured_batch.prj",
        mfront_project=ROOT / "beacon_1a01_notebook_mcc_inflow_unstructured_batch.prj",
    ),
    BeaconCase(
        case_id="1b",
        native_project=ROOT / "beacon_1b_unstructured_batch.prj",
        mfront_project=ROOT / "beacon_1b_notebook_mcc_unstructured_batch.prj",
    ),
]


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the comparison workflow."""
    parser = argparse.ArgumentParser(
        description=(
            "Run BEACON native-vs-MFront comparison with dry-density calibrated "
            "vdW multipliers, then merge with ANCHORS DD calibration rows."
        )
    )
    parser.add_argument(
        "--native-ogs",
        type=Path,
        default=Path("/Users/vinaykumar/git/build/release-native-transition-mfront/bin/ogs"),
    )
    parser.add_argument(
        "--mfront-ogs",
        type=Path,
        default=Path("/Users/vinaykumar/git/build/release-mfront-tpm/bin/ogs"),
    )
    parser.add_argument(
        "--out-csv",
        type=Path,
        default=OUTPUT_ROOT / "calibrated_benchmark_runs_summary.csv",
    )
    parser.add_argument(
        "--out-json",
        type=Path,
        default=OUTPUT_ROOT / "calibrated_benchmark_runs_summary.json",
    )
    parser.add_argument(
        "--out-deltas-csv",
        type=Path,
        default=OUTPUT_ROOT / "calibrated_benchmark_pairwise_deltas.csv",
    )
    parser.add_argument(
        "--skip-beacon-runs",
        action="store_true",
        help="Only merge existing ANCHORS calibration rows.",
    )
    return parser.parse_args()


def parse_csv_rows(path: Path) -> list[dict[str, str]]:
    """Load a CSV file into a list of row dictionaries."""
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def load_multiplier_curve(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Load the dry-density calibration curve as sorted arrays."""
    rows = parse_csv_rows(path)
    dd = np.array([float(r["dry_density_kg_m3"]) for r in rows], dtype=float)
    mult = np.array([float(r["vdw_multiplier"]) for r in rows], dtype=float)
    order = np.argsort(dd)
    return dd[order], mult[order]


def interpolate_multiplier(dd: float, dd_curve: np.ndarray, mult_curve: np.ndarray) -> float:
    """Interpolate the effective multiplier for one dry density."""
    return float(np.interp(dd, dd_curve, mult_curve, left=mult_curve[0], right=mult_curve[-1]))


def find_parameter_value(root: ET.Element, name: str) -> str | None:
    node = root.find(f"./parameters/parameter[name='{name}']/value")
    if node is None or node.text is None:
        return None
    return node.text.strip()


def extract_dd_from_project(project: Path) -> float:
    """Infer the dry density from a project file's solid density and phi0."""
    root = ET.parse(project).getroot()
    phi0_text = find_parameter_value(root, "phi0")
    if phi0_text is None:
        raise RuntimeError(f"Could not find parameter phi0 in {project}")
    rho_text = root.findtext("./media/medium/phases/phase[type='Solid']/properties/property[name='density']/value")
    if rho_text is None:
        raise RuntimeError(f"Could not find solid density in {project}")
    phi0 = float(phi0_text)
    rho_s = float(rho_text.strip())
    return rho_s * (1.0 - phi0)


def absolutize_mesh_and_geometry(root: ET.Element, base_dir: Path) -> None:
    """Resolve mesh and geometry paths relative to the source project."""
    mesh = root.find("./mesh")
    if mesh is not None and mesh.text:
        mesh.text = str((base_dir / mesh.text.strip()).resolve())
    geometry = root.find("./geometry")
    if geometry is not None and geometry.text:
        geometry.text = str((base_dir / geometry.text.strip()).resolve())


def set_output_prefix(root: ET.Element, prefix: str) -> None:
    """Rewrite the VTK output prefix in a copied project file."""
    node = root.find("./time_loop/output/prefix")
    if node is None:
        raise RuntimeError("Could not find <time_loop><output><prefix> in project.")
    node.text = prefix


def set_biot_coefficient(root: ET.Element, value: float) -> None:
    """Force all medium biot_coefficient properties to a consistent scalar."""
    found = False
    for prop in root.findall(".//property"):
        name = (prop.findtext("name") or "").strip()
        if name != "biot_coefficient":
            continue
        val = prop.find("value")
        if val is None:
            val = ET.SubElement(prop, "value")
        val.text = f"{value:.16g}"
        found = True
    if not found:
        raise RuntimeError("Could not find any biot_coefficient property in project.")


def set_nonlinear_max_iter(root: ET.Element, value: int) -> None:
    """Set Newton max_iter in copied benchmark projects."""
    nodes = root.findall("./nonlinear_solvers/nonlinear_solver/max_iter")
    if not nodes:
        raise RuntimeError("Could not find nonlinear solver max_iter in project.")
    for node in nodes:
        node.text = str(value)


def relax_first_component_abstol(root: ET.Element, value: float) -> None:
    """Relax first PerComponentDeltaX absolute tolerance in copied decks."""
    nodes = root.findall(".//convergence_criterion/abstols")
    if not nodes:
        raise RuntimeError("Could not find convergence_criterion/abstols in project.")
    for node in nodes:
        if node.text is None:
            continue
        parts = node.text.strip().split()
        if not parts:
            continue
        parts[0] = f"{value:.16g}"
        node.text = " ".join(parts)


def set_native_hamaker(root: ET.Element, hamaker_j: float) -> None:
    """Inject the calibrated Hamaker constant into a native project copy."""
    node = root.find("./processes/process/potential_exchange/hamaker_constant")
    if node is None:
        raise RuntimeError("Could not find potential_exchange/hamaker_constant in native project.")
    node.text = f"{hamaker_j:.16g}"


def set_mfront_hamaker(root: ET.Element, hamaker_j: float) -> None:
    """Inject the calibrated Hamaker constant into an MFront project copy."""
    node = root.find("./parameters/parameter[name='HamakerConstant']/value")
    if node is None:
        raise RuntimeError("Could not find parameter HamakerConstant in MFront project.")
    node.text = f"{hamaker_j:.16g}"


def write_project_copy(
    source: Path,
    target: Path,
    prefix: str,
    hamaker_j: float,
    implementation: str,
) -> None:
    """Clone a project and patch its output prefix and vdW parameter."""
    root = ET.parse(source).getroot()
    absolutize_mesh_and_geometry(root, source.parent)
    set_output_prefix(root, prefix)
    set_biot_coefficient(root, 1.0)
    # Current-density constitutive updates in the dense benchmark shell can
    # require more global Newton iterations on early inflow steps.
    set_nonlinear_max_iter(root, 60)
    relax_first_component_abstol(root, 5e-7)
    if implementation == "native":
        set_native_hamaker(root, hamaker_j)
    elif implementation == "mfront":
        set_mfront_hamaker(root, hamaker_j)
    else:
        raise ValueError(implementation)
    target.parent.mkdir(parents=True, exist_ok=True)
    ET.ElementTree(root).write(target, encoding="ISO-8859-1", xml_declaration=True)


def run_ogs(ogs_bin: Path, project: Path, output_dir: Path) -> None:
    """Run OGS for one temporary project and fail with captured logs."""
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
            f"ogs failed for {project.name} with {ogs_bin}\n"
            f"stdout:\n{result.stdout}\n"
            f"stderr:\n{result.stderr}"
        )


def summarize_error(exc: Exception) -> str:
    """Extract a compact, informative error line from a raised exception."""
    text = str(exc)
    lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
    # Prefer explicit OGS physics-failure lines.
    for ln in reversed(lines):
        low = ln.lower()
        if "biot-coefficient" in low or ("porosity" in low and low.startswith("error:")):
            return ln[:500]
    # Then prefer explicit OGS error lines.
    for ln in reversed(lines):
        low = ln.lower()
        if low.startswith("error:") or "terminated with error" in low:
            return ln[:500]
    # Fall back to any line mentioning stderr/stdout details.
    for ln in reversed(lines):
        if "stderr" in ln.lower() or "stdout" in ln.lower():
            continue
        return ln[:500]
    return text.strip()[:500]


def parse_output_time(path: Path) -> float:
    marker = "_t_"
    if marker not in path.stem:
        return -math.inf
    return float(path.stem.split(marker, 1)[1])


def latest_vtu(output_dir: Path, prefix: str) -> Path:
    files = sorted(output_dir.glob(f"{prefix}_t_*.vtu"), key=parse_output_time)
    if not files:
        raise FileNotFoundError(f"No VTU output found for prefix {prefix} in {output_dir}")
    return files[-1]


def load_vtu(path: Path) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(path))
    reader.Update()
    grid = reader.GetOutput()
    points = np.array([grid.GetPoint(i) for i in range(grid.GetNumberOfPoints())])
    data: dict[str, np.ndarray] = {}
    for i in range(grid.GetPointData().GetNumberOfArrays()):
        arr = grid.GetPointData().GetArray(i)
        name = grid.GetPointData().GetArrayName(i)
        if arr is None or name is None:
            continue
        data[name] = vtk_to_numpy(arr)
    return points, data


def boundary_mean(points: np.ndarray, values: np.ndarray, axis: int, side: str) -> np.ndarray:
    target = points[:, axis].max() if side == "max" else points[:, axis].min()
    mask = np.isclose(points[:, axis], target, atol=1e-10)
    return np.asarray(values[mask]).mean(axis=0)


def value_component(arr: np.ndarray, index: int, default: float = math.nan) -> float:
    flat = np.asarray(arr).reshape(-1)
    if flat.size == 0:
        return default
    return float(flat[index] if index < flat.size else flat[0])


def evaluate_beacon_metrics(vtu_path: Path) -> dict[str, float]:
    """Extract the final stresses, displacements, and density statistics."""
    points, data = load_vtu(vtu_path)
    sigma = data.get("sigma")
    swelling = data.get("swelling_stress")
    disp = data.get("displacement")
    dry_density = data.get("dry_density_solid")

    axial_sigma_kpa = radial_sigma_kpa = math.nan
    axial_swelling_kpa = radial_swelling_kpa = math.nan
    mean_total_stress_mpa = math.nan
    if sigma is not None:
        sigma_top = boundary_mean(points, sigma, axis=1, side="max")
        sigma_side = boundary_mean(points, sigma, axis=0, side="max")
        axial_sigma_kpa = abs(value_component(sigma_top, 1)) / 1e3
        radial_sigma_kpa = abs(value_component(sigma_side, 0)) / 1e3

        sig = np.asarray(sigma)
        if sig.ndim == 2 and sig.shape[1] >= 3:
            p_mean = (-(sig[:, 0] + sig[:, 1] + sig[:, 2]) / 3.0).mean()
            mean_total_stress_mpa = float(p_mean / 1e6)

    if swelling is not None:
        swelling_top = boundary_mean(points, swelling, axis=1, side="max")
        swelling_side = boundary_mean(points, swelling, axis=0, side="max")
        axial_swelling_kpa = abs(value_component(swelling_top, 1)) / 1e3
        radial_swelling_kpa = abs(value_component(swelling_side, 0)) / 1e3

    axial_disp_mm = radial_disp_mm = math.nan
    if disp is not None:
        disp_top = boundary_mean(points, disp, axis=1, side="max")
        disp_side = boundary_mean(points, disp, axis=0, side="max")
        axial_disp_mm = value_component(disp_top, 1) * 1e3
        radial_disp_mm = value_component(disp_side, 0) * 1e3

    return {
        "axial_sigma_kpa": axial_sigma_kpa,
        "radial_sigma_kpa": radial_sigma_kpa,
        "axial_swelling_stress_kpa": axial_swelling_kpa,
        "radial_swelling_stress_kpa": radial_swelling_kpa,
        "axial_displacement_mm": axial_disp_mm,
        "radial_displacement_mm": radial_disp_mm,
        "dry_density_mean_kg_m3": float(np.asarray(dry_density).mean()) if dry_density is not None else math.nan,
        "final_total_stress_mpa": mean_total_stress_mpa,
        "final_swelling_pressure_mpa": axial_sigma_kpa / 1e3 if not math.isnan(axial_sigma_kpa) else math.nan,
    }


def nan_beacon_metrics() -> dict[str, float]:
    """Return placeholder metrics for a failed BEACON run."""
    return {
        "axial_sigma_kpa": math.nan,
        "radial_sigma_kpa": math.nan,
        "axial_swelling_stress_kpa": math.nan,
        "radial_swelling_stress_kpa": math.nan,
        "axial_displacement_mm": math.nan,
        "radial_displacement_mm": math.nan,
        "dry_density_mean_kg_m3": math.nan,
        "final_total_stress_mpa": math.nan,
        "final_swelling_pressure_mpa": math.nan,
    }


def anchors_rows_from_calibrations() -> list[dict[str, object]]:
    """Convert the ANCHORS calibration CSVs into the common comparison schema."""
    mfront_rows = parse_csv_rows(MFRONT_CALIBRATION_CSV)
    native_rows = parse_csv_rows(NATIVE_CALIBRATION_CSV)
    native_by_dd = {float(r["dry_density_kg_m3"]): r for r in native_rows}

    out: list[dict[str, object]] = []
    for mr in mfront_rows:
        dd = float(mr["dry_density_kg_m3"])
        nr = native_by_dd.get(dd)
        if nr is None:
            continue

        mfront_multiplier = float(mr["vdw_multiplier"])
        native_multiplier = float(nr["vdw_multiplier"])
        mfront_ps = float(mr["mfront_calibrated_MPa"])
        native_ps = float(nr["native_calibrated_MPa"])

        out.append(
            {
                "benchmark_group": "ANCHORS_MS33_ModelI",
                "case_id": f"dd{int(dd)}",
                "implementation": "mfront",
                "dry_density_kg_m3": dd,
                "multiplier_used": mfront_multiplier,
                "hamaker_constant_J": float(mr["hamaker_effective_J"]),
                "final_total_stress_mpa": mfront_ps,
                "final_swelling_pressure_mpa": mfront_ps,
                "axial_sigma_kpa": mfront_ps * 1e3,
                "radial_sigma_kpa": math.nan,
                "axial_swelling_stress_kpa": math.nan,
                "radial_swelling_stress_kpa": math.nan,
                "axial_displacement_mm": math.nan,
                "radial_displacement_mm": math.nan,
                "dry_density_mean_kg_m3": dd,
                "run_status": "success",
                "source": MFRONT_CALIBRATION_CSV.name,
            }
        )
        out.append(
            {
                "benchmark_group": "ANCHORS_MS33_ModelI",
                "case_id": f"dd{int(dd)}",
                "implementation": "native",
                "dry_density_kg_m3": dd,
                "multiplier_used": native_multiplier,
                "hamaker_constant_J": float(nr["hamaker_effective_J"]),
                "final_total_stress_mpa": native_ps,
                "final_swelling_pressure_mpa": native_ps,
                "axial_sigma_kpa": native_ps * 1e3,
                "radial_sigma_kpa": math.nan,
                "axial_swelling_stress_kpa": math.nan,
                "radial_swelling_stress_kpa": math.nan,
                "axial_displacement_mm": math.nan,
                "radial_displacement_mm": math.nan,
                "dry_density_mean_kg_m3": dd,
                "run_status": "success",
                "source": NATIVE_CALIBRATION_CSV.name,
            }
        )
    return out


def run_beacon_rows(
    native_ogs: Path,
    mfront_ogs: Path,
    native_curve: tuple[np.ndarray, np.ndarray],
    mfront_curve: tuple[np.ndarray, np.ndarray],
) -> list[dict[str, object]]:
    """Execute BEACON cases with calibrated multipliers and collect metrics."""
    dd_native, mult_native = native_curve
    dd_mfront, mult_mfront = mfront_curve

    rows: list[dict[str, object]] = []
    tmpdir = Path(tempfile.mkdtemp(prefix="calibrated-beacon-runs-"))
    try:
        for case in BEACON_CASES:
            dd_native_case = extract_dd_from_project(case.native_project)
            dd_mfront_case = extract_dd_from_project(case.mfront_project)

            native_multiplier = interpolate_multiplier(dd_native_case, dd_native, mult_native)
            mfront_multiplier = interpolate_multiplier(dd_mfront_case, dd_mfront, mult_mfront)
            native_hamaker = HAMAKER_REFERENCE_J * native_multiplier
            mfront_hamaker = HAMAKER_REFERENCE_J * mfront_multiplier

            native_prefix = f"beacon_{case.case_id}_native_ddcal"
            mfront_prefix = f"beacon_{case.case_id}_mfront_ddcal"
            native_project = tmpdir / f"{native_prefix}.prj"
            mfront_project = tmpdir / f"{mfront_prefix}.prj"
            native_out = tmpdir / f"out_native_{case.case_id}"
            mfront_out = tmpdir / f"out_mfront_{case.case_id}"

            write_project_copy(
                source=case.native_project,
                target=native_project,
                prefix=native_prefix,
                hamaker_j=native_hamaker,
                implementation="native",
            )
            write_project_copy(
                source=case.mfront_project,
                target=mfront_project,
                prefix=mfront_prefix,
                hamaker_j=mfront_hamaker,
                implementation="mfront",
            )

            native_status = "success"
            native_error = ""
            mfront_status = "success"
            mfront_error = ""

            try:
                run_ogs(native_ogs, native_project, native_out)
                native_metrics = evaluate_beacon_metrics(latest_vtu(native_out, native_prefix))
            except Exception as exc:
                native_status = "failed"
                native_error = summarize_error(exc)
                native_metrics = nan_beacon_metrics()

            try:
                run_ogs(mfront_ogs, mfront_project, mfront_out)
                mfront_metrics = evaluate_beacon_metrics(latest_vtu(mfront_out, mfront_prefix))
            except Exception as exc:
                mfront_status = "failed"
                mfront_error = summarize_error(exc)
                mfront_metrics = nan_beacon_metrics()

            rows.append(
                {
                    "benchmark_group": "BEACON_report",
                    "case_id": case.case_id,
                    "implementation": "native",
                    "dry_density_kg_m3": dd_native_case,
                    "multiplier_used": native_multiplier,
                    "hamaker_constant_J": native_hamaker,
                    "run_status": native_status,
                    "error_message": native_error,
                    "source": case.native_project.name,
                    **native_metrics,
                }
            )
            rows.append(
                {
                    "benchmark_group": "BEACON_report",
                    "case_id": case.case_id,
                    "implementation": "mfront",
                    "dry_density_kg_m3": dd_mfront_case,
                    "multiplier_used": mfront_multiplier,
                    "hamaker_constant_J": mfront_hamaker,
                    "run_status": mfront_status,
                    "error_message": mfront_error,
                    "source": case.mfront_project.name,
                    **mfront_metrics,
                }
            )
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)
    return rows


def pairwise_deltas(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    """Compute native-minus-MFront deltas for each benchmark case."""
    by_key: dict[tuple[str, str], dict[str, dict[str, object]]] = {}
    for row in rows:
        key = (str(row["benchmark_group"]), str(row["case_id"]))
        impl = str(row["implementation"])
        by_key.setdefault(key, {})[impl] = row

    out: list[dict[str, object]] = []
    for (group, case_id), impls in sorted(by_key.items()):
        if "native" not in impls or "mfront" not in impls:
            continue
        n = impls["native"]
        m = impls["mfront"]

        def f(row: dict[str, object], key: str) -> float:
            value = row.get(key, math.nan)
            try:
                return float(value)
            except Exception:
                return math.nan

        out.append(
            {
                "benchmark_group": group,
                "case_id": case_id,
                "dry_density_kg_m3": f(n, "dry_density_kg_m3"),
                "native_multiplier": f(n, "multiplier_used"),
                "mfront_multiplier": f(m, "multiplier_used"),
                "native_to_mfront_multiplier_ratio": (
                    f(n, "multiplier_used") / f(m, "multiplier_used")
                    if abs(f(m, "multiplier_used")) > 0.0
                    else math.nan
                ),
                "native_final_swelling_pressure_mpa": f(n, "final_swelling_pressure_mpa"),
                "mfront_final_swelling_pressure_mpa": f(m, "final_swelling_pressure_mpa"),
                "delta_swelling_pressure_mpa_native_minus_mfront": (
                    f(n, "final_swelling_pressure_mpa") - f(m, "final_swelling_pressure_mpa")
                ),
                "native_axial_sigma_kpa": f(n, "axial_sigma_kpa"),
                "mfront_axial_sigma_kpa": f(m, "axial_sigma_kpa"),
                "delta_axial_sigma_kpa_native_minus_mfront": (
                    f(n, "axial_sigma_kpa") - f(m, "axial_sigma_kpa")
                ),
                "native_radial_sigma_kpa": f(n, "radial_sigma_kpa"),
                "mfront_radial_sigma_kpa": f(m, "radial_sigma_kpa"),
                "delta_radial_sigma_kpa_native_minus_mfront": (
                    f(n, "radial_sigma_kpa") - f(m, "radial_sigma_kpa")
                ),
                "native_axial_displacement_mm": f(n, "axial_displacement_mm"),
                "mfront_axial_displacement_mm": f(m, "axial_displacement_mm"),
                "delta_axial_displacement_mm_native_minus_mfront": (
                    f(n, "axial_displacement_mm") - f(m, "axial_displacement_mm")
                ),
                "native_radial_displacement_mm": f(n, "radial_displacement_mm"),
                "mfront_radial_displacement_mm": f(m, "radial_displacement_mm"),
                "delta_radial_displacement_mm_native_minus_mfront": (
                    f(n, "radial_displacement_mm") - f(m, "radial_displacement_mm")
                ),
            }
        )
    return out


def write_rows_csv(rows: list[dict[str, object]], path: Path) -> None:
    """Write a list of heterogeneous row dictionaries as a CSV table."""
    keys = sorted({k for row in rows for k in row.keys()})
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def git_hash(repo: Path) -> str:
    """Return a short Git hash for provenance tracking."""
    try:
        return (
            subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], cwd=repo, text=True)
            .strip()
        )
    except Exception:
        return ""


def main() -> None:
    """Run the full comparison workflow and write CSV/JSON artifacts."""
    args = parse_args()

    if not MFRONT_CALIBRATION_CSV.exists() or not NATIVE_CALIBRATION_CSV.exists():
        raise FileNotFoundError(
            "Missing calibration CSVs.\n"
            f"Expected: {MFRONT_CALIBRATION_CSV}\n"
            f"Expected: {NATIVE_CALIBRATION_CSV}"
        )

    rows = anchors_rows_from_calibrations()
    native_curve = load_multiplier_curve(NATIVE_CALIBRATION_CSV)
    mfront_curve = load_multiplier_curve(MFRONT_CALIBRATION_CSV)

    if not args.skip_beacon_runs:
        beacon_rows = run_beacon_rows(
            native_ogs=args.native_ogs.resolve(),
            mfront_ogs=args.mfront_ogs.resolve(),
            native_curve=native_curve,
            mfront_curve=mfront_curve,
        )
        rows.extend(beacon_rows)

    rows.sort(key=lambda r: (str(r["benchmark_group"]), str(r["case_id"]), str(r["implementation"])))
    deltas = pairwise_deltas(rows)
    deltas.sort(key=lambda r: (str(r["benchmark_group"]), str(r["case_id"])))

    write_rows_csv(rows, args.out_csv)
    write_rows_csv(deltas, args.out_deltas_csv)

    summary = {
        "reference_hamaker_J": HAMAKER_REFERENCE_J,
        "row_count": len(rows),
        "pairwise_delta_count": len(deltas),
        "sources": {
            "anchors_mfront_calibration_csv": str(MFRONT_CALIBRATION_CSV),
            "anchors_native_calibration_csv": str(NATIVE_CALIBRATION_CSV),
            "beacon_cases": [case.case_id for case in BEACON_CASES],
        },
        "git_hashes": {
            "ogs_repo_hash": git_hash(ROOT.parents[2]),
            "materialmodels_repo_hash": git_hash(Path("/Users/vinaykumar/Documents/GitHub/materialmodels")),
        },
        "notes": [
            (
                "Dry-density calibrated multipliers are implementation-specific: "
                "native uses villar_dense_dd_native_notebook_calibration.csv and "
                "MFront uses villar_dense_dd_calibration.csv."
            ),
            (
                "BEACON dry density is inferred from rho_s * (1 - phi0) in each "
                "project file; multiplier is linearly interpolated on DD and "
                "clamped to curve bounds."
            ),
            (
                "For ANCHORS, swelling pressure values are taken from calibrated DD runs; "
                "for BEACON, values come from fresh reruns with calibrated Hamaker values."
            ),
            (
                "Biot coefficient is enforced as 1.0 in generated BEACON run copies "
                "for both native and MFront implementations."
            ),
        ],
        "rows": rows,
        "pairwise_deltas": deltas,
    }
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    print(f"Wrote rows CSV: {args.out_csv}")
    print(f"Wrote summary JSON: {args.out_json}")
    print(f"Wrote pairwise deltas CSV: {args.out_deltas_csv}")


if __name__ == "__main__":
    main()
