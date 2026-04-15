#!/usr/bin/env python3
"""DD=1700 comparison: linear swelling (MCC) vs native DSM vs MFront DSM.

This script performs one focused comparison requested for the transition note:
1) Calibrate native and MFront DSM vdW multipliers at DD=1700 kg/m^3 using
   the existing Villar-curve calibration drivers.
2) Run three one-element confined hydration variants with MCC carrier:
   - linear SaturationDependentSwelling with max swelling pressure = 11e6 Pa,
   - native DSM + MCC with calibrated multiplier,
   - MFront DSM + MCC with calibrated multiplier.
3) Write CSV/JSON/PNG artifacts for documentation.
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import math
import os
import re
import shutil
import subprocess
import sys
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


ROOT = Path(__file__).resolve().parent
ANCHORS_DIR = ROOT / "ANCHORS_MS33_ModelI"
OUTPUT_ROOT = ROOT / "_outputs" / "linear_mcc_dd1700_dsm_comparison"


@dataclass(frozen=True)
class RunSpec:
    implementation: str
    project_path: Path
    ogs_bin: Path
    output_prefix: str


def load_module(module_path: Path, module_name: str):
    spec = importlib.util.spec_from_file_location(module_name, module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load module: {module_path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run DD=1700 linear(MCC, p_sw,max=11e6) vs native DSM vs MFront DSM "
            "comparison with calibrated multipliers."
        )
    )
    parser.add_argument(
        "--native-ogs",
        type=Path,
        default=Path("/Users/vinaykumar/git/build/release-native-transition2/bin/ogs"),
    )
    parser.add_argument(
        "--mfront-ogs",
        type=Path,
        default=Path("/Users/vinaykumar/git/build/release-mfront-tpm/bin/ogs"),
    )
    parser.add_argument("--dry-density", type=float, default=1700.0)
    parser.add_argument(
        "--max-swelling-pressure-pa",
        type=float,
        default=11e6,
        help="Linear SaturationDependentSwelling isotropic max pressure (Pa).",
    )
    parser.add_argument(
        "--rel-tol",
        type=float,
        default=0.02,
        help="Relative tolerance used in the two calibration routines.",
    )
    parser.add_argument(
        "--out-root",
        type=Path,
        default=OUTPUT_ROOT,
        help="Directory for CSV/JSON/PNG outputs.",
    )
    parser.add_argument(
        "--keep-runtime",
        action="store_true",
        help="Keep temporary project files in ANCHORS folder.",
    )
    return parser.parse_args()


def git_short_hash(repo: Path) -> str:
    try:
        return (
            subprocess.check_output(
                ["git", "rev-parse", "--short", "HEAD"], cwd=repo, text=True
            )
            .strip()
        )
    except Exception:
        return ""


def set_output_prefix(root: ET.Element, prefix: str) -> None:
    node = root.find("./time_loop/output/prefix")
    if node is None:
        raise RuntimeError("Output prefix node not found in project.")
    node.text = prefix


def limit_outputs_to_linear_core(root: ET.Element) -> None:
    allowed = {
        "pressure",
        "displacement",
        "sigma",
        "swelling_stress",
        "saturation",
        "porosity",
        "dry_density_solid",
    }
    sec = root.find("./processes/process/secondary_variables")
    if sec is not None:
        for child in list(sec):
            name = child.get("name", "").strip()
            if name not in allowed:
                sec.remove(child)

    out_vars = root.find("./time_loop/output/variables")
    if out_vars is not None:
        for child in list(out_vars):
            name = (child.text or "").strip()
            if name not in allowed:
                out_vars.remove(child)


def configure_linear_mcc_project(
    source_project: Path, target_project: Path, prefix: str, max_swelling_pressure_pa: float
) -> None:
    root = ET.parse(source_project).getroot()
    process = root.find("./processes/process")
    if process is None:
        raise RuntimeError("Missing process node in source project.")

    pe = process.find("potential_exchange")
    if pe is not None:
        process.remove(pe)
    mp = process.find("micro_porosity")
    if mp is not None:
        process.remove(mp)

    sw = root.find(
        "./media/medium/phases/phase[type='Solid']/properties/property[name='swelling_stress_rate']"
    )
    if sw is None:
        raise RuntimeError("Solid swelling_stress_rate property not found.")
    sw_type = sw.find("type")
    if sw_type is None:
        sw_type = ET.SubElement(sw, "type")
    sw_type.text = "SaturationDependentSwelling"
    pressures = sw.find("swelling_pressures")
    if pressures is None:
        pressures = ET.SubElement(sw, "swelling_pressures")
    pressures.text = (
        f"{max_swelling_pressure_pa:.16g} "
        f"{max_swelling_pressure_pa:.16g} "
        f"{max_swelling_pressure_pa:.16g}"
    )

    exponents = sw.find("exponents")
    if exponents is None:
        exponents = ET.SubElement(sw, "exponents")
    exponents.text = "1 1 1"
    smin = sw.find("lower_saturation_limit")
    if smin is None:
        smin = ET.SubElement(sw, "lower_saturation_limit")
    smin.text = "0"
    smax = sw.find("upper_saturation_limit")
    if smax is None:
        smax = ET.SubElement(sw, "upper_saturation_limit")
    smax.text = "1"

    set_output_prefix(root, prefix)
    limit_outputs_to_linear_core(root)
    ET.ElementTree(root).write(target_project, encoding="ISO-8859-1", xml_declaration=True)


def parse_time_from_vtu(path: Path) -> float:
    stem = path.stem
    match = re.search(r"_t_([-0-9eE+.]+)$", stem)
    if not match:
        raise RuntimeError(f"Cannot parse time from output filename {path.name}")
    return float(match.group(1))


def latest_vtu(output_dir: Path, prefix: str) -> Path:
    files = sorted(output_dir.glob(f"{prefix}_ts_*_t_*.vtu"), key=parse_time_from_vtu)
    if not files:
        raise FileNotFoundError(f"No VTU outputs for prefix {prefix} in {output_dir}")
    return files[-1]


def load_vtu(path: Path) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(path))
    reader.Update()
    grid = reader.GetOutput()
    points = np.array([grid.GetPoint(i) for i in range(grid.GetNumberOfPoints())])
    arrays: dict[str, np.ndarray] = {}

    point_data = grid.GetPointData()
    for i in range(point_data.GetNumberOfArrays()):
        arr = point_data.GetArray(i)
        if arr is None:
            continue
        name = point_data.GetArrayName(i)
        if name is None:
            continue
        arrays[name] = vtk_to_numpy(arr)

    cell_data = grid.GetCellData()
    for i in range(cell_data.GetNumberOfArrays()):
        arr = cell_data.GetArray(i)
        if arr is None:
            continue
        name = cell_data.GetArrayName(i)
        if name is None or name in arrays:
            continue
        arrays[name] = vtk_to_numpy(arr)
    return points, arrays


def boundary_mask(points: np.ndarray, axis: int, side: str) -> np.ndarray:
    val = points[:, axis].max() if side == "max" else points[:, axis].min()
    return np.isclose(points[:, axis], val, atol=1e-10)


def mean_component(values: np.ndarray, mask: np.ndarray, idx: int) -> float:
    arr = np.asarray(values)
    if arr.ndim == 1:
        return float(arr[mask].mean())
    return float(arr[mask, idx].mean())


def mean_isotropic_pressure_mpa(tensor: np.ndarray) -> float:
    arr = np.asarray(tensor)
    p_pa = float((-arr[:, 0] - arr[:, 1] - arr[:, 2]).mean() / 3.0)
    return p_pa / 1e6


def evaluate_metrics(vtu_path: Path) -> dict[str, float]:
    points, arrays = load_vtu(vtu_path)
    sigma = arrays.get("sigma")
    swelling = arrays.get("swelling_stress")
    displacement = arrays.get("displacement")
    saturation = arrays.get("saturation")
    porosity = arrays.get("porosity")
    dry_density = arrays.get("dry_density_solid")

    top = boundary_mask(points, axis=1, side="max")
    side = boundary_mask(points, axis=0, side="max")

    metrics: dict[str, float] = {
        "final_total_stress_mpa": mean_isotropic_pressure_mpa(sigma)
        if sigma is not None
        else math.nan,
        "final_swelling_pressure_mpa": mean_isotropic_pressure_mpa(swelling)
        if swelling is not None
        else math.nan,
        "axial_sigma_kpa": abs(mean_component(sigma, top, 1)) / 1000.0
        if sigma is not None
        else math.nan,
        "radial_sigma_kpa": abs(mean_component(sigma, side, 0)) / 1000.0
        if sigma is not None
        else math.nan,
        "axial_swelling_stress_kpa": abs(mean_component(swelling, top, 1)) / 1000.0
        if swelling is not None
        else math.nan,
        "radial_swelling_stress_kpa": abs(mean_component(swelling, side, 0)) / 1000.0
        if swelling is not None
        else math.nan,
        "axial_displacement_mm": mean_component(displacement, top, 1) * 1000.0
        if displacement is not None
        else math.nan,
        "radial_displacement_mm": mean_component(displacement, side, 0) * 1000.0
        if displacement is not None
        else math.nan,
        "saturation_mean": float(np.asarray(saturation).mean())
        if saturation is not None
        else math.nan,
        "porosity_mean": float(np.asarray(porosity).mean())
        if porosity is not None
        else math.nan,
        "dry_density_mean_kg_m3": float(np.asarray(dry_density).mean())
        if dry_density is not None
        else math.nan,
    }
    return metrics


def collect_history(output_dir: Path, prefix: str) -> list[dict[str, float]]:
    files = sorted(output_dir.glob(f"{prefix}_ts_*_t_*.vtu"), key=parse_time_from_vtu)
    history: list[dict[str, float]] = []
    for path in files:
        t = parse_time_from_vtu(path)
        metrics = evaluate_metrics(path)
        history.append(
            {
                "time_s": t,
                "time_d": t / 86400.0,
                "total_stress_mpa": metrics["final_total_stress_mpa"],
                "swelling_pressure_mpa": metrics["final_swelling_pressure_mpa"],
            }
        )
    return history


def run_ogs(ogs_bin: Path, project_path: Path, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    result = subprocess.run(
        [str(ogs_bin), "-o", str(output_dir), str(project_path)],
        cwd=project_path.parent,
        text=True,
        capture_output=True,
        check=False,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"OGS failed for {project_path.name} with {ogs_bin}\n"
            f"stdout:\n{result.stdout}\n"
            f"stderr:\n{result.stderr}"
        )


def write_rows_csv(rows: list[dict[str, object]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    keys: list[str] = []
    seen: set[str] = set()
    for row in rows:
        for key in row:
            if key not in seen:
                seen.add(key)
                keys.append(key)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def pairwise_deltas(rows: list[dict[str, object]]) -> list[dict[str, float]]:
    by_impl = {str(r["implementation"]): r for r in rows}
    lin = by_impl["linear_mcc"]
    nat = by_impl["native_dsm_mcc"]
    mfr = by_impl["mfront_dsm_mcc"]

    def val(r: dict[str, object], key: str) -> float:
        return float(r.get(key, math.nan))

    p_lin = val(lin, "final_swelling_pressure_mpa")
    p_nat = val(nat, "final_swelling_pressure_mpa")
    p_mfr = val(mfr, "final_swelling_pressure_mpa")
    s_lin = val(lin, "final_total_stress_mpa")
    s_nat = val(nat, "final_total_stress_mpa")
    s_mfr = val(mfr, "final_total_stress_mpa")

    return [
        {
            "quantity": "swelling_pressure_mpa",
            "linear_mcc": p_lin,
            "native_dsm_mcc": p_nat,
            "mfront_dsm_mcc": p_mfr,
            "delta_native_minus_linear": p_nat - p_lin,
            "delta_mfront_minus_linear": p_mfr - p_lin,
            "delta_native_minus_mfront": p_nat - p_mfr,
            "rel_native_vs_linear": abs(p_nat - p_lin) / max(abs(p_lin), 1e-12),
            "rel_mfront_vs_linear": abs(p_mfr - p_lin) / max(abs(p_lin), 1e-12),
            "rel_native_vs_mfront": abs(p_nat - p_mfr) / max(abs(p_mfr), 1e-12),
        },
        {
            "quantity": "total_stress_mpa",
            "linear_mcc": s_lin,
            "native_dsm_mcc": s_nat,
            "mfront_dsm_mcc": s_mfr,
            "delta_native_minus_linear": s_nat - s_lin,
            "delta_mfront_minus_linear": s_mfr - s_lin,
            "delta_native_minus_mfront": s_nat - s_mfr,
            "rel_native_vs_linear": abs(s_nat - s_lin) / max(abs(s_lin), 1e-12),
            "rel_mfront_vs_linear": abs(s_mfr - s_lin) / max(abs(s_lin), 1e-12),
            "rel_native_vs_mfront": abs(s_nat - s_mfr) / max(abs(s_mfr), 1e-12),
        },
    ]


def plot_results(rows: list[dict[str, object]], histories: dict[str, list[dict[str, float]]], out_png: Path) -> None:
    impl_order = ["linear_mcc", "native_dsm_mcc", "mfront_dsm_mcc"]
    labels = ["Linear+MCC (11 MPa)", "Native DSM+MCC", "MFront DSM+MCC"]
    p = [float(next(r for r in rows if r["implementation"] == impl)["final_swelling_pressure_mpa"]) for impl in impl_order]

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))

    x = np.arange(len(labels))
    axes[0].bar(x, p, color=["#4e79a7", "#f28e2b", "#59a14f"])
    axes[0].set_xticks(x, labels, rotation=12, ha="right")
    axes[0].set_ylabel("Final swelling pressure (MPa)")
    axes[0].set_title("DD=1700 final comparison")
    axes[0].grid(True, axis="y", ls=":", alpha=0.4)

    for impl, color, label in zip(
        impl_order,
        ["#4e79a7", "#f28e2b", "#59a14f"],
        labels,
    ):
        hist = histories.get(impl, [])
        if not hist:
            continue
        t = np.array([r["time_d"] for r in hist], dtype=float)
        y = np.array([r["total_stress_mpa"] for r in hist], dtype=float)
        axes[1].plot(t, y, marker="o", ms=3, lw=1.2, color=color, label=label)
    axes[1].set_xlabel("Time (days)")
    axes[1].set_ylabel("Mean total stress (MPa)")
    axes[1].set_title("Stress evolution")
    axes[1].grid(True, ls=":", alpha=0.4)
    axes[1].legend(loc="best")

    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=220)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    args.out_root.mkdir(parents=True, exist_ok=True)

    mfront_mod = load_module(
        ANCHORS_DIR / "run_villar_dense_dd_calibration.py",
        "run_villar_dense_dd_calibration_dd1700_triplet",
    )
    native_mod = load_module(
        ANCHORS_DIR / "run_villar_dense_dd_native_dsm_micromacro_calibration.py",
        "run_villar_dense_dd_native_dsm_micromacro_calibration_dd1700_triplet",
    )

    case_m = mfront_mod.Case(float(args.dry_density))
    case_n = native_mod.Case(float(args.dry_density))
    target_mpa = float(case_m.villar_target_swelling_mpa)

    mfront_cal = mfront_mod.calibrate_multiplier_for_case(
        args.mfront_ogs, case_m, target_mpa, rel_tol=args.rel_tol
    )
    n_l0_native = native_mod.n_l0_from_micro_suction(case_n.phi0, native_mod.HAMAKER_LITERATURE)
    native_cal = native_mod.calibrate_multiplier_for_case(
        args.native_ogs, case_n, n_l0_native, target_mpa, rel_tol=args.rel_tol
    )

    runtime_prefix = f"dd{int(args.dry_density)}_linear_mcc_triplet"
    native_dsm_project = ANCHORS_DIR / f"{runtime_prefix}_native_dsm.prj"
    mfront_dsm_project = ANCHORS_DIR / f"{runtime_prefix}_mfront_dsm.prj"
    linear_mcc_project = ANCHORS_DIR / f"{runtime_prefix}_linear_mcc.prj"

    native_prefix = f"{runtime_prefix}_native_dsm"
    mfront_prefix = f"{runtime_prefix}_mfront_dsm"
    linear_prefix = f"{runtime_prefix}_linear_mcc"

    native_mod.write_native_dsm_micromacro_project(
        case_n,
        native_dsm_project,
        float(native_cal["multiplier"]),
        float(n_l0_native),
    )
    mfront_mod.write_mfront_project(
        case_m,
        float(mfront_cal["multiplier"]),
        float(mfront_cal.get("n_l0", mfront_mod.n_l0_from_micro_suction(case_m.phi0, mfront_mod.HAMAKER_LITERATURE))),
        mfront_dsm_project,
    )
    configure_linear_mcc_project(
        native_dsm_project,
        linear_mcc_project,
        linear_prefix,
        float(args.max_swelling_pressure_pa),
    )

    # Ensure output prefixes are unique.
    native_root = ET.parse(native_dsm_project).getroot()
    set_output_prefix(native_root, native_prefix)
    ET.ElementTree(native_root).write(
        native_dsm_project, encoding="ISO-8859-1", xml_declaration=True
    )

    mfront_root = ET.parse(mfront_dsm_project).getroot()
    set_output_prefix(mfront_root, mfront_prefix)
    ET.ElementTree(mfront_root).write(
        mfront_dsm_project, encoding="ISO-8859-1", xml_declaration=True
    )

    run_specs = [
        RunSpec("linear_mcc", linear_mcc_project, args.native_ogs, linear_prefix),
        RunSpec("native_dsm_mcc", native_dsm_project, args.native_ogs, native_prefix),
        RunSpec("mfront_dsm_mcc", mfront_dsm_project, args.mfront_ogs, mfront_prefix),
    ]

    rows: list[dict[str, object]] = []
    histories: dict[str, list[dict[str, float]]] = {}
    run_root = args.out_root / "runs"
    for spec in run_specs:
        out_dir = run_root / spec.implementation
        run_ogs(spec.ogs_bin, spec.project_path, out_dir)
        final_vtu = latest_vtu(out_dir, spec.output_prefix)
        metrics = evaluate_metrics(final_vtu)
        histories[spec.implementation] = collect_history(out_dir, spec.output_prefix)
        rows.append(
            {
                "implementation": spec.implementation,
                "dry_density_kg_m3": float(args.dry_density),
                "max_swelling_pressure_pa_linear": float(args.max_swelling_pressure_pa)
                if spec.implementation == "linear_mcc"
                else math.nan,
                "vdw_multiplier": math.nan
                if spec.implementation == "linear_mcc"
                else (
                    float(native_cal["multiplier"])
                    if spec.implementation == "native_dsm_mcc"
                    else float(mfront_cal["multiplier"])
                ),
                "hamaker_effective_j": math.nan
                if spec.implementation == "linear_mcc"
                else (
                    float(native_cal["hamaker_effective_J"])
                    if spec.implementation == "native_dsm_mcc"
                    else float(mfront_cal["hamaker_effective"])
                ),
                "calibration_target_mpa": target_mpa,
                "output_vtu": str(final_vtu.resolve()),
                **metrics,
            }
        )

    deltas = pairwise_deltas(rows)

    rows_csv = args.out_root / "dd1700_linear_mcc_dsm_runs_summary.csv"
    deltas_csv = args.out_root / "dd1700_linear_mcc_dsm_pairwise_deltas.csv"
    history_csv = args.out_root / "dd1700_linear_mcc_dsm_history.csv"
    figure_png = args.out_root / "dd1700_linear_mcc_dsm_comparison.png"
    summary_json = args.out_root / "dd1700_linear_mcc_dsm_summary.json"

    write_rows_csv(rows, rows_csv)
    write_rows_csv(deltas, deltas_csv)
    flat_history = []
    for impl, hist in histories.items():
        for row in hist:
            flat_history.append({"implementation": impl, **row})
    write_rows_csv(flat_history, history_csv)
    plot_results(rows, histories, figure_png)

    summary = {
        "mode": "dd1700_linear_mcc_vs_dsm_triplet",
        "dry_density_kg_m3": float(args.dry_density),
        "max_swelling_pressure_pa_linear": float(args.max_swelling_pressure_pa),
        "calibration_target_mpa_villar": target_mpa,
        "native_calibration": {
            "multiplier": float(native_cal["multiplier"]),
            "pressure_mpa": float(native_cal["pressure_mpa"]),
            "hamaker_effective_j": float(native_cal["hamaker_effective_J"]),
            "n_l0": float(native_cal["n_l0"]),
        },
        "mfront_calibration": {
            "multiplier": float(mfront_cal["multiplier"]),
            "pressure_mpa": float(mfront_cal["pressure_mpa"]),
            "hamaker_effective_j": float(mfront_cal["hamaker_effective"]),
            "n_l0": float(mfront_cal["n_l0"]),
        },
        "runs": rows,
        "pairwise_deltas": deltas,
        "git_hashes": {
            "mfront_repo_hash": git_short_hash(ROOT.parent.parent.parent),
            "native_repo_hash": git_short_hash(Path("/Users/vinaykumar/git/ogs-native-dsm-transition")),
        },
        "artifacts": {
            "rows_csv": str(rows_csv.resolve()),
            "deltas_csv": str(deltas_csv.resolve()),
            "history_csv": str(history_csv.resolve()),
            "figure_png": str(figure_png.resolve()),
        },
    }
    summary_json.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    print(f"Wrote rows CSV: {rows_csv}")
    print(f"Wrote pairwise deltas CSV: {deltas_csv}")
    print(f"Wrote history CSV: {history_csv}")
    print(f"Wrote plot: {figure_png}")
    print(f"Wrote summary JSON: {summary_json}")

    if not args.keep_runtime:
        for project in (native_dsm_project, mfront_dsm_project, linear_mcc_project):
            project.unlink(missing_ok=True)
        for p in ANCHORS_DIR.glob(f"{runtime_prefix}_*.pvd"):
            p.unlink(missing_ok=True)
        for p in ANCHORS_DIR.glob(f"{runtime_prefix}_*_ts_*_t_*.vtu"):
            p.unlink(missing_ok=True)
        shutil.rmtree(run_root, ignore_errors=True)


if __name__ == "__main__":
    main()
