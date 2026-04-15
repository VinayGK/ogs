#!/usr/bin/env python3
"""Calibrate DSM to OGS linear swelling law, then compare benchmarks.

Workflow:
1) Build ANCHORS dry-density calibration curves for DSM native and DSM MFront
   against the OGS linear SaturationDependentSwelling reference response.
2) Rerun BEACON + EPFL with three variants:
   - linear swelling-law baseline (native project, potential_exchange disabled),
   - DSM native (potential_exchange enabled, linear swelling disabled),
   - DSM MFront (bridge behaviour with calibrated Hamaker).
3) Merge ANCHORS, BEACON, EPFL into one comparison table and summary.
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path

import numpy as np

os.environ.setdefault("MPLCONFIGDIR", tempfile.mkdtemp(prefix="mplconfig-"))
os.environ.setdefault("XDG_CACHE_HOME", tempfile.mkdtemp(prefix="xdg-cache-"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parent
ANCHORS_DIR = ROOT / "ANCHORS_MS33_ModelI"
EPFL_DIR = ROOT / "BEACON_WP3_BGR_EPFL"
OUTPUT_ROOT = ROOT / "_outputs" / "linear_law_calibrated_dsm_comparison"
HAMAKER_REFERENCE_J = 5.1e-21


@dataclass(frozen=True)
class BenchmarkCase:
    benchmark_group: str
    case_id: str
    native_project: Path
    mfront_project: Path


BEACON_AND_EPFL_CASES = [
    BenchmarkCase(
        benchmark_group="BEACON_report",
        case_id="1a01",
        native_project=ROOT / "beacon_1a01_inflow_unstructured_batch.prj",
        mfront_project=ROOT / "beacon_1a01_dsm_micromacro_mcc_inflow_unstructured_batch.prj",
    ),
    BenchmarkCase(
        benchmark_group="BEACON_report",
        case_id="1b",
        native_project=ROOT / "beacon_1b_unstructured_batch.prj",
        mfront_project=ROOT / "beacon_1b_dsm_micromacro_mcc_unstructured_batch.prj",
    ),
    BenchmarkCase(
        benchmark_group="EPFL_BGR_WP3",
        case_id="p2_1_abprime",
        native_project=EPFL_DIR / "bgr_wp3_p2_1_abprime_native.prj",
        mfront_project=EPFL_DIR / "bgr_wp3_p2_1_abprime_mfront.prj",
    ),
]


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
            "Calibrate DSM native + MFront against the OGS linear swelling-law "
            "reference, then compare ANCHORS + BEACON + EPFL."
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
    parser.add_argument(
        "--linear-ref-ogs",
        type=Path,
        default=Path("/Users/vinaykumar/git/build/release-mfront-tpm/bin/ogs"),
        help="OGS binary used for linear-law ANCHORS reference runs.",
    )
    parser.add_argument("--dd-min", type=float, default=1400.0)
    parser.add_argument("--dd-max", type=float, default=1800.0)
    parser.add_argument("--dd-step", type=float, default=25.0)
    parser.add_argument("--rel-tol", type=float, default=0.02)
    parser.add_argument(
        "--out-root",
        type=Path,
        default=OUTPUT_ROOT,
        help="Output directory for CSV/JSON/PNG artifacts.",
    )
    parser.add_argument(
        "--assert-thresholds",
        action="store_true",
        help="Fail if ANCHORS calibration relative errors exceed rel-tol.",
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


def write_rows_csv(rows: list[dict[str, object]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    keys: list[str] = []
    seen: set[str] = set()
    for row in rows:
        for k in row.keys():
            if k not in seen:
                seen.add(k)
                keys.append(k)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def interpolate_multiplier(dd: float, dd_curve: np.ndarray, mult_curve: np.ndarray) -> float:
    return float(np.interp(dd, dd_curve, mult_curve, left=mult_curve[0], right=mult_curve[-1]))


def remove_potential_exchange(root: ET.Element) -> None:
    process = root.find("./processes/process")
    if process is None:
        return
    pe = process.find("potential_exchange")
    if pe is not None:
        process.remove(pe)


def zero_swelling_stress_rate(root: ET.Element) -> None:
    for prop in root.findall(".//property"):
        name = (prop.findtext("name") or "").strip()
        if name != "swelling_stress_rate":
            continue
        s = prop.find("swelling_pressures")
        if s is not None:
            s.text = "0 0 0"


def find_existing_parameter(root: ET.Element, names: list[str], fallback: str) -> str:
    existing = {
        (node.findtext("name") or "").strip()
        for node in root.findall("./parameters/parameter")
    }
    for n in names:
        if n in existing:
            return n
    return fallback


def set_linear_elastic_constitutive(root: ET.Element) -> None:
    """Replace constitutive relation by native linear elastic isotropic law."""
    cr = root.find("./processes/process/constitutive_relation")
    if cr is None:
        return
    young_name = find_existing_parameter(root, ["YoungModulus", "E"], "YoungModulus")
    poisson_name = find_existing_parameter(root, ["PoissonRatio", "nu"], "PoissonRatio")
    for child in list(cr):
        cr.remove(child)
    ET.SubElement(cr, "type").text = "LinearElasticIsotropic"
    ET.SubElement(cr, "youngs_modulus").text = young_name
    ET.SubElement(cr, "poissons_ratio").text = poisson_name


def configure_linear_reference_outputs(root: ET.Element) -> None:
    """Limit outputs/secondary variables to fields available in linear reference."""
    sec = root.find("./processes/process/secondary_variables")
    if sec is not None:
        for child in list(sec):
            sec.remove(child)
        for name in [
            "sigma",
            "swelling_stress",
            "saturation",
            "porosity",
            "transport_porosity",
            "dry_density_solid",
            "epsilon",
            "velocity",
        ]:
            sv = ET.SubElement(sec, "secondary_variable")
            sv.set("name", name)

    out_vars = root.find("./time_loop/output/variables")
    if out_vars is not None:
        for child in list(out_vars):
            out_vars.remove(child)
        for name in [
            "pressure",
            "displacement",
            "saturation",
            "swelling_stress",
            "sigma",
            "porosity",
            "transport_porosity",
            "dry_density_solid",
        ]:
            v = ET.SubElement(out_vars, "variable")
            v.text = name


def run_project_variant(
    *,
    root: ET.Element,
    source_project: Path,
    prefix: str,
    project_path: Path,
    output_dir: Path,
    ogs_bin: Path,
    compare_mod,
    max_iter: int = 80,
    first_abstol: float = 5e-7,
) -> tuple[dict[str, float], str, str]:
    compare_mod.absolutize_mesh_and_geometry(root, source_project.parent)
    compare_mod.set_output_prefix(root, prefix)
    compare_mod.set_biot_coefficient(root, 1.0)
    compare_mod.set_nonlinear_max_iter(root, max_iter)
    compare_mod.relax_first_component_abstol(root, first_abstol)
    ET.ElementTree(root).write(project_path, encoding="ISO-8859-1", xml_declaration=True)
    try:
        compare_mod.run_ogs(ogs_bin, project_path, output_dir)
        metrics = compare_mod.evaluate_beacon_metrics(compare_mod.latest_vtu(output_dir, prefix))
        return metrics, "success", ""
    except Exception as exc:
        return compare_mod.nan_beacon_metrics(), "failed", compare_mod.summarize_error(exc)


def calibration_rows_and_curves(args: argparse.Namespace, mfront_mod, native_mod) -> tuple[
    list[dict[str, object]],
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
]:
    dd_values = np.arange(args.dd_min, args.dd_max + 0.5 * args.dd_step, args.dd_step)
    rows: list[dict[str, object]] = []
    for dd in dd_values:
        mcase = mfront_mod.Case(float(dd))
        ncase = native_mod.Case(float(dd))

        linear_target = mfront_mod.run_native_case(args.linear_ref_ogs, mcase)

        mfront_cal = mfront_mod.calibrate_multiplier_for_case(
            args.mfront_ogs, mcase, linear_target, rel_tol=args.rel_tol
        )

        n_l0_fixed = native_mod.n_l0_from_micro_suction(
            ncase.phi0, native_mod.HAMAKER_LITERATURE
        )
        native_cal = native_mod.calibrate_multiplier_for_case(
            args.native_ogs,
            ncase,
            n_l0_fixed,
            linear_target,
            rel_tol=args.rel_tol,
        )

        mfront_rel_err = abs(mfront_cal["pressure_mpa"] - linear_target) / max(
            abs(linear_target), 1e-12
        )
        native_rel_err = abs(native_cal["pressure_mpa"] - linear_target) / max(
            abs(linear_target), 1e-12
        )

        rows.append(
            {
                "dry_density_kg_m3": float(dd),
                "linear_target_MPa": float(linear_target),
                "mfront_calibrated_MPa": float(mfront_cal["pressure_mpa"]),
                "native_calibrated_MPa": float(native_cal["pressure_mpa"]),
                "mfront_vdw_multiplier": float(mfront_cal["multiplier"]),
                "native_vdw_multiplier": float(native_cal["multiplier"]),
                "mfront_hamaker_effective_J": float(mfront_cal["hamaker_effective"]),
                "native_hamaker_effective_J": float(native_cal["hamaker_effective_J"]),
                "mfront_relative_error": float(mfront_rel_err),
                "native_relative_error": float(native_rel_err),
            }
        )
        print(
            f"ANCHORS linear calibration dd={dd:.0f}: "
            f"linear={linear_target:.6f} MPa, "
            f"native={native_cal['pressure_mpa']:.6f} MPa (m={native_cal['multiplier']:.6e}), "
            f"mfront={mfront_cal['pressure_mpa']:.6f} MPa (m={mfront_cal['multiplier']:.6e})"
        )

    rows = sorted(rows, key=lambda r: float(r["dry_density_kg_m3"]))
    dd_curve = np.array([float(r["dry_density_kg_m3"]) for r in rows], dtype=float)
    native_curve = np.array([float(r["native_vdw_multiplier"]) for r in rows], dtype=float)
    mfront_curve = np.array([float(r["mfront_vdw_multiplier"]) for r in rows], dtype=float)
    return rows, dd_curve, native_curve, dd_curve, mfront_curve


def anchors_triplet_rows(calibration_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for r in calibration_rows:
        dd = float(r["dry_density_kg_m3"])
        linear = float(r["linear_target_MPa"])
        native = float(r["native_calibrated_MPa"])
        mfront = float(r["mfront_calibrated_MPa"])
        rows.extend(
            [
                {
                    "benchmark_group": "ANCHORS_MS33_ModelI",
                    "case_id": f"dd{int(dd)}",
                    "implementation": "linear",
                    "dry_density_kg_m3": dd,
                    "multiplier_used": math.nan,
                    "hamaker_constant_J": math.nan,
                    "final_total_stress_mpa": linear,
                    "final_swelling_pressure_mpa": linear,
                    "axial_sigma_kpa": linear * 1e3,
                    "radial_sigma_kpa": math.nan,
                    "axial_swelling_stress_kpa": math.nan,
                    "radial_swelling_stress_kpa": math.nan,
                    "axial_displacement_mm": math.nan,
                    "radial_displacement_mm": math.nan,
                    "dry_density_mean_kg_m3": dd,
                    "run_status": "success",
                    "error_message": "",
                    "source": "linear_reference_native_project",
                },
                {
                    "benchmark_group": "ANCHORS_MS33_ModelI",
                    "case_id": f"dd{int(dd)}",
                    "implementation": "native_dsm",
                    "dry_density_kg_m3": dd,
                    "multiplier_used": float(r["native_vdw_multiplier"]),
                    "hamaker_constant_J": float(r["native_hamaker_effective_J"]),
                    "final_total_stress_mpa": native,
                    "final_swelling_pressure_mpa": native,
                    "axial_sigma_kpa": native * 1e3,
                    "radial_sigma_kpa": math.nan,
                    "axial_swelling_stress_kpa": math.nan,
                    "radial_swelling_stress_kpa": math.nan,
                    "axial_displacement_mm": math.nan,
                    "radial_displacement_mm": math.nan,
                    "dry_density_mean_kg_m3": dd,
                    "run_status": "success",
                    "error_message": "",
                    "source": "anchors_linear_calibration_native",
                },
                {
                    "benchmark_group": "ANCHORS_MS33_ModelI",
                    "case_id": f"dd{int(dd)}",
                    "implementation": "mfront_dsm",
                    "dry_density_kg_m3": dd,
                    "multiplier_used": float(r["mfront_vdw_multiplier"]),
                    "hamaker_constant_J": float(r["mfront_hamaker_effective_J"]),
                    "final_total_stress_mpa": mfront,
                    "final_swelling_pressure_mpa": mfront,
                    "axial_sigma_kpa": mfront * 1e3,
                    "radial_sigma_kpa": math.nan,
                    "axial_swelling_stress_kpa": math.nan,
                    "radial_swelling_stress_kpa": math.nan,
                    "axial_displacement_mm": math.nan,
                    "radial_displacement_mm": math.nan,
                    "dry_density_mean_kg_m3": dd,
                    "run_status": "success",
                    "error_message": "",
                    "source": "anchors_linear_calibration_mfront",
                },
            ]
        )
    return rows


def mfront_defaults_for_epfl(root: ET.Element, parity_mod) -> dict[str, float]:
    return {
        "ReferenceDensitySolid": parity_mod.get_parameter_value(
            root, "ReferenceDensitySolid", default=2780.0
        ),
        "SpecificSurface": parity_mod.get_parameter_value(root, "SpecificSurface", default=523.0),
        "MicroSwellingStrainSlope": parity_mod.get_parameter_value(
            root, "MicroSwellingStrainSlope", default=0.0
        ),
        "MassExchangeCoefficient": parity_mod.get_parameter_value(
            root, "MassExchangeCoefficient", default=5e-15
        ),
        "n_l0": parity_mod.get_parameter_value(root, "n_l0", default=1e-12),
    }


def run_beacon_epfl_triplet_rows(
    *,
    args: argparse.Namespace,
    compare_mod,
    parity_mod,
    dd_curve_native: np.ndarray,
    mult_curve_native: np.ndarray,
    dd_curve_mfront: np.ndarray,
    mult_curve_mfront: np.ndarray,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    tmpdir = Path(tempfile.mkdtemp(prefix="linear-law-triplet-runs-"))
    try:
        for case in BEACON_AND_EPFL_CASES:
            dd_native = compare_mod.extract_dd_from_project(case.native_project)
            dd_mfront = compare_mod.extract_dd_from_project(case.mfront_project)
            native_multiplier = interpolate_multiplier(
                dd_native, dd_curve_native, mult_curve_native
            )
            mfront_multiplier = interpolate_multiplier(
                dd_mfront, dd_curve_mfront, mult_curve_mfront
            )
            hamaker_native = HAMAKER_REFERENCE_J * native_multiplier
            hamaker_mfront = HAMAKER_REFERENCE_J * mfront_multiplier

            # 1) Linear baseline on native project with potential_exchange disabled.
            linear_prefix = f"{case.benchmark_group}_{case.case_id}_linear".lower()
            linear_project = tmpdir / f"{linear_prefix}.prj"
            linear_out = tmpdir / f"out_{linear_prefix}"
            linear_root = ET.parse(case.native_project).getroot()
            remove_potential_exchange(linear_root)
            set_linear_elastic_constitutive(linear_root)
            configure_linear_reference_outputs(linear_root)
            linear_metrics, linear_status, linear_error = run_project_variant(
                root=linear_root,
                source_project=case.native_project,
                prefix=linear_prefix,
                project_path=linear_project,
                output_dir=linear_out,
                ogs_bin=args.native_ogs,
                compare_mod=compare_mod,
                max_iter=200,
                first_abstol=1e-5,
            )

            # 2) Native DSM with linear swelling disabled.
            native_prefix = f"{case.benchmark_group}_{case.case_id}_native_dsm".lower()
            native_project = tmpdir / f"{native_prefix}.prj"
            native_out = tmpdir / f"out_{native_prefix}"
            native_root = ET.parse(case.native_project).getroot()
            pe = native_root.find("./processes/process/potential_exchange")
            if pe is None:
                epfl_mfront_root = ET.parse(case.mfront_project).getroot()
                parity_mod.ensure_native_epfl_potential_exchange(
                    native_root,
                    hamaker_native,
                    mfront_parameter_defaults=mfront_defaults_for_epfl(
                        epfl_mfront_root, parity_mod
                    ),
                )
            else:
                compare_mod.set_native_hamaker(native_root, hamaker_native)
            zero_swelling_stress_rate(native_root)
            native_metrics, native_status, native_error = run_project_variant(
                root=native_root,
                source_project=case.native_project,
                prefix=native_prefix,
                project_path=native_project,
                output_dir=native_out,
                ogs_bin=args.native_ogs,
                compare_mod=compare_mod,
            )

            # 3) MFront DSM.
            mfront_prefix = f"{case.benchmark_group}_{case.case_id}_mfront_dsm".lower()
            mfront_project = tmpdir / f"{mfront_prefix}.prj"
            mfront_out = tmpdir / f"out_{mfront_prefix}"
            mfront_root = ET.parse(case.mfront_project).getroot()
            compare_mod.set_mfront_hamaker(mfront_root, hamaker_mfront)
            mfront_metrics, mfront_status, mfront_error = run_project_variant(
                root=mfront_root,
                source_project=case.mfront_project,
                prefix=mfront_prefix,
                project_path=mfront_project,
                output_dir=mfront_out,
                ogs_bin=args.mfront_ogs,
                compare_mod=compare_mod,
            )

            common = {
                "benchmark_group": case.benchmark_group,
                "case_id": case.case_id,
                "dry_density_kg_m3": float(0.5 * (dd_native + dd_mfront)),
            }
            rows.extend(
                [
                    {
                        **common,
                        "implementation": "linear",
                        "multiplier_used": math.nan,
                        "hamaker_constant_J": math.nan,
                        "run_status": linear_status,
                        "error_message": linear_error,
                        "source": case.native_project.name,
                        **linear_metrics,
                    },
                    {
                        **common,
                        "implementation": "native_dsm",
                        "multiplier_used": native_multiplier,
                        "hamaker_constant_J": hamaker_native,
                        "run_status": native_status,
                        "error_message": native_error,
                        "source": case.native_project.name,
                        **native_metrics,
                    },
                    {
                        **common,
                        "implementation": "mfront_dsm",
                        "multiplier_used": mfront_multiplier,
                        "hamaker_constant_J": hamaker_mfront,
                        "run_status": mfront_status,
                        "error_message": mfront_error,
                        "source": case.mfront_project.name,
                        **mfront_metrics,
                    },
                ]
            )
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)
    return rows


def triplet_deltas(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    by_case: dict[tuple[str, str], dict[str, dict[str, object]]] = {}
    for row in rows:
        by_case.setdefault((str(row["benchmark_group"]), str(row["case_id"])), {})[
            str(row["implementation"])
        ] = row

    out: list[dict[str, object]] = []
    for (group, case_id), impl in sorted(by_case.items()):
        linear = impl.get("linear")
        native = impl.get("native_dsm")
        mfront = impl.get("mfront_dsm")
        if linear is None or native is None or mfront is None:
            continue

        def val(r: dict[str, object], key: str) -> float:
            try:
                return float(r.get(key, math.nan))
            except Exception:
                return math.nan

        p_lin = val(linear, "final_swelling_pressure_mpa")
        p_nat = val(native, "final_swelling_pressure_mpa")
        p_mfr = val(mfront, "final_swelling_pressure_mpa")
        out.append(
            {
                "benchmark_group": group,
                "case_id": case_id,
                "dry_density_kg_m3": val(linear, "dry_density_kg_m3"),
                "linear_final_swelling_pressure_mpa": p_lin,
                "native_dsm_final_swelling_pressure_mpa": p_nat,
                "mfront_dsm_final_swelling_pressure_mpa": p_mfr,
                "delta_native_minus_linear_mpa": p_nat - p_lin,
                "delta_mfront_minus_linear_mpa": p_mfr - p_lin,
                "delta_native_minus_mfront_mpa": p_nat - p_mfr,
                "rel_native_vs_linear": abs(p_nat - p_lin) / max(abs(p_lin), 1e-12),
                "rel_mfront_vs_linear": abs(p_mfr - p_lin) / max(abs(p_lin), 1e-12),
                "native_multiplier": val(native, "multiplier_used"),
                "mfront_multiplier": val(mfront, "multiplier_used"),
                "native_to_mfront_multiplier_ratio": val(native, "multiplier_used")
                / max(val(mfront, "multiplier_used"), 1e-300),
                "native_status": str(native.get("run_status", "")),
                "mfront_status": str(mfront.get("run_status", "")),
                "linear_status": str(linear.get("run_status", "")),
            }
        )
    return out


def plot_triplet_summary(deltas: list[dict[str, object]], out_png: Path) -> None:
    if not deltas:
        return
    rows = sorted(deltas, key=lambda r: (r["benchmark_group"], r["case_id"]))
    labels = [f"{r['benchmark_group']}:{r['case_id']}" for r in rows]
    y_lin = np.array([float(r["linear_final_swelling_pressure_mpa"]) for r in rows])
    y_nat = np.array([float(r["native_dsm_final_swelling_pressure_mpa"]) for r in rows])
    y_mfr = np.array([float(r["mfront_dsm_final_swelling_pressure_mpa"]) for r in rows])

    x = np.arange(len(rows))
    w = 0.28
    fig, ax = plt.subplots(figsize=(1.1 * max(8, len(rows)), 4.8))
    ax.bar(x - w, y_lin, width=w, label="Linear reference")
    ax.bar(x, y_nat, width=w, label="DSM native")
    ax.bar(x + w, y_mfr, width=w, label="DSM MFront")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=30, ha="right")
    ax.set_ylabel("Final swelling pressure (MPa)")
    ax.set_title("Linear-law-calibrated comparison: linear vs DSM native vs DSM MFront")
    ax.grid(True, axis="y", ls=":", alpha=0.5)
    ax.legend(loc="best")
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=200)
    plt.close(fig)


def fail_if_needed(
    *,
    args: argparse.Namespace,
    calibration_rows: list[dict[str, object]],
    all_rows: list[dict[str, object]],
) -> None:
    if not args.assert_thresholds:
        return
    failures: list[str] = []
    for r in calibration_rows:
        dd = float(r["dry_density_kg_m3"])
        if float(r["native_relative_error"]) > args.rel_tol:
            failures.append(
                f"ANCHORS dd={dd:.0f}: native calibration rel err "
                f"{float(r['native_relative_error']):.6g} > {args.rel_tol:.6g}"
            )
        if float(r["mfront_relative_error"]) > args.rel_tol:
            failures.append(
                f"ANCHORS dd={dd:.0f}: mfront calibration rel err "
                f"{float(r['mfront_relative_error']):.6g} > {args.rel_tol:.6g}"
            )
    for r in all_rows:
        if str(r.get("run_status", "success")) != "success":
            failures.append(
                f"{r.get('benchmark_group')}:{r.get('case_id')}:{r.get('implementation')} failed: "
                f"{r.get('error_message', '')}"
            )
    if failures:
        raise RuntimeError(
            "Linear-law-calibrated comparison thresholds failed:\n- "
            + "\n- ".join(failures)
        )


def main() -> None:
    args = parse_args()
    args.out_root.mkdir(parents=True, exist_ok=True)

    compare_mod = load_module(
        ROOT / "run_calibrated_beacon_anchors_comparison.py",
        "run_calibrated_beacon_anchors_comparison_linear_mode",
    )
    parity_mod = load_module(
        ROOT / "run_identical_parameter_benchmark_parity.py",
        "run_identical_parameter_benchmark_parity_linear_mode",
    )
    mfront_mod = load_module(
        ANCHORS_DIR / "run_villar_dense_dd_calibration.py",
        "run_villar_dense_dd_calibration_linear_mode",
    )
    native_mod = load_module(
        ANCHORS_DIR / "run_villar_dense_dd_native_dsm_micromacro_calibration.py",
        "run_villar_dense_dd_native_dsm_micromacro_calibration_linear_mode",
    )

    calibration_rows, dd_native, mult_native, dd_mfront, mult_mfront = (
        calibration_rows_and_curves(args, mfront_mod, native_mod)
    )
    anchors_rows = anchors_triplet_rows(calibration_rows)
    benchmark_rows = run_beacon_epfl_triplet_rows(
        args=args,
        compare_mod=compare_mod,
        parity_mod=parity_mod,
        dd_curve_native=dd_native,
        mult_curve_native=mult_native,
        dd_curve_mfront=dd_mfront,
        mult_curve_mfront=mult_mfront,
    )

    all_rows = anchors_rows + benchmark_rows
    deltas = triplet_deltas(all_rows)

    calibration_csv = args.out_root / "anchors_linear_law_calibration.csv"
    rows_csv = args.out_root / "linear_law_triplet_runs_summary.csv"
    deltas_csv = args.out_root / "linear_law_triplet_pairwise_deltas.csv"
    summary_json = args.out_root / "linear_law_triplet_summary.json"
    plot_png = args.out_root / "linear_law_triplet_swelling_pressure_comparison.png"

    write_rows_csv(calibration_rows, calibration_csv)
    write_rows_csv(all_rows, rows_csv)
    write_rows_csv(deltas, deltas_csv)
    plot_triplet_summary(deltas, plot_png)

    summary = {
        "mode": "linear_law_calibrated_triplet_comparison",
        "native_ogs": str(args.native_ogs.resolve()),
        "mfront_ogs": str(args.mfront_ogs.resolve()),
        "linear_ref_ogs": str(args.linear_ref_ogs.resolve()),
        "dd_range": {
            "min": args.dd_min,
            "max": args.dd_max,
            "step": args.dd_step,
        },
        "rel_tol": args.rel_tol,
        "row_count": len(all_rows),
        "delta_count": len(deltas),
        "max_native_calibration_rel_error": max(
            (float(r["native_relative_error"]) for r in calibration_rows), default=math.nan
        ),
        "max_mfront_calibration_rel_error": max(
            (float(r["mfront_relative_error"]) for r in calibration_rows), default=math.nan
        ),
        "git_hashes": {
            "ogs_repo_hash": git_short_hash(ROOT.parent.parent.parent),
            "native_repo_hash": git_short_hash(
                Path("/Users/vinaykumar/git/ogs-native-dsm-transition")
            ),
        },
        "artifacts": {
            "calibration_csv": str(calibration_csv.resolve()),
            "rows_csv": str(rows_csv.resolve()),
            "deltas_csv": str(deltas_csv.resolve()),
            "plot_png": str(plot_png.resolve()),
        },
        "deltas": deltas,
    }
    summary_json.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    fail_if_needed(args=args, calibration_rows=calibration_rows, all_rows=all_rows)

    print(f"Wrote calibration CSV: {calibration_csv}")
    print(f"Wrote rows CSV: {rows_csv}")
    print(f"Wrote deltas CSV: {deltas_csv}")
    print(f"Wrote summary JSON: {summary_json}")
    print(f"Wrote plot: {plot_png}")


if __name__ == "__main__":
    main()
