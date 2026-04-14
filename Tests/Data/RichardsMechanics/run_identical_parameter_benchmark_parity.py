#!/usr/bin/env python3
"""Run native-vs-MFront parity with one identical shared parameter set.

This workflow enforces a single shared vdW multiplier curve for BOTH
implementations, then reruns:
1) BEACON benchmark cases (1a01, 1b) via the calibrated comparison runner.
2) ANCHORS dense dry-density points via the existing one-point calibration
   drivers (single run per dry density, no fitting).
3) EPFL/BGR WP3 AB' benchmark with MCC + dsm_micromacro micro-potential path
   for both native and MFront implementations.

Outputs:
- rows CSV
- pairwise-deltas CSV
- summary JSON
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import math
import shutil
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
ANCHORS_DIR = ROOT / "ANCHORS_MS33_ModelI"
EPFL_DIR = ROOT / "BEACON_WP3_BGR_EPFL"
EPFL_NATIVE_PROJECT = EPFL_DIR / "bgr_wp3_p2_1_abprime_native.prj"
EPFL_MFRONT_PROJECT = EPFL_DIR / "bgr_wp3_p2_1_abprime_mfront.prj"
DEFAULT_SHARED_CURVE_CSV = ANCHORS_DIR / "villar_dense_dd_calibration.csv"
OUTPUT_ROOT = ROOT / "_outputs" / "identical_parameter_native_mfront_comparison_with_epfl"


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
            "Run BEACON+ANCHORS native-vs-MFront parity with identical "
            "dry-density multiplier curve."
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
        "--shared-curve-csv",
        type=Path,
        default=DEFAULT_SHARED_CURVE_CSV,
        help="CSV with columns dry_density_kg_m3, vdw_multiplier.",
    )
    parser.add_argument(
        "--out-csv",
        type=Path,
        default=OUTPUT_ROOT / "identical_parameter_benchmark_runs_summary.csv",
    )
    parser.add_argument(
        "--out-deltas-csv",
        type=Path,
        default=OUTPUT_ROOT / "identical_parameter_benchmark_pairwise_deltas.csv",
    )
    parser.add_argument(
        "--out-json",
        type=Path,
        default=OUTPUT_ROOT / "identical_parameter_benchmark_runs_summary.json",
    )
    return parser.parse_args()


def load_multiplier_curve(path: Path) -> tuple[np.ndarray, np.ndarray]:
    with path.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    dd = np.array([float(r["dry_density_kg_m3"]) for r in rows], dtype=float)
    mult = np.array([float(r["vdw_multiplier"]) for r in rows], dtype=float)
    order = np.argsort(dd)
    return dd[order], mult[order]


def interpolate_multiplier(dd: float, dd_curve: np.ndarray, mult_curve: np.ndarray) -> float:
    return float(np.interp(dd, dd_curve, mult_curve, left=mult_curve[0], right=mult_curve[-1]))


def get_parameter_value(root: ET.Element, name: str, default: float | None = None) -> float:
    node = root.find(f"./parameters/parameter[name='{name}']/value")
    if node is None or node.text is None:
        if default is None:
            raise RuntimeError(f"Missing parameter '{name}' in project.")
        return default
    return float(node.text.strip())


def _get_property_node(root: ET.Element, name: str) -> ET.Element | None:
    return root.find(f"./media/medium/properties/property[name='{name}']")


def _get_property_float(node: ET.Element, tag: str, default: float) -> float:
    child = node.find(tag)
    if child is None or child.text is None:
        return default
    return float(child.text.strip())


def estimate_initial_micro_water_content(root: ET.Element) -> float:
    """Estimate n_l0 from EPFL-style native inputs when n_l0 is not explicit."""
    phi0 = get_parameter_value(root, "phi0")
    phi_tr0 = get_parameter_value(root, "phi_tr0", default=0.0)
    micro_porosity0 = max(1e-12, phi0 - phi_tr0)
    p_l0 = get_parameter_value(root, "pressure_ic", default=0.0)
    p_cap = max(0.0, -p_l0)

    sat_prop = _get_property_node(root, "saturation_micro")
    if sat_prop is None:
        return 0.99 * micro_porosity0

    sat_type = sat_prop.find("type")
    sat_type_text = sat_type.text.strip() if sat_type is not None and sat_type.text else ""
    if sat_type_text != "SaturationVanGenuchten":
        return 0.99 * micro_porosity0

    s_lr = _get_property_float(sat_prop, "residual_liquid_saturation", 0.0)
    s_gr = _get_property_float(sat_prop, "residual_gas_saturation", 0.0)
    m = _get_property_float(sat_prop, "exponent", 0.5)
    p_b = _get_property_float(sat_prop, "p_b", 1.0)
    if p_cap <= 0.0:
        s_eff = 1.0
    else:
        m = min(max(m, 1e-12), 1.0 - 1e-12)
        p_b = max(p_b, 1e-12)
        n = 1.0 / (1.0 - m)
        s_eff = (1.0 + (p_cap / p_b) ** n) ** (-m)

    sat = s_lr + (1.0 - s_lr - s_gr) * s_eff
    sat = min(max(sat, 0.0), 1.0)
    return micro_porosity0 * sat


def set_or_create(parent: ET.Element, tag: str, value: str) -> ET.Element:
    node = parent.find(tag)
    if node is None:
        node = ET.SubElement(parent, tag)
    node.text = value
    return node


def ensure_native_epfl_potential_exchange(
    root: ET.Element,
    hamaker_j: float,
    mfront_parameter_defaults: dict[str, float] | None = None,
) -> None:
    """Inject a native potential_exchange block mapped from EPFL MFront params."""
    process = root.find("./processes/process")
    if process is None:
        raise RuntimeError("Missing process block in EPFL native project.")

    defaults = mfront_parameter_defaults or {}

    # Pull shared constitutive constants from parameters to keep a single set.
    phi0 = get_parameter_value(root, "phi0")
    rho_sr = get_parameter_value(
        root,
        "ReferenceDensitySolid",
        default=defaults.get("ReferenceDensitySolid", 2780.0),
    )
    specific_surface = get_parameter_value(
        root,
        "SpecificSurface",
        default=defaults.get("SpecificSurface", 523.0),
    )
    try:
        n_l0 = get_parameter_value(root, "n_l0")
    except RuntimeError:
        n_l0 = defaults.get("n_l0", estimate_initial_micro_water_content(root))
    micro_swelling_slope = get_parameter_value(
        root,
        "MicroSwellingStrainSlope",
        default=defaults.get("MicroSwellingStrainSlope", 0.0),
    )
    mass_exchange_coeff = get_parameter_value(
        root,
        "MassExchangeCoefficient",
        default=defaults.get("MassExchangeCoefficient", 5e-15),
    )

    # Keep mass exchange coefficient explicitly aligned with the shared set.
    micro_porosity = process.find("micro_porosity")
    if micro_porosity is None:
        micro_porosity = ET.SubElement(process, "micro_porosity")
    set_or_create(
        micro_porosity,
        "mass_exchange_coefficient",
        f"{mass_exchange_coeff:.16g}",
    )

    potential_exchange = process.find("potential_exchange")
    if potential_exchange is None:
        potential_exchange = ET.SubElement(process, "potential_exchange")

    set_or_create(potential_exchange, "enabled", "true")
    set_or_create(potential_exchange, "mode", "full_potential")
    set_or_create(potential_exchange, "pressure_tolerance", "1e-12")
    set_or_create(potential_exchange, "hamaker_constant", f"{hamaker_j:.16g}")
    set_or_create(potential_exchange, "specific_surface", f"{specific_surface:.16g}")
    set_or_create(
        potential_exchange,
        "micro_solid_density_reference",
        f"{rho_sr:.16g}",
    )
    # Reference split for native exchange path; use total-solid fraction.
    set_or_create(
        potential_exchange,
        "micro_solid_volume_fraction_reference",
        f"{max(1e-12, 1.0 - phi0):.16g}",
    )
    set_or_create(
        potential_exchange,
        "initial_micro_water_content",
        f"{max(1e-12, n_l0):.16g}",
    )
    set_or_create(
        potential_exchange,
        "local_nonlinear_solve_mode",
        "scalar_microstate_storage_mode",
    )
    set_or_create(potential_exchange, "fd_jacobian_for_exchange", "false")
    set_or_create(
        potential_exchange,
        "micro_potential_convention",
        "negative_attractive",
    )
    set_or_create(
        potential_exchange,
        "micro_solid_volume_fraction_mode",
        "reference",
    )
    set_or_create(
        potential_exchange,
        "micro_water_content_swelling_slope",
        f"{max(0.0, micro_swelling_slope):.16g}",
    )


def run_epfl_rows(
    compare_mod,
    native_ogs: Path,
    mfront_ogs: Path,
    dd_curve: np.ndarray,
    mult_curve: np.ndarray,
) -> list[dict[str, object]]:
    """Run EPFL native and MFront decks with one shared DD-based multiplier."""
    if not EPFL_NATIVE_PROJECT.exists() or not EPFL_MFRONT_PROJECT.exists():
        raise FileNotFoundError("EPFL benchmark projects are missing.")

    dd_native = compare_mod.extract_dd_from_project(EPFL_NATIVE_PROJECT)
    dd_mfront = compare_mod.extract_dd_from_project(EPFL_MFRONT_PROJECT)
    dd_epfl = 0.5 * (dd_native + dd_mfront)
    multiplier = interpolate_multiplier(dd_epfl, dd_curve, mult_curve)
    hamaker_j = compare_mod.HAMAKER_REFERENCE_J * multiplier

    tmpdir = Path(tempfile.mkdtemp(prefix="identical-epfl-runs-"))
    try:
        native_prefix = "epfl_abprime_native_identical"
        mfront_prefix = "epfl_abprime_mfront_identical"

        native_prj = tmpdir / f"{native_prefix}.prj"
        mfront_prj = tmpdir / f"{mfront_prefix}.prj"
        native_out = tmpdir / native_prefix
        mfront_out = tmpdir / mfront_prefix

        mfront_root = ET.parse(EPFL_MFRONT_PROJECT).getroot()
        mfront_defaults = {
            "ReferenceDensitySolid": get_parameter_value(
                mfront_root, "ReferenceDensitySolid", default=2780.0
            ),
            "SpecificSurface": get_parameter_value(
                mfront_root, "SpecificSurface", default=523.0
            ),
            "MicroSwellingStrainSlope": get_parameter_value(
                mfront_root, "MicroSwellingStrainSlope", default=0.0
            ),
            "MassExchangeCoefficient": get_parameter_value(
                mfront_root, "MassExchangeCoefficient", default=5e-15
            ),
            "n_l0": get_parameter_value(mfront_root, "n_l0", default=1e-12),
        }

        native_root = ET.parse(EPFL_NATIVE_PROJECT).getroot()
        compare_mod.absolutize_mesh_and_geometry(native_root, EPFL_NATIVE_PROJECT.parent)
        compare_mod.set_output_prefix(native_root, native_prefix)
        compare_mod.set_biot_coefficient(native_root, 1.0)
        compare_mod.set_nonlinear_max_iter(native_root, 80)
        compare_mod.relax_first_component_abstol(native_root, 5e-7)
        ensure_native_epfl_potential_exchange(
            native_root, hamaker_j, mfront_parameter_defaults=mfront_defaults
        )
        ET.ElementTree(native_root).write(
            native_prj, encoding="ISO-8859-1", xml_declaration=True
        )

        compare_mod.absolutize_mesh_and_geometry(mfront_root, EPFL_MFRONT_PROJECT.parent)
        compare_mod.set_output_prefix(mfront_root, mfront_prefix)
        compare_mod.set_biot_coefficient(mfront_root, 1.0)
        compare_mod.set_nonlinear_max_iter(mfront_root, 80)
        compare_mod.relax_first_component_abstol(mfront_root, 5e-7)
        compare_mod.set_mfront_hamaker(mfront_root, hamaker_j)
        ET.ElementTree(mfront_root).write(
            mfront_prj, encoding="ISO-8859-1", xml_declaration=True
        )

        compare_mod.run_ogs(native_ogs, native_prj, native_out)
        compare_mod.run_ogs(mfront_ogs, mfront_prj, mfront_out)

        native_metrics = compare_mod.evaluate_beacon_metrics(
            compare_mod.latest_vtu(native_out, native_prefix)
        )
        mfront_metrics = compare_mod.evaluate_beacon_metrics(
            compare_mod.latest_vtu(mfront_out, mfront_prefix)
        )

        common = {
            "benchmark_group": "EPFL_BGR_WP3",
            "case_id": "p2_1_abprime",
            "dry_density_kg_m3": dd_epfl,
            "multiplier_used": multiplier,
            "hamaker_constant_J": hamaker_j,
            "source": "BEACON_WP3_BGR_EPFL",
            "run_status": "success",
        }
        native_row = {
            **common,
            "implementation": "native",
            **native_metrics,
        }
        mfront_row = {
            **common,
            "implementation": "mfront",
            **mfront_metrics,
        }
        return [native_row, mfront_row]
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def write_rows_csv(rows: list[dict[str, object]], path: Path) -> None:
    keys = sorted({k for row in rows for k in row.keys()})
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def git_hash(repo: Path) -> str:
    try:
        return (
            subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], cwd=repo, text=True)
            .strip()
        )
    except Exception:
        return ""


def anchors_rows_identical_multiplier(
    native_mod,
    mfront_mod,
    native_ogs: Path,
    mfront_ogs: Path,
    dd_curve: np.ndarray,
    mult_curve: np.ndarray,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for i, dd in enumerate(dd_curve):
        multiplier = interpolate_multiplier(float(dd), dd_curve, mult_curve)

        mfront_case = mfront_mod.Case(float(dd))
        native_case = native_mod.Case(float(dd))
        mfront_nl0 = mfront_mod.n_l0_from_micro_suction(
            mfront_case.phi0, mfront_mod.HAMAKER_LITERATURE
        )
        native_nl0 = native_mod.n_l0_from_micro_suction(
            native_case.phi0, native_mod.HAMAKER_LITERATURE
        )

        mfront_run = mfront_mod.run_mfront_case(
            mfront_ogs,
            mfront_case,
            multiplier,
            mfront_nl0,
            900 + i,
        )
        native_run = native_mod.run_native_dsm_micromacro_case(
            native_ogs,
            native_case,
            multiplier,
            native_nl0,
            900 + i,
        )

        mfront_ps = float(mfront_run["pressure_mpa"])
        native_ps = float(native_run["pressure_mpa"])

        rows.append(
            {
                "benchmark_group": "ANCHORS_MS33_ModelI",
                "case_id": f"dd{int(dd)}",
                "implementation": "mfront",
                "dry_density_kg_m3": float(dd),
                "multiplier_used": multiplier,
                "hamaker_constant_J": float(mfront_run["hamaker_effective"]),
                "final_total_stress_mpa": mfront_ps,
                "final_swelling_pressure_mpa": mfront_ps,
                "axial_sigma_kpa": mfront_ps * 1e3,
                "radial_sigma_kpa": math.nan,
                "axial_swelling_stress_kpa": math.nan,
                "radial_swelling_stress_kpa": math.nan,
                "axial_displacement_mm": math.nan,
                "radial_displacement_mm": math.nan,
                "dry_density_mean_kg_m3": float(dd),
                "run_status": "success",
                "source": "anchors_identical_parameter_rerun",
            }
        )
        rows.append(
            {
                "benchmark_group": "ANCHORS_MS33_ModelI",
                "case_id": f"dd{int(dd)}",
                "implementation": "native",
                "dry_density_kg_m3": float(dd),
                "multiplier_used": multiplier,
                "hamaker_constant_J": float(native_run["hamaker_effective_J"]),
                "final_total_stress_mpa": native_ps,
                "final_swelling_pressure_mpa": native_ps,
                "axial_sigma_kpa": native_ps * 1e3,
                "radial_sigma_kpa": math.nan,
                "axial_swelling_stress_kpa": math.nan,
                "radial_swelling_stress_kpa": math.nan,
                "axial_displacement_mm": math.nan,
                "radial_displacement_mm": math.nan,
                "dry_density_mean_kg_m3": float(dd),
                "run_status": "success",
                "source": "anchors_identical_parameter_rerun",
            }
        )
    return rows


def max_abs_rel_delta(deltas: list[dict[str, object]], group: str) -> dict[str, object]:
    filtered = [d for d in deltas if str(d["benchmark_group"]) == group]
    if not filtered:
        return {}

    def abs_delta(row: dict[str, object]) -> float:
        return abs(float(row["delta_swelling_pressure_mpa_native_minus_mfront"]))

    def rel_delta(row: dict[str, object]) -> float:
        ref = abs(float(row["mfront_final_swelling_pressure_mpa"]))
        return abs_delta(row) / max(ref, 1e-12)

    max_abs_row = max(filtered, key=abs_delta)
    max_rel_row = max(filtered, key=rel_delta)
    return {
        "max_abs_case_id": max_abs_row["case_id"],
        "max_abs_delta_mpa": abs_delta(max_abs_row),
        "max_rel_case_id": max_rel_row["case_id"],
        "max_rel_delta": rel_delta(max_rel_row),
    }


def main() -> None:
    args = parse_args()

    compare_mod = load_module(
        ROOT / "run_calibrated_beacon_anchors_comparison.py",
        "run_calibrated_beacon_anchors_comparison_mod",
    )
    mfront_mod = load_module(
        ANCHORS_DIR / "run_villar_dense_dd_calibration.py",
        "run_villar_dense_dd_calibration_mod",
    )
    native_mod = load_module(
        ANCHORS_DIR / "run_villar_dense_dd_native_dsm_micromacro_calibration.py",
        "run_villar_dense_dd_native_dsm_micromacro_calibration_mod",
    )

    if not args.shared_curve_csv.exists():
        raise FileNotFoundError(args.shared_curve_csv)
    dd_curve, mult_curve = load_multiplier_curve(args.shared_curve_csv)

    # BEACON: enforce identical multiplier curve in both implementations.
    beacon_rows = compare_mod.run_beacon_rows(
        native_ogs=args.native_ogs.resolve(),
        mfront_ogs=args.mfront_ogs.resolve(),
        native_curve=(dd_curve, mult_curve),
        mfront_curve=(dd_curve, mult_curve),
    )

    # ANCHORS: rerun both implementations with the same multiplier per DD.
    anchors_rows = anchors_rows_identical_multiplier(
        native_mod=native_mod,
        mfront_mod=mfront_mod,
        native_ogs=args.native_ogs.resolve(),
        mfront_ogs=args.mfront_ogs.resolve(),
        dd_curve=dd_curve,
        mult_curve=mult_curve,
    )

    # EPFL: rerun native and MFront with the same DD-derived multiplier and MCC setup.
    epfl_rows = run_epfl_rows(
        compare_mod=compare_mod,
        native_ogs=args.native_ogs.resolve(),
        mfront_ogs=args.mfront_ogs.resolve(),
        dd_curve=dd_curve,
        mult_curve=mult_curve,
    )

    rows = anchors_rows + beacon_rows + epfl_rows
    rows.sort(key=lambda r: (str(r["benchmark_group"]), str(r["case_id"]), str(r["implementation"])))
    deltas = compare_mod.pairwise_deltas(rows)
    deltas.sort(key=lambda d: (str(d["benchmark_group"]), str(d["case_id"])))

    write_rows_csv(rows, args.out_csv)
    write_rows_csv(deltas, args.out_deltas_csv)

    summary = {
        "mode": "identical_parameter_shared_multiplier_curve",
        "shared_curve_csv": str(args.shared_curve_csv.resolve()),
        "row_count": len(rows),
        "pairwise_delta_count": len(deltas),
        "beacon_delta_stats": max_abs_rel_delta(deltas, "BEACON_report"),
        "anchors_delta_stats": max_abs_rel_delta(deltas, "ANCHORS_MS33_ModelI"),
        "epfl_delta_stats": max_abs_rel_delta(deltas, "EPFL_BGR_WP3"),
        "git_hashes": {
            "ogs_repo_hash": git_hash(ROOT.parents[2]),
            "native_repo_hash": git_hash(Path("/Users/vinaykumar/git/ogs-native-dsm-transition")),
            "materialmodels_repo_hash": "",
        },
        "rows": rows,
        "pairwise_deltas": deltas,
    }
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    print(f"Wrote rows CSV: {args.out_csv}")
    print(f"Wrote pairwise deltas CSV: {args.out_deltas_csv}")
    print(f"Wrote summary JSON: {args.out_json}")


if __name__ == "__main__":
    main()
