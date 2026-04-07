#!/usr/bin/env python3
"""Run native-vs-MFront parity with identical dry-density multipliers.

This workflow enforces a single shared vdW multiplier curve for BOTH
implementations, then reruns:
1) BEACON benchmark cases (1a01, 1b) via the calibrated comparison runner.
2) ANCHORS dense dry-density points via the existing one-point calibration
   drivers (single run per dry density, no fitting).

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
import subprocess
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
ANCHORS_DIR = ROOT / "ANCHORS_MS33_ModelI"
DEFAULT_SHARED_CURVE_CSV = ANCHORS_DIR / "villar_dense_dd_calibration.csv"
OUTPUT_ROOT = ROOT / "_outputs" / "identical_parameter_native_mfront_comparison"


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
        default=Path("/Users/vinaykumar/git/build/release-native-transition-mfront/bin/ogs"),
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
        native_run = native_mod.run_native_notebook_case(
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
        ANCHORS_DIR / "run_villar_dense_dd_native_notebook_calibration.py",
        "run_villar_dense_dd_native_notebook_calibration_mod",
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

    rows = anchors_rows + beacon_rows
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
        "git_hashes": {
            "ogs_repo_hash": git_hash(ROOT.parents[2]),
            "native_repo_hash": git_hash(Path("/Users/vinaykumar/git/ogs-native-dsm-transition")),
            "materialmodels_repo_hash": git_hash(
                Path("/Users/vinaykumar/Documents/GitHub/materialmodels")
            ),
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
