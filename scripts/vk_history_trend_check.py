#!/usr/bin/env python3
"""Check trend-level agreement between two notebook-role history CSV files.

The notebook-role history extraction snapshot is meant to be a calibration
baseline, not a brittle exact-value regression. This helper therefore verifies
the invariants that matter for the current comparison step:
- identical time stamps,
- monotone growth of the micro water content / porosity proxy,
- positive micro exchange source throughout,
- growing stress magnitude,
- and matching trend directions between the committed snapshot and the current
  extracted history.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def load_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def as_float_rows(rows: list[dict[str, str]], key: str) -> list[float]:
    return [float(row[key]) for row in rows]


def is_monotone_non_decreasing(values: list[float], tol: float = 0.0) -> bool:
    return all(b >= a - tol for a, b in zip(values, values[1:]))


def trend_label(values: list[float], tol: float = 0.0) -> str:
    delta = values[-1] - values[0]
    if abs(delta) <= tol:
        return "flat"
    return "inc" if delta > 0.0 else "dec"


def compare_times(reference: list[dict[str, str]], actual: list[dict[str, str]]) -> None:
    if len(reference) != len(actual):
        raise SystemExit(
            f"Row count mismatch: reference has {len(reference)} rows, actual has "
            f"{len(actual)} rows."
        )

    ref_times = as_float_rows(reference, "time_s")
    act_times = as_float_rows(actual, "time_s")
    if len(ref_times) < 2:
        raise SystemExit("Need at least two time points for trend comparison.")

    for idx, (ref_t, act_t) in enumerate(zip(ref_times, act_times)):
        if not math.isclose(ref_t, act_t, rel_tol=0.0, abs_tol=1e-12):
            raise SystemExit(
                f"Time stamp mismatch at row {idx}: reference={ref_t}, actual={act_t}."
            )


def compare_trends(reference: list[dict[str, str]], actual: list[dict[str, str]]) -> None:
    trend_columns = [
        ("pressure_mean", "inc"),
        ("micro_pressure_mean", "dec"),
        ("micro_saturation_mean", "inc"),
        ("vk_micro_water_content_mean", "inc"),
        ("vk_micro_porosity_mean", "inc"),
        ("vk_micro_exchange_source_mean", "dec"),
        ("swelling_stress_c0_mean", "dec"),
        ("sigma_c0_mean", "dec"),
    ]

    for column, expected_reference_trend in trend_columns:
        ref_values = as_float_rows(reference, column)
        act_values = as_float_rows(actual, column)
        ref_trend = trend_label(ref_values, tol=1e-15)
        act_trend = trend_label(act_values, tol=1e-15)
        print(
            f"{column}: ref {ref_values[0]:.16g} -> {ref_values[-1]:.16g} "
            f"({ref_trend}), act {act_values[0]:.16g} -> {act_values[-1]:.16g} "
            f"({act_trend})"
        )
        if ref_trend != expected_reference_trend:
            raise SystemExit(
                f"Reference CSV trend for '{column}' changed: expected "
                f"{expected_reference_trend}, got {ref_trend}."
            )
        if act_trend != ref_trend:
            raise SystemExit(
                f"Trend mismatch for '{column}': reference {ref_trend}, actual {act_trend}."
            )

    actual_n_l = as_float_rows(actual, "vk_micro_water_content_mean")
    actual_phi_m = as_float_rows(actual, "vk_micro_porosity_mean")
    actual_exchange = as_float_rows(actual, "vk_micro_exchange_source_mean")
    actual_sigma = as_float_rows(actual, "sigma_c0_mean")

    if not is_monotone_non_decreasing(actual_n_l, tol=1e-15):
        raise SystemExit("Actual vk_micro_water_content_mean is not monotone non-decreasing.")
    if not is_monotone_non_decreasing(actual_phi_m, tol=1e-15):
        raise SystemExit("Actual vk_micro_porosity_mean is not monotone non-decreasing.")
    if any(value <= 0.0 for value in actual_exchange):
        raise SystemExit("Actual vk_micro_exchange_source_mean is not strictly positive.")
    if not (abs(actual_sigma[-1]) > abs(actual_sigma[0])):
        raise SystemExit("Actual sigma_c0_mean magnitude did not grow.")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Compare notebook-role history CSVs on trend rather than exact values."
    )
    parser.add_argument("--reference", required=True, help="Committed reference CSV.")
    parser.add_argument("--actual", required=True, help="Extracted history CSV.")
    args = parser.parse_args()

    reference = load_csv(Path(args.reference))
    actual = load_csv(Path(args.actual))
    compare_times(reference, actual)
    compare_trends(reference, actual)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
