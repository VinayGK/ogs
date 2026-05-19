#!/usr/bin/env python3
"""Compare BEACON 1a01 VTU stress outputs against benchmark targets.

Reads sigma from all VTU files matching:
  <prefix>_ts_*_t_*.vtu
and reports domain-mean axial/radial stresses in kPa.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from pathlib import Path

import vtk
from vtk.util.numpy_support import vtk_to_numpy


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--prefix", required=True, help="VTU prefix (without _ts_...)")
    p.add_argument(
        "--root",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="Directory containing VTU files (default: script directory)",
    )
    p.add_argument("--intermediate-time-s", type=float, default=864000.0)
    p.add_argument("--final-time-s", type=float, default=10368000.0)
    p.add_argument("--intermediate-axial-target-kpa", type=float, default=8604.0)
    p.add_argument("--intermediate-radial-target-kpa", type=float, default=9994.0)
    p.add_argument("--final-axial-target-kpa", type=float, default=2566.0)
    p.add_argument("--final-radial-target-kpa", type=float, default=3240.0)
    p.add_argument("--output-json", type=Path)
    p.add_argument("--output-csv", type=Path)
    return p.parse_args()


def read_mean_axial_radial_kpa(vtu_path: Path) -> tuple[float, float]:
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(vtu_path))
    reader.Update()
    grid = reader.GetOutput()
    sigma = vtk_to_numpy(grid.GetPointData().GetArray("sigma"))
    axial = float((-sigma[:, 1]).mean() / 1e3)
    radial = float((-sigma[:, 0]).mean() / 1e3)
    return axial, radial


def choose_row_by_time(rows: list[dict], target_time_s: float) -> dict:
    return min(rows, key=lambda r: abs(r["time_s"] - target_time_s))


def rel_err(observed: float, target: float) -> float:
    if math.isclose(target, 0.0):
        return float("nan")
    return (observed - target) / target


def main() -> int:
    args = parse_args()
    root = args.root.resolve()

    pattern = re.compile(r"_ts_(\d+)_t_([-0-9eE+.]+)\.vtu$")
    rows: list[dict] = []
    for vtu in sorted(root.glob(f"{args.prefix}_ts_*_t_*.vtu")):
        m = pattern.search(vtu.name)
        if not m:
            continue
        ts = int(m.group(1))
        t = float(m.group(2))
        axial, radial = read_mean_axial_radial_kpa(vtu)
        rows.append(
            {
                "file": vtu.name,
                "timestep": ts,
                "time_s": t,
                "time_days": t / 86400.0,
                "axial_kpa": axial,
                "radial_kpa": radial,
            }
        )

    if not rows:
        raise FileNotFoundError(f"No VTU files found for prefix: {args.prefix}")

    rows.sort(key=lambda r: (r["time_s"], r["timestep"]))

    inter = choose_row_by_time(rows, args.intermediate_time_s)
    final = choose_row_by_time(rows, args.final_time_s)

    targets = {
        "intermediate": {
            "axial_kpa": args.intermediate_axial_target_kpa,
            "radial_kpa": args.intermediate_radial_target_kpa,
            "time_s": args.intermediate_time_s,
        },
        "final": {
            "axial_kpa": args.final_axial_target_kpa,
            "radial_kpa": args.final_radial_target_kpa,
            "time_s": args.final_time_s,
        },
    }

    summary = {
        "prefix": args.prefix,
        "targets": targets,
        "comparison": {
            "intermediate": {
                "matched_output": inter,
                "axial_rel_error": rel_err(inter["axial_kpa"], args.intermediate_axial_target_kpa),
                "radial_rel_error": rel_err(inter["radial_kpa"], args.intermediate_radial_target_kpa),
            },
            "final": {
                "matched_output": final,
                "axial_rel_error": rel_err(final["axial_kpa"], args.final_axial_target_kpa),
                "radial_rel_error": rel_err(final["radial_kpa"], args.final_radial_target_kpa),
            },
        },
        "time_series": rows,
    }

    json_path = args.output_json or (root / f"{args.prefix}_targets_summary.json")
    csv_path = args.output_csv or (root / f"{args.prefix}_timeseries.csv")

    json_path.write_text(json.dumps(summary, indent=2) + "\n")

    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["file", "timestep", "time_s", "time_days", "axial_kpa", "radial_kpa"],
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote JSON: {json_path}")
    print(f"Wrote CSV: {csv_path}")
    print(
        "Intermediate: "
        f"t={inter['time_s']:.0f}s ({inter['time_days']:.2f}d), "
        f"axial={inter['axial_kpa']:.3f} kPa, radial={inter['radial_kpa']:.3f} kPa"
    )
    print(
        "Final: "
        f"t={final['time_s']:.0f}s ({final['time_days']:.2f}d), "
        f"axial={final['axial_kpa']:.3f} kPa, radial={final['radial_kpa']:.3f} kPa"
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
