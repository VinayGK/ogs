#!/usr/bin/env python3
"""
DSM parity runner and reporter.

Usage:
    python3 run_dsm_parity.py [--suite SUITE] [--outdir DIR] [--no-run]

Options:
    --suite SUITE   Run only the named suite (e.g. "dsm_ms33"). Default: all.
    --outdir DIR    Override output directory (default: /tmp/dsm_parity_runs).
    --no-run        Skip simulation; compare existing output in --outdir.

Suites are defined in PARITY_SUITES at the bottom of this file.
To add a new model pair, append an entry there — no other changes needed.
"""

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

# ── VTK reader (requires vtk Python package, same one OGS ships with) ────────
try:
    import vtk
    from vtk.util.numpy_support import vtk_to_numpy
except ImportError:
    sys.exit("ERROR: vtk Python package not found. "
             "Install with: pip install vtk  or use the OGS venv.")


# ─────────────────────────────────────────────────────────────────────────────
# Build / repo roots — edit these if your layout differs
# ─────────────────────────────────────────────────────────────────────────────
REPO_ROOT   = Path(__file__).resolve().parents[1]
NATIVE_OGS  = Path("/Users/vinaykumar/git/build/native-release-omp-sharedcache/bin/ogs")
MFRONT_OGS  = Path("/Users/vinaykumar/git/build/mfront-release-omp-sharedcache/bin/ogs")
PRJ_DIR     = REPO_ROOT / "Tests/Data/RichardsMechanics/ANCHORS_MS33_StrictParity"
DEFAULT_OUT = Path("/tmp/dsm_parity_runs")


# ─────────────────────────────────────────────────────────────────────────────
# Parity suite registry
#
# Each entry is a dict with:
#   name         – short identifier used with --suite
#   description  – one-line human description
#   native_prj   – path to native PRJ (relative to PRJ_DIR, or absolute)
#   mfront_prj   – path to mfront PRJ
#   native_bin   – OGS binary for native run  (Path or None → NATIVE_OGS)
#   mfront_bin   – OGS binary for mfront run  (Path or None → MFRONT_OGS)
#   field_map    – {native_field: mfront_field} pairs to compare
#   notes        – optional string shown in the report header
#
# ─────────────────────────────────────────────────────────────────────────────
PARITY_SUITES = [
    {
        "name": "dsm_ms33",
        "description": "MS33 DSM micro-macro bridge — vdW potential exchange parity",
        "native_prj":  PRJ_DIR / "ms33_dsm_parity_native.prj",
        "mfront_prj":  PRJ_DIR / "ms33_dsm_parity_mfront.prj",
        "native_bin":  None,   # → NATIVE_OGS
        "mfront_bin":  None,   # → MFRONT_OGS
        "field_map": {
            "micro_water_content":  "n_l",
            "micro_porosity":       "phi_m",
            "micro_exchange_source":"rho_l_hat",
            "sigma":                "sigma",
            "displacement":         "displacement",
            "pressure":             "pressure",
            "saturation":           "saturation",
        },
        "notes": (
            "n_l/phi_m should reach machine-epsilon parity at ts≥60 (suction≤50 MPa). "
            "sigma early-time gap (~20%) is a known micro-EOS initialisation difference."
        ),
    },
    # ── ADD NEW SUITES HERE ───────────────────────────────────────────────────
    # {
    #     "name": "dsm_ms33_mcc",
    #     "description": "MS33 DSM MCC variant — ModCamClay bridge parity",
    #     "native_prj":  PRJ_DIR / "ms33_dsm_mcc_parity_native.prj",
    #     "mfront_prj":  PRJ_DIR / "ms33_dsm_mcc_parity_mfront.prj",
    #     "native_bin":  None,
    #     "mfront_bin":  None,
    #     "field_map": {
    #         "micro_water_content":  "n_l",
    #         "micro_porosity":       "phi_m",
    #         "micro_exchange_source":"rho_l_hat",
    #         "sigma":                "sigma",
    #         "displacement":         "displacement",
    #     },
    #     "notes": "MCC variant — also compare pre_consolidation_pressure.",
    # },
]


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def run_ogs(binary: Path, prj: Path, outdir: Path) -> float:
    """Run OGS and return elapsed wall-time in seconds."""
    outdir.mkdir(parents=True, exist_ok=True)
    cmd = [str(binary), str(prj), "-o", str(outdir)]
    print(f"  Running: {binary.name}  {prj.name}  → {outdir}")
    t0 = time.perf_counter()
    result = subprocess.run(cmd, capture_output=True, text=True)
    elapsed = time.perf_counter() - t0
    if result.returncode != 0:
        print("  STDOUT:", result.stdout[-2000:])
        print("  STDERR:", result.stderr[-2000:])
        raise RuntimeError(f"OGS run failed (rc={result.returncode}): {prj.name}")
    # Extract accepted/rejected step counts from OGS log
    for line in result.stdout.splitlines():
        if "accepted steps" in line:
            print(f"  {line.strip()}")
    print(f"  Completed in {elapsed:.2f}s")
    return elapsed


def read_vtu(path: Path) -> dict:
    """Return {field_name: np.ndarray} for all point and cell data arrays."""
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(path))
    reader.Update()
    ug = reader.GetOutput()
    fields = {}
    for src in [ug.GetPointData(), ug.GetCellData()]:
        for i in range(src.GetNumberOfArrays()):
            arr = src.GetArray(i)
            fields[arr.GetName()] = vtk_to_numpy(arr).astype(float)
    return fields


def vtu_files(directory: Path) -> list[Path]:
    """Sorted list of .vtu files (timestep outputs only, not pvd/pvtu)."""
    files = sorted(directory.glob("*.vtu"))
    # Exclude mesh/geo files that don't follow the ts_NNN naming convention
    return [f for f in files if "_ts_" in f.name]


def extract_ts(path: Path) -> int:
    """Extract integer timestep number from OGS VTU filename."""
    stem = path.stem
    ts_part = stem.split("_ts_")[1].split("_")[0]
    return int(ts_part)


def compare_suite(suite: dict, native_outdir: Path, mfront_outdir: Path) -> dict:
    """
    Compare native vs mfront VTU outputs for one suite.

    Returns:
        {field_native: {"ts": [...], "mae": [...], "max_abs": [...], "max_rel": [...]}}
    """
    native_vtus  = vtu_files(native_outdir)
    mfront_vtus  = vtu_files(mfront_outdir)

    if len(native_vtus) != len(mfront_vtus):
        print(f"  WARNING: file count mismatch "
              f"({len(native_vtus)} native, {len(mfront_vtus)} mfront)")

    field_map = suite["field_map"]
    results = {nf: {"ts": [], "mae": [], "max_abs": [], "max_rel": []}
               for nf in field_map}

    for nf_path, mf_path in zip(native_vtus, mfront_vtus):
        ts = extract_ts(nf_path)
        n_data = read_vtu(nf_path)
        m_data = read_vtu(mf_path)

        for nfield, mfield in field_map.items():
            na = n_data.get(nfield)
            ma = m_data.get(mfield)
            if na is None or ma is None:
                continue
            diff = np.abs(na.ravel() - ma.ravel())
            scale = np.abs(na).max()
            results[nfield]["ts"].append(ts)
            results[nfield]["mae"].append(diff.mean())
            results[nfield]["max_abs"].append(diff.max())
            results[nfield]["max_rel"].append(
                diff.max() / scale if scale > 1e-30 else 0.0
            )

    return results


def print_report(suite: dict, results: dict,
                 native_time: float, mfront_time: float):
    """Print a formatted parity report for one suite."""
    width = 78
    print()
    print("=" * width)
    print(f"  Suite : {suite['name']}")
    print(f"  Desc  : {suite['description']}")
    if suite.get("notes"):
        print(f"  Notes : {suite['notes']}")
    print(f"  Times : native={native_time:.1f}s  mfront={mfront_time:.1f}s")
    print("=" * width)

    field_map = suite["field_map"]
    col = 24

    # Per-field summary: worst early (ts≤50) and worst late (ts≥60)
    print(f"\n  {'Field (native→mfront)':<{col}}  "
          f"{'MAE early':>12}  {'MAE late':>10}  "
          f"{'maxRel early':>13}  {'maxRel late':>11}")
    print("  " + "-" * (col + 54))

    all_ok = True
    for nfield, mfield in field_map.items():
        r = results.get(nfield)
        if not r or not r["ts"]:
            print(f"  {nfield:<{col}}  (no data)")
            continue

        tsa   = np.array(r["ts"])
        mae   = np.array(r["mae"])
        mxrel = np.array(r["max_rel"])

        early = tsa <= 50
        late  = tsa >= 60

        mae_e   = mae[early].max()   if early.any() else float("nan")
        mae_l   = mae[late].max()    if late.any()  else float("nan")
        mxrel_e = mxrel[early].max() if early.any() else float("nan")
        mxrel_l = mxrel[late].max()  if late.any()  else float("nan")

        label = f"{nfield} → {mfield}"
        flag = ""
        # Flag anomalies: late max_rel > 1e-6 for primary physics fields
        if nfield in ("micro_water_content", "micro_porosity") and mxrel_l > 1e-6:
            flag = "  ← UNEXPECTED"
            all_ok = False

        print(f"  {label:<{col}}  "
              f"{mae_e:>12.3e}  {mae_l:>10.3e}  "
              f"{mxrel_e:>13.3e}  {mxrel_l:>11.3e}{flag}")

    print()
    if all_ok:
        print("  RESULT: ✓  Primary physics variables within expected tolerance.")
    else:
        print("  RESULT: ✗  One or more fields exceed expected tolerance — investigate.")
    print("=" * width)


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--suite",  default=None,
                        help="Run only this suite by name.")
    parser.add_argument("--outdir", default=str(DEFAULT_OUT),
                        help=f"Root output directory (default: {DEFAULT_OUT}).")
    parser.add_argument("--no-run", action="store_true",
                        help="Skip simulation; compare existing output.")
    args = parser.parse_args()

    outroot = Path(args.outdir)
    suites  = PARITY_SUITES
    if args.suite:
        suites = [s for s in PARITY_SUITES if s["name"] == args.suite]
        if not suites:
            names = [s["name"] for s in PARITY_SUITES]
            sys.exit(f"Unknown suite '{args.suite}'. Available: {names}")

    for suite in suites:
        print(f"\n{'─'*78}")
        print(f"  Suite: {suite['name']}  —  {suite['description']}")
        print(f"{'─'*78}")

        native_bin  = suite.get("native_bin")  or NATIVE_OGS
        mfront_bin  = suite.get("mfront_bin") or MFRONT_OGS
        native_prj  = Path(suite["native_prj"])
        mfront_prj  = Path(suite["mfront_prj"])

        native_outdir  = outroot / suite["name"] / "native"
        mfront_outdir  = outroot / suite["name"] / "mfront"

        native_time  = 0.0
        mfront_time  = 0.0

        if not args.no_run:
            for binary, prj, label in [
                (native_bin,  native_prj,  "native"),
                (mfront_bin,  mfront_prj,  "mfront"),
            ]:
                if not binary.exists():
                    sys.exit(f"Binary not found: {binary}")
                if not prj.exists():
                    sys.exit(f"PRJ not found: {prj}")

            native_time  = run_ogs(native_bin,  native_prj,  native_outdir)
            mfront_time  = run_ogs(mfront_bin,  mfront_prj,  mfront_outdir)
        else:
            print("  --no-run: skipping simulation, using existing output.")

        results = compare_suite(suite, native_outdir, mfront_outdir)
        print_report(suite, results, native_time, mfront_time)


if __name__ == "__main__":
    main()
