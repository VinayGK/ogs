#!/usr/bin/env python3
"""Run the BEACON report benchmarks for native and bridge implementations.

The script executes one paired case, captures the final VTU files, and
delegates comparison metrics to `analyze_beacon_unstructured_batch`.
"""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

from analyze_beacon_unstructured_batch import summarise


CASE_CONFIGS = {
    "1a01": {
        "native_project": "beacon_1a01_vk_inflow_unstructured_batch.prj",
        "bridge_project": "beacon_1a01_vk_notebook_mcc_inflow_unstructured_batch.prj",
        "native_output_prefix": "beacon_1a01_vk_inflow_unstructured_batch",
        "bridge_output_prefix": "beacon_1a01_vk_notebook_mcc_inflow_unstructured_batch",
        "final_time": "100000.000000",
    },
    "1b": {
        "native_project": "beacon_1b_vk_unstructured_batch.prj",
        "bridge_project": "beacon_1b_vk_notebook_mcc_unstructured_batch.prj",
        "native_output_prefix": "beacon_1b_vk_unstructured_batch",
        "bridge_output_prefix": "beacon_1b_vk_notebook_mcc_unstructured_batch",
        "final_time": "43200000.000000",
    },
}


def run_ogs(executable: Path, project: Path, output_dir: Path, source_dir: Path) -> None:
    """Execute OGS with a project file and capture a useful failure report."""
    output_dir.mkdir(parents=True, exist_ok=True)
    result = subprocess.run(
        [str(executable), "-o", str(output_dir), str(project)],
        cwd=source_dir,
        text=True,
        capture_output=True,
        check=False,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"ogs failed for '{project.name}' using '{executable}'.\n"
            f"stdout:\n{result.stdout}\n"
            f"stderr:\n{result.stderr}"
        )


def expected_output_file(output_dir: Path, prefix: str, final_time: str) -> Path:
    """Resolve the VTU file that should exist at the end of a case run."""
    path = output_dir / f"{prefix}_t_{final_time}.vtu"
    if not path.exists():
        raise FileNotFoundError(f"Expected output file not found: {path}")
    return path


def main() -> None:
    """Run one BEACON case and emit the summary JSON."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", choices=sorted(CASE_CONFIGS), required=True)
    parser.add_argument("--native-ogs", type=Path, required=True)
    parser.add_argument("--bridge-ogs", type=Path, required=True)
    parser.add_argument("--source-dir", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--json-out", type=Path)
    args = parser.parse_args()

    config = CASE_CONFIGS[args.case]
    source_dir = args.source_dir.resolve()
    native_project = source_dir / config["native_project"]
    bridge_project = source_dir / config["bridge_project"]

    native_output_dir = (args.output_root / args.case / "native").resolve()
    bridge_output_dir = (args.output_root / args.case / "bridge").resolve()

    run_ogs(args.native_ogs.resolve(), native_project, native_output_dir, source_dir)
    run_ogs(args.bridge_ogs.resolve(), bridge_project, bridge_output_dir, source_dir)

    native_vtu = expected_output_file(
        native_output_dir, config["native_output_prefix"], config["final_time"]
    )
    bridge_vtu = expected_output_file(
        bridge_output_dir, config["bridge_output_prefix"], config["final_time"]
    )

    result = summarise(args.case, native_vtu, bridge_vtu)
    result["artifacts"] = {
        "native_project": str(native_project),
        "bridge_project": str(bridge_project),
        "native_vtu": str(native_vtu),
        "bridge_vtu": str(bridge_vtu),
    }
    result["executables"] = {
        "native_ogs": str(args.native_ogs.resolve()),
        "bridge_ogs": str(args.bridge_ogs.resolve()),
    }

    text = json.dumps(result, indent=2, sort_keys=True)
    if args.json_out:
        args.json_out.write_text(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
