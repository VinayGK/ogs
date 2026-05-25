# RichardsMechanics Benchmark Script Inventory

Last updated: 2026-05-25 (CEST)  
Repository branch: `dsm_mfront_hierarchical`

This file documents and consolidates the scripts used for DSM/native/MFront
calibration and benchmark exercises (ANCHORS, BEACON, EPFL/BGR).

## Main End-to-End Workflows

| Script | Purpose | Primary outputs |
|---|---|---|
| `../../scripts/run_dsm_parity.py` | **DSM native-vs-mfront parity runner and MAE reporter.** Runs both OGS binaries side-by-side, compares every registered field, prints structured MAE table (early ts≤50 / late ts≥60). Suites registered in `PARITY_SUITES` at the bottom of the file — append a dict to add a new model pair without touching the runner logic. PRJ files live in `ANCHORS_MS33_StrictParity/`. | `/tmp/dsm_parity_runs/<suite>/{native,mfront}/` |
| `run_linear_law_calibrated_dsm_comparison.py` | ANCHORS linear-law calibration + triplet comparison (`linear`, `native_dsm`, `mfront_dsm`) for ANCHORS/BEACON/EPFL. | `_outputs/linear_law_calibrated_dsm_comparison/*` |
| `run_identical_parameter_benchmark_parity.py` | Same-parameter parity run for ANCHORS + BEACON + EPFL (`native` vs `mfront`) using shared dry-density multiplier curve. | `_outputs/identical_parameter_native_mfront_comparison_with_epfl/*` |
| `run_calibrated_beacon_anchors_comparison.py` | Dry-density-calibrated comparison for ANCHORS + BEACON (`native` vs `mfront`) with threshold assertions. | `_outputs/calibrated_native_mfront_comparison/*` |

## Calibration Scripts (ANCHORS)

| Script | Scope | Notes |
|---|---|---|
| `ANCHORS_MS33_ModelI/run_villar_dense_dd_calibration.py` | MFront DSM calibration across dry-density sweep. | Produces `villar_dense_dd_calibration.csv`. |
| `ANCHORS_MS33_ModelI/run_villar_dense_dd_native_dsm_micromacro_calibration.py` | Native DSM calibration across dry-density sweep. | Produces `villar_dense_dd_native_dsm_micromacro_calibration.csv`. |
| `ANCHORS_MS33_ModelI/run_villar_dense_dd_native_dsm_micromacro_branch.py` | Native DSM branch sweep with fixed settings (non-iterative branch response). | Useful for sanity checks vs target curves. |

## BEACON / EPFL Analysis and Plotting

| Script | Scope | Notes |
|---|---|---|
| `plot_calibrated_benchmark_comparison.py` | Plots calibrated ANCHORS/BEACON native-vs-MFront comparisons. | Reads calibrated delta CSVs. |
| `plot_beacon_report_figures.py` | BEACON report-style figure generation (`1a01`, `1b`, `1c`). | Produces report-panel PNGs. |
| `analyze_beacon_unstructured_batch.py` | Postprocessing for BEACON unstructured runs. | Metrics extraction helper. |
| `calibrate_beacon_1b_ddcurve_probe.py` | Dry-density / multiplier sensitivity probing for BEACON `1b`. | Diagnostic calibration probe. |
| `BEACON_WP3_BGR_EPFL/plot_bgr_wp3_epfl_report_figures.py` | EPFL/BGR report-style figure generation. | Uses native + MFront EPFL decks. |
| `BEACON_WP3_BGR_EPFL/analyze_bgr_wp3_abprime.py` | EPFL/BGR AB' quantitative analysis helper. | CSV/JSON support utilities. |

## Mesh/Geometry Helper Scripts

| Script | Scope |
|---|---|
| `generate_beacon_unstructured_meshes.py` | BEACON unstructured mesh generation. |
| `BEACON_WP3_BGR_EPFL/generate_bgr_epfl_mesh.py` | EPFL/BGR mesh and geometry generation. |

## Ancillary / Reference Utilities

| Script | Scope |
|---|---|
| `run_beacon_native_bridge_report_batch.py` | Legacy native/bridge BEACON batch driver. |
| `run_native_branch_dsm_micromacro_comparison.py` | Native-branch comparison helper. |
| `ANCHORS_MS33_ModelI/plot_ms33_vs_ciemat_swelling_pressure_dd.py` | MS33 vs CIEMAT/Villar pressure-density plotting utility. |
| `ANCHORS_MS33_ModelI/summarize_ms33_model_i.py` | ANCHORS Model I summary helper. |
| `Mathematica_DSM_Driver/run_dsm_driver_demo.wl` | Mathematica DSM constitutive demo driver. |
| `Mathematica_DSM_Driver/DsmMicromacroDriver.wl` | Mathematica DSM constitutive equations/driver. |

## Storage Convention

- Script sources are stored under:
  - `Tests/Data/RichardsMechanics/`
  - `Tests/Data/RichardsMechanics/ANCHORS_MS33_ModelI/`
  - `Tests/Data/RichardsMechanics/BEACON_WP3_BGR_EPFL/`
  - `Tests/Data/RichardsMechanics/Mathematica_DSM_Driver/`
- Generated run artifacts are stored under `_outputs/` and are intentionally git-ignored.
- If a new benchmark workflow is added, update this inventory in the same commit.
