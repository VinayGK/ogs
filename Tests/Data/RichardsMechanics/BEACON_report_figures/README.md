# BEACON Report Figures

This folder stores report-style native vs MFront comparison assets for:
- `1a01`
- `1b`
- `1c`

Tracked artefacts here are:
- project files (`*.prj`) used to run the report plots,
- plot summaries (`*.csv`, `*.json`),
- figure panels (`*.png`) used in the transition note.
- native-branch replay comparison outputs
  (`native_branch_dsm_micromacro_comparison_report.json`,
  `native_branch_dsm_micromacro_comparison_summary.csv`,
  `native_branch_dsm_micromacro_stress_compare.png`).

Runtime-only simulation outputs (`*.vtu`, `*.pvd`) are intentionally not
versioned in this folder.

Run provenance (including git hashes) is tracked in:
- `Tests/Data/RichardsMechanics/test_run_table.csv`

## Native DSMMicroMacro Branch Validation

The native dsm_micromacro implementation was replayed from branch hash
`d46e11ac00` (`/Users/vinaykumar/git/ogs-native-dsm-transition`) and
compared against the current MFront bridge branch hash `f0a453fb89`
(`ogs-TPM_Swelling_MCC_Coupled`).

Run command:

```bash
python3 Tests/Data/RichardsMechanics/run_native_branch_dsm_micromacro_comparison.py \
  --output-root /tmp/ogs_native_branch_validation \
  --json-out Tests/Data/RichardsMechanics/BEACON_report_figures/native_branch_dsm_micromacro_comparison_report.json \
  --csv-out Tests/Data/RichardsMechanics/BEACON_report_figures/native_branch_dsm_micromacro_comparison_summary.csv \
  --plot-out Tests/Data/RichardsMechanics/BEACON_report_figures/native_branch_dsm_micromacro_stress_compare.png
```

Summary (final-output comparison):

| case | native axial sigma [kPa] | MFront axial sigma [kPa] | native radial sigma [kPa] | MFront radial sigma [kPa] |
| --- | ---: | ---: | ---: | ---: |
| `1a01_smoke` | 0.000 | 0.000 | 0.000 | 0.000 |
| `1b_smoke` | 0.000 | 0.000 | 0.000 | 0.000 |
| `1c_smoke` | 0.000 | 0.000 | 0.000 | 0.000 |
| `1a01_inflow` | 3.563 | 607.172 | 3.562 | 603.935 |

Observed behavior:
- `1a01/1b/1c` smoke decks are stress-parity checks and pass with effectively
  zero stress mismatch.
- `1a01_inflow` is not parity-calibrated between branches and shows a large
  stress gap (`~604` kPa) and pressure-state mismatch.
