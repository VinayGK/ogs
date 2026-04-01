# ANCHORS MS33 Model-I

This folder contains a simplified MS33 Model-I confined swelling sweep at
three dry-density points.

Tracked artefacts include:
- project files (`ms33_model_i_dd1400.prj`, `ms33_model_i_dd1600.prj`,
  `ms33_model_i_dd1800.prj`),
- summarized outputs (`ms33_model_i_history.csv`,
  `ms33_model_i_summary.json`),
- the summarizer script (`summarize_ms33_model_i.py`),
- the CIEMAT/Villar overlay script and outputs
  (`plot_ms33_vs_ciemat_swelling_pressure_dd.py`,
  `ms33_vs_ciemat_swelling_pressure_vs_dry_density.csv`,
  `ms33_vs_ciemat_swelling_pressure_vs_dry_density.png`),
- the dense dry-density MFront-vs-native calibration artifacts
  (`run_villar_dense_dd_calibration.py`,
  `villar_dense_dd_calibration.csv`,
  `villar_dense_dd_calibration_summary.json`,
  `villar_dense_dd_swelling_pressure_comparison.png`,
  `villar_dense_dd_vdw_multiplier.png`),
- the dense dry-density native-notebook-branch Villar-style sweep
  (`run_villar_dense_dd_native_notebook_branch.py`,
  `villar_dense_dd_native_notebook_branch.csv`,
  `villar_dense_dd_native_notebook_branch_summary.json`,
  `villar_dense_dd_native_notebook_branch_vs_villar.png`),
- the dense dry-density native-notebook-branch effective-vdW calibration
  attempt (`run_villar_dense_dd_native_notebook_calibration.py`,
  `villar_dense_dd_native_notebook_calibration.csv`,
  `villar_dense_dd_native_notebook_calibration_summary.json`,
  `villar_dense_dd_native_notebook_calibration_comparison.png`,
  `villar_dense_dd_native_notebook_vdw_multiplier.png`).

Runtime-only outputs (`*.vtu`, `*.pvd`) are intentionally not versioned.

Run provenance (including git hashes) is tracked in:
- `Tests/Data/RichardsMechanics/test_run_table.csv`

Dense calibration run command (17 points from 1400 to 1800 kg/m³):

```bash
HOME=/tmp MPLCONFIGDIR=/tmp/mplcache XDG_CACHE_HOME=/tmp \
python3 run_villar_dense_dd_calibration.py --dd-step 25 --rel-tol 0.02
```

Native notebook branch Villar-style sweep command (same density grid):

```bash
HOME=/tmp MPLCONFIGDIR=/tmp/mplcache XDG_CACHE_HOME=/tmp \
python3 run_villar_dense_dd_native_notebook_branch.py --dd-step 25
```

Native notebook branch dense calibration attempt command
(pointwise effective-vdW multiplier per dry density):

```bash
HOME=/tmp MPLCONFIGDIR=/tmp/mplcache XDG_CACHE_HOME=/tmp \
python3 run_villar_dense_dd_native_notebook_calibration.py --dd-step 25 --rel-tol 0.02
```

Key result from the native notebook branch baseline run (`d46e11ac00`):
- simulated swelling pressure stays around `0.050` to `0.062` MPa for
  `1400` to `1800` kg/m³, while Villar Eq.(7) rises from about `1.50` to
  `22.56` MPa.
- mean mismatch over the 17 dry-density points is about `-7.99` MPa
  (native below target).

Key result from the native notebook branch calibration attempt:
- baseline mean relative error vs Villar: `98.74%`.
- calibrated mean relative error vs Villar: `62.02%`.
- calibrated max relative error vs Villar: `94.59%`.
- despite large fitted multipliers, final pressure plateaus near
  `1.2` to `1.8` MPa on this setup.

Native vs MFront comparison notes (same notebook-derived equations target):
- both paths implement the same reduced notebook mass-exchange core
  (`n_l` update from `mu_LR - mu_lR` and swelling from `Delta phi_m`),
  but they currently differ in practical constraints used during this
  calibration workflow.
- native notebook storage enforces `n_l <= phi` via porosity split update,
  while the current MFront scalar notebook-storage path does not apply an
  equivalent hard upper clamp in the Newton iterate path.
- calibration scripts also differ in initialization policy:
  MFront dense calibration keeps `n_l0` fixed per dry density from the
  literature Hamaker baseline, while the native calibration script recomputes
  `n_l0` from the effective Hamaker value per multiplier trial.
