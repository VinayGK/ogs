# ANCHORS MS33 Model-I

## 2026-05-27 — Swap to DSM-native hierarchical formulation

The three project files (`ms33_model_i_dd1400.prj`, `ms33_model_i_dd1600.prj`,
`ms33_model_i_dd1800.prj`) have been replaced with the EURAD-2 MS33 spec-aligned
DSM-native hierarchical formulation promoted from
`tex/cc2024/VK_SB_EURAD_DSM/MS33_LE_RERUN/ANCHORS_MS33_ModelI/`.

What changed:
- Constitutive law: phenomenological `LinearElasticIsotropic` +
  `SaturationDependentSwelling` (with a single 2 MPa swelling pressure shared
  by all three dry densities) replaced by DSM-native hierarchical model with
  per-density vdW augmentation prefactors `K` calibrated against Villar
  swelling-pressure targets (1.504 / 5.824 / 22.556 MPa at
  1400 / 1600 / 1800 kg/m³).
- Mesh: spec-aligned single-element cell `../square_1e-2_quad_1e0.vtu`
  (1 cm × 1 cm) instead of the previous `square_1x1_quad_1e0.vtu` 1 m cell.
- Saturation BC: applied to both top and bottom faces so all 4 nodes track
  the suction ramp (proper suction-controlled single-element benchmark).
- Output prefixes: `ms33_modelI_dd*` (camelCase) — produced files therefore no
  longer collide with the previous `ms33_model_i_dd*` snapshot artefacts in
  this directory. The in-tree PRJ filenames remain `ms33_model_i_dd*.prj`
  (snake_case) so that pre-existing references (CTest, scripts) keep working.
- Output variables now include `intrinsic_permeability` and
  `relative_permeability` (sibling C++ change enables emitting these).

Why the swap:
- The previous in-tree PRJs used the same 2 MPa swelling pressure for all
  three dry densities, which is physically wrong — Villar measurements span
  ~1.5 to ~22.6 MPa across this density range. The promoted spec-aligned
  PRJs are the only versions reproducing those targets.

## Historical content (pre-2026-05-27)

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
- the dense dry-density native-dsm_micromacro-branch Villar-style sweep
  (`run_villar_dense_dd_native_dsm_micromacro_branch.py`,
  `villar_dense_dd_native_dsm_micromacro_branch.csv`,
  `villar_dense_dd_native_dsm_micromacro_branch_summary.json`,
  `villar_dense_dd_native_dsm_micromacro_branch_vs_villar.png`),
- the dense dry-density native-dsm_micromacro-branch effective-vdW calibration
  run (`run_villar_dense_dd_native_dsm_micromacro_calibration.py`,
  `villar_dense_dd_native_dsm_micromacro_calibration.csv`,
  `villar_dense_dd_native_dsm_micromacro_calibration_summary.json`,
  `villar_dense_dd_native_dsm_micromacro_calibration_comparison.png`,
  `villar_dense_dd_native_dsm_micromacro_vdw_multiplier.png`),
- direct calibrated native-vs-MFront parity summary
  (`villar_dense_dd_native_vs_mfront_calibrated_comparison.csv`,
  `villar_dense_dd_native_vs_mfront_calibrated_comparison_summary.json`).

Runtime-only outputs (`*.vtu`, `*.pvd`) are intentionally not versioned.

Run provenance (including git hashes) is tracked in:
- `Tests/Data/RichardsMechanics/test_run_table.csv`

Dense calibration run command (17 points from 1400 to 1800 kg/m³):

```bash
HOME=/tmp MPLCONFIGDIR=/tmp/mplcache XDG_CACHE_HOME=/tmp \
python3 run_villar_dense_dd_calibration.py --dd-step 25 --rel-tol 0.02
```

Native dsm_micromacro branch Villar-style sweep command (same density grid):

```bash
HOME=/tmp MPLCONFIGDIR=/tmp/mplcache XDG_CACHE_HOME=/tmp \
python3 run_villar_dense_dd_native_dsm_micromacro_branch.py --dd-step 25
```

Native dsm_micromacro branch dense calibration command
(pointwise effective-vdW multiplier per dry density):

```bash
HOME=/tmp MPLCONFIGDIR=/tmp/mplcache XDG_CACHE_HOME=/tmp \
python3 run_villar_dense_dd_native_dsm_micromacro_calibration.py \
  --native-ogs /Users/vinaykumar/git/build/release-native-beacon/bin/ogs \
  --native-source /Users/vinaykumar/Documents/GitHub/ogs \
  --dd-step 25 --rel-tol 0.02
```

Key result from the native dsm_micromacro branch baseline run (`d46e11ac00`):
- simulated swelling pressure stays around `0.050` to `0.062` MPa for
  `1400` to `1800` kg/m³, while Villar Eq.(7) rises from about `1.50` to
  `22.56` MPa.
- mean mismatch over the 17 dry-density points is about `-7.99` MPa
  (native below target).

Key result from the updated native dsm_micromacro branch calibration run:
- baseline mean relative error vs Villar: `98.69%`.
- calibrated mean relative error vs Villar: `1.26%`.
- calibrated max relative error vs Villar: `1.83%`.
- calibrated pressure follows Villar closely over `1400` to `1800` kg/m³.

Native vs MFront comparison notes (same dsm_micromacro-derived equations target):
- native dense run now uses `scalar_micro_macro_mass_storage_mode` with
  EOS-driven `rho_lR` carry-over and robust multiplier bracketing in log-space.
- calibrated native and calibrated MFront curves are now close on the dense
  grid (mean absolute pressure difference about `0.135 MPa`, max about
  `0.467 MPa`).
- strict parameter-level parity is still open: fitted `vdw_multiplier`
  remains branch-dependent (MFront/native ratio roughly
  `3.1e5` to `8.2e5`).
