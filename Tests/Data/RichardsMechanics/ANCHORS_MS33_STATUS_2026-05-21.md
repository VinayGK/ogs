# ANCHORS MS33 status after nS consistency change (2026-05-21)

## Code change
- `ProcessLib/RichardsMechanics/RichardsMechanicsFEM-impl.h`
  - `computeActiveMicroSolidVolumeFraction`: use `n_S = 1 - phi_M` (instead of total-solid `1 - phi`).
  - `computePreviousMicroSolidVolumeFraction`: use previous `phi_M` only.

## What works
- Model I (dd1400, dd1600, dd1800): converged, 308 accepted / 0 rejected each.
- Model III (gap2mm): converged, 284 accepted / 0 rejected.
  - Deformed-gap estimate from displacement at 200 d: ~1.689 mm (from 2.000 mm initial).
- Model IV (pellets): converged, 268 accepted / 0 rejected.
  - `transport_porosity` remains non-negative in output (`min = 0.0`).
- Targeted RM regression checks passed:
  - `RichardsMechanics_double_porosity_swelling_dsm_micromacro_constbc_reference`
  - `RichardsMechanics_beacon_1c_*`

## What does not work
- Model VII (freeswelling): fails at first step due to nonlinear divergence; timestep cutbacks exhausted.
- Failure mode appears hydraulic/micro-exchange dominated (pressure component diverges first), not MCC plasticity (suite currently uses `LinearElasticIsotropic`).

## Calibration updates (after rerun)
- `K_opt(dd1400) = 3624.0887 J/kg`
- `K_opt(dd1600) = 16787.1049 J/kg`
- `K_opt(dd1800) = 75877.7890 J/kg`
- Saved in `ANCHORS_MS33_ModelIV/ms33_calibrate_K_results.txt` and propagated to Model I PRJs.
