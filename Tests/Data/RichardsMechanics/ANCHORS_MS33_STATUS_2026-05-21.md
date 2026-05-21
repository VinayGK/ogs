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

## TODO — Model V salinity characterization (literature-based)
- Current Model V runs are saline surrogates; they are not yet validated against chemistry-informed HM response.
- Required next step: perform a focused literature review to quantify how salinity affects bentonite HM behavior (at minimum: swelling pressure reduction, permeability evolution, retention/capillary response, and stiffness/compressibility trends where supported).
- After extracting defensible ranges/trends, parameterize salinity effects explicitly in the benchmark model and re-run Model V (LE and MCC variants) with documented provenance for each changed parameter.
- Acceptance target: saline Model V must show physically justified deviation from freshwater reference (especially swelling pressure path), not only numerical convergence.

## Mandatory spec-compliance policy (effective 2026-05-21)
- All ANCHORS MS33 simulations must remain strictly benchmark-spec compliant at all times.
- Any deviation from specified ICs, BCs, geometry, material parameters, loading path, or solver-relevant benchmark settings is not allowed in committed/run configurations.
- No workaround edits (for convergence, convenience, or sensitivity exploration) may overwrite the spec baseline inputs.
- Non-spec alternatives may only be proposed/discussed in chat as optional investigations and must be clearly labeled non-spec.
- If a non-spec test is explicitly requested, it must be isolated in clearly named scratch files and must not replace or mutate spec baseline files.
