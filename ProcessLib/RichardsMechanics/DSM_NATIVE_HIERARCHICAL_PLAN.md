# DSM Native Hierarchical Porosity Swap Plan

## Goal
Replace the additive/naive micro-macro porosity split in the RM DSM-native path with the hierarchical split used in the PINION derivation.

## Target model

- State variable kept as-is:
  - `n_l` = micro water content (`MicroWaterContent`), aggregate-referenced.
- Porosity split law (REV-referenced):
  - `phi = phi_M + (1 - phi_M) * n_l`
  - `phi_M = (phi - n_l) / (1 - n_l)`
  - `phi_m = (1 - phi_M) * n_l`

## Implementation scope

File:
- `ProcessLib/RichardsMechanics/RichardsMechanicsFEM-impl.h`

Changes:
1. Replace `computeTransportPorosityUpdate(...)` with hierarchical split evaluation.
2. Keep `macro_porosity_update_mode="additive_macro_porosity_rate_mode"` as a config alias, but evaluate the hierarchical law.
3. Use previous micro-porosity (`phi_m_prev`) consistently in local solve context (instead of previous `n_l`).
4. Initialize `n_l` from `(phi, phi_M)` via hierarchical inversion, not via `phi - phi_M`.
5. Initialize `phi_m` from hierarchical relation `phi_m = (1 - phi_M) * n_l`.
6. Update split state consistently for all DSM-local modes through `computeTransportPorosityUpdate(...)`.

## Expected behavior changes

- `micro_porosity` now represents REV micro porosity from the hierarchical split, not directly `n_l`.
- Macro/micro split remains bounded and consistent with total porosity.
- Legacy project files using `additive_macro_porosity_rate_mode` run on this branch with hierarchical split behavior.

## Validation checklist

1. Compile RM process.
2. Run DSM micromacro Richards tests, especially:
   - `beacon_1a01_dsm_micromacro_inflowprobe`
   - `beacon_1a01_dsm_micromacro_stressprobe`
3. Verify:
   - no negative porosity states,
   - `phi_m <= phi`,
   - stable local implicit `n_l` solve,
   - expected trend shifts vs additive baseline.
