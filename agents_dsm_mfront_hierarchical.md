# DSM MFront Hierarchical Port Log

## Objective
Create `dsm_mfront_hierarchical` as an increment over `dsm_mfront` that reproduces the hierarchical DSM capability of `dsm_native_hierarchical` in an MFront-idiomatic way.

## Branch set
- `master`
- `dsm_native`
- `dsm_native_hierarchical`
- `dsm_mfront`
- `dsm_mfront_hierarchical`

## Core semantic requirement
Reproduce `dsm_native_hierarchical` functionality via MFront-side implementation in the same architectural style by which `dsm_mfront` implements `dsm_native`.

## Native patch-recipe usage statement
`DSM_NATIVE_HIERARCHICAL_PATCH_RECIPE.md` from `dsm_native_hierarchical` was used as analytical input for semantic deltas and sequencing, not replayed mechanically into MFront code.

## Initial repo state
- Active branch at start: `dsm_native_hierarchical`.
- `dsm_mfront_hierarchical` already existed with one scaffold commit.
- Unexpected untracked files were present; stashed before work.

## Delta mapping
- `dsm_native -> dsm_native_hierarchical`:
  - hierarchical porosity split (`phi = phi_M + (1-phi_M) n_l`)
  - REV-consistent micro mass storage (`phi_m * rho_lR`)
  - swelling and potential-exchange consistency refinements.
- `dsm_native -> dsm_mfront`:
  - same DSM semantics expressed primarily via `RichardsMechanicsDSMMicroMacroBridge_MCC.mfront` and bridge plumbing.
- Target `dsm_mfront -> dsm_mfront_hierarchical`:
  - implement hierarchical split and REV mass terms in the MFront bridge where prior assumptions used flat `phi_m == n_l` semantics.

## Chronological engineering log
1. Stash workspace noise to ensure clean baseline.
   - Command: `git stash push -u -m "codex-temp-pre-dsm_mfront_hierarchical-2026-05-25"`
   - Result: clean working tree.
2. Analyze branch deltas and native patch recipe.
   - Commands: `git diff`, `git show dsm_native_hierarchical:...PATCH_RECIPE.md`, targeted `rg`.
   - Result: mapped native hierarchical semantics to MFront bridge loci.
3. Port step: hierarchical REV micro-porosity mapping in MFront bridge.
   - File: `MaterialLib/SolidModels/MFront/RichardsMechanicsDSMMicroMacroBridge_MCC.mfront`
   - Changes:
     - Introduced `phi_m_from_n_l(...)` under fixed `phi_total_trial`.
     - Replaced previous-mass baseline from `n_l_prev * rho_lR_prev` to `phi_m_prev * rho_lR_prev`.
     - Replaced residual mass terms from `n_l * rho_lR` to `phi_m(n_l) * rho_lR`.
     - Replaced trial-state assignment `phi_m_trial_value = n_l_trial_value` with hierarchical mapping and consistent `phi_M_trial_value`.
   - Reversibility note: isolated in one file and one logical commit.

## Replay instructions from clean `master`
1. `git checkout master`
2. `git checkout dsm_mfront`
3. `git checkout -b dsm_mfront_hierarchical`
4. Replay patches:
   - `git format-patch dsm_mfront..dsm_mfront_hierarchical -o /tmp/dsm_mfront_hierarchical_patches`
   - `git am /tmp/dsm_mfront_hierarchical_patches/*.patch`

## Patch-series generation
- `git format-patch dsm_mfront..dsm_mfront_hierarchical -o /tmp/dsm_mfront_hierarchical_patches`

## Validation matrix
- Build: not yet executed in this pass.
- ctests: not yet executed in this pass.
- Benchmarks/models: not yet executed in this pass.
- Parity assessment: partial semantic port applied (micro-mass/porosity mapping); full parity pending broader run matrix.
- Known deviations/limitations:
  - Remaining hierarchical parity points may still exist in MFront bridge and/or RM plumbing.
  - Full benchmark and model replay pending.
