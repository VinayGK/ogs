# AGENTS.md — ANCHORS MS33 DSM (streamlined)

## Scope
- Models: I / III / IV / VII (plus V_LE/V_MCC where present).
- Core implementation: `ProcessLib/RichardsMechanics/RichardsMechanicsFEM-impl.h`.
- Branch: `dsm_native_hierarchical`.
- Binary: `/Users/vinaykumar/git/build/release-omp-mfront/bin/ogs`.

## Roadmap (one-line commit refs)
- Step 1 REV-scale storage + split consistency: `0d7a9edd64`.
- Step 2 thermodynamic swelling stress + K recalibration: `88d42c98fd`.
- Step 3 Pi-path Gibbs–Duhem consistency + flag cleanup: `c4888b6db4`, `ce9178fa96`.
- Step 5 vdW dimensional fix (`/rho_lR`) + literature A lock: `0d579e8aeb`.
- Step 6 DSM hardening (viscosity guards, micro-pressure density default true): `66b782afa1`.
- Step 7 dead-code removal (compatibility overload/unused flag): `4d47efff55`, `ce9178fa96`.
- Step 8 DSM micro-macro test refactor (13/13 passing): `3ac6b7de1f`.

## Key physics/implementation invariants
- Porosity split: `phi = phi_M + phi_m`, with micro state carried by `n_l`.
- Storage is REV-scale: `phi_m * rho_lR`.
- Swelling stress uses thermodynamic Pi-path tied to `rho_LR` for Gibbs–Duhem consistency.
- vdW potential terms are additive; never replace additive update with assignment.
- `hamaker_constant` is literature/material-fixed; calibration target is K (not A).

## Execution instructions
- Keep committed runs benchmark-spec compliant.
- After physics changes, require:
  1. Model-I Villar target check within tolerance,
  2. canonical LE reruns with zero rejected steps,
  3. `./bin/testrunner --gtest_filter='*DSMMicroMacro*'` passing.

## Current summary
- Production path stable under latest DSM fixes.
- Canonical LE outcomes unchanged in accepted/rejected-step sense.
- Open benchmark-quality work is primarily calibration/interpretation side (not immediate solver-break state).
