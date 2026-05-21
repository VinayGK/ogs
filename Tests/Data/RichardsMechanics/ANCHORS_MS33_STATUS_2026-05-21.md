# ANCHORS MS33 status after EOS Pi-path fix (2026-05-21)

## Code changes completed
- `ProcessLib/RichardsMechanics/RichardsMechanicsFEM-impl.h`
  - `computeReferenceMicroPorositySwellingStressIncrement`: Pi-path now uses
    `rho_LR` (bulk liquid density) instead of `rho_lR` (micro-liquid EOS state).
  - `computeSwellingStressIncrement`: propagated `rho_LR` argument.
  - `updateSwellingState`: propagated `rho_LR` argument and call sites updated.
- `Tests/ProcessLib/RichardsMechanics/DSMMicroMacroSingleIntegrationPoint.cpp`
  - Updated for new function signature and Pi-path expectation.

## Recalibration (post-fix)
- `K_opt(dd1400) = 7656.5016 J/kg`
- `K_opt(dd1600) = 29999.2513 J/kg`
- `K_opt(dd1800) = 118585.8600 J/kg`
- Stored in `ANCHORS_MS33_ModelIV/ms33_calibrate_K_results.txt`.
- Propagated to baseline PRJs:
  - Model I dd1400/dd1600/dd1800
  - Model III gap2mm
  - Model IV pellets
  - Model VII freeswelling

## Verification results (spec-compliant rerun)
- Build: `cmake --build /Users/vinaykumar/git/build/release-omp-mfront -j8` passed.
- Model I Villar verification:
  - dd1400: `p_sw=1.50378 MPa`, err `-0.0020%`
  - dd1600: `p_sw=5.82306 MPa`, err `-0.0174%`
  - dd1800: `p_sw=22.56225 MPa`, err `+0.0278%`
  - mean absolute error: `0.0158%`
- Full MS33 suite:
  - I dd1400: `308 accepted / 0 rejected`
  - I dd1600: `308 accepted / 0 rejected`
  - I dd1800: `308 accepted / 0 rejected`
  - III gap2mm: `299 accepted / 0 rejected`
  - IV pellets: `279 accepted / 0 rejected`
  - VII freeswelling: `458 accepted / 0 rejected`
- Final-time swelling pressure sign check:
  - III: `5.493395 MPa` (>0)
  - IV: `4.294508 MPa` (>0)
  - VII: `4.746589 MPa` (>0)

## Mandatory spec-compliance policy (effective 2026-05-21)
- All ANCHORS MS33 simulations must remain strictly benchmark-spec compliant at all times.
- Any deviation from specified ICs, BCs, geometry, material parameters, loading path, or solver-relevant benchmark settings is not allowed in committed/run configurations.
- No workaround edits (for convergence, convenience, or sensitivity exploration) may overwrite the spec baseline inputs.
- Non-spec alternatives may only be proposed/discussed in chat as optional investigations and must be clearly labeled non-spec.
- If a non-spec test is explicitly requested, it must be isolated in clearly named scratch files and must not replace or mutate spec baseline files.
