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

## DSM_NATIVE_HIERARCHICAL_PATCH_RECIPE.md maintenance rule

`ProcessLib/RichardsMechanics/DSM/DSM_NATIVE_HIERARCHICAL_PATCH_RECIPE.md` is the reconstruction
recipe for this branch from a fresh `master`. It must stay current.

**Update DSM_NATIVE_HIERARCHICAL_PATCH_RECIPE.md before committing whenever:**
- Any hunk in `RichardsMechanicsFEM-impl.h` or `PotentialExchange.h` changes.
- The DSMMicroMacro unit tests change (step 8 section + passing count).
- Any PRJ `hamaker_constant` or `vdw_augmentation_prefactor` (K) value changes.
- A new benchmark model is added to the canonical LE set.
- Build flags or the verification `ctest` invocation changes.
- A new step beyond Step 8 is added (add a new numbered section).

Do not mark a step done in AGENTS.md unless DSM_NATIVE_HIERARCHICAL_PATCH_RECIPE.md already reflects it.

## Known limitations (logged 2026-05-27)

### Hydraulic-side double-counting of suction in Darcy flux

**Location:** `ProcessLib/RichardsMechanics/RichardsMechanicsFEM-impl.h:2627` (and
the parallel branches at 3368 and 3991) — assembly of the macro Darcy flux:

```cpp
laplace_p += dNdx_p^T · (k_intr·k_rel·ρ_LR/μ) · dNdx_p · w
// drives:  q_L = -(k_intr·k_rel/μ) · ∇p_L
```

**Issue:** `p_L` is the primary `pressure` process variable, which the boundary
conditions impose as *total* suction (lab-measured, capillary + molecular).
The full ∇p_L therefore drives the macro Darcy flux, including the molecular
(disjoining-pressure Π) component. Physically the molecular component should
drive the micro→macro mass exchange via the DSM source term ρ̇_micro→macro,
NOT bulk advection in the macro pore.

**Manifestation:** boundary suction ramps of order 100 MPa create huge ∇p_L
gradients that propagate the wetting front much faster than the material
physically would. Affects all transient III/IV/VII results at finite t. Does
NOT affect Model I (no spatial gradient) and does NOT affect t→∞ equilibria
of III/IV/VII (saturated swelling pressure, asymptotic gap closure, equilibrium
void ratio).

**Mechanical-side companion (already fixed):** the same total p_L was being
fed into Bishop's effective stress via the `BishopsPowerLaw(exponent=1)` path
(χ = S_L → χ·p_L = molecular component leaks into σ_eff during unsaturated
phase). Fixed by switching all 4 MS33 PRJs to `BishopsSaturationCutoff(cutoff=1)`
so χ=0 below S_L=1; DSM Π then carries all swelling-source work in the
unsaturated regime.

**Right fix (open):** split p_L into p_L_macro (capillary, bounded by
Young–Laplace ~3 MPa for compacted bentonite) and p_L_micro (disjoining, Π).
Use ∇p_L_macro only in Darcy; route the (p_L − p_L_macro) residual through
the existing DSM micro↔macro exchange. Requires a new constitutive law for
the split and a process-variable refactor in `RichardsMechanicsProcess`.
Estimate: 1–2 weeks careful work + tests.

**Pragmatic interim:** submit current results with explicit caveat in the
deliverable (deck frame "Known Limitation — Hydraulic Side of OGS RM-DSM").
