# Agent instructions for OpenGeoSys-6

**Project**: Scientific THMC simulation framework in modern C++. See [README.md](./README.md) for overview.

## Mandatory code standards

**C++**: C++23 standard (required). Use `ranges-v3` (not `std::ranges`), Eigen (not raw loops). Follow [style guide](https://ufz.github.io/styleguide/cppguide.xml).

**CMake**: 3.31+, target-based commands, no global variables.

**Python**: black-based formatting. ruff check linter.

**Documentation**: British English spelling (colour, behaviour).

## Architecture layers

Foundation:  BaseLib → MathLib → NumLib
Geometry:    GeoLib, MeshLib, MeshGeoToolsLib, MeshToolsLib
Materials:   MaterialLib (MPL), ParameterLib
Processes:   ProcessLib (20+ process implementations)
Apps:        CLI, ApplicationsLib, FileIO, Utils, DataExplorer(Qt)

## Process implementation pattern (REQUIRED)

Every process must have:

1. `{Name}Process.h` - Inherits Process, manages assembly/timestepping
2. `{Name}ProcessData.h` - Material properties, parameters, solver configuration
3. `{Name}LocalAssembler.h` - Element-level assembly (M, K, b matrices)
4. `Create{Name}Process.h` - Factory from XML configuration

## Testing & validation

- Unit tests: `Tests/{LibName}/` (Google Test)
- Integration tests: `Tests/Data/{ProcessName}/` (`.prj` files with reference outputs)
- Always run ctests from release build.
- Check `.clang-format`, `.clang-tidy` for linting rules

## DSM native transition checkpoint

- Active native DSM source tree:
  `/Users/vinaykumar/git/ogs-native-dsm-transition`
- Verified build tree:
  `/Users/vinaykumar/git/build/release-native-transition2`
- As of 2026-05-06, the focused DSM native gate is green:
  `ctest -j 18 --output-on-failure -R "ogs-RichardsMechanics(/DoubleStructureBenchmark/double_porosity_swelling_RM|_.*dsm_micromacro)"`
  passed `32/32`.
- The failing DSM native checks were stale `vtkdiff` references, not solver
  failures and not a tolerance-only issue. Serial and OpenMP regenerated VTUs
  were byte-identical for the failing reference cases.
- The refreshed reference files are:
  - `Tests/Data/RichardsMechanics/beacon_1a01_reference_t_1000.000000.vtu`
  - `Tests/Data/RichardsMechanics/beacon_1a01_dsm_micromacro_inflow_reference_t_100000.000000.vtu`
  - `Tests/Data/RichardsMechanics/beacon_1c_reference_t_1000.000000.vtu`
- DSM source-audit guardrail: keep
  `ProcessLib/RichardsMechanics/ConstitutiveRelations/PotentialExchange.h`
  aligned with the MFront vdW regularisation: use
  `omega3 + omega_min_vdw3` with `h_min_vdw = 5e-11`, and keep the analytic
  derivatives consistent with that denominator.

## 2026-05-11 — Appendix A/B sync protocol (paper-code-presentation)

- Canonical derivation source is Appendix A in the paper:
  `/Users/vinaykumar/tex/dsm-bgr-paper/draft/paper_DSM.tex`
  section *Porosity split and reference-volume consistency*.
- Canonical comment-status source is Appendix B:
  `/Users/vinaykumar/tex/dsm-bgr-paper/draft/appendix_comment_status.tex`.
- Keep DSM-native implementation narrative aligned with Appendix A replacement logic:
  - macro porosity rate follows total-minus-micro coupling,
  - micro evolution is exchange-driven in the local residual,
  - exchange uses `\mu^Macro-\mu^Micro` at consistent REV scale.
- Do not reintroduce ad-hoc `A_H` scaling or ad-hoc `\left(\alpha_\text{Biot}-\phi\right)`
  split statements in docs/scripts when Appendix A equations are the active basis.
- If implementation changes affect Nagel comments #5--#14 scope (reference
  volume, prefactors, exchange scale, missing micro-density-rate term), update
  Appendix B status rows in the same work cycle.
- For Appendix B edits, preserve status taxonomy and validate:
  row count, per-status counts, and `N + M + K = total`.

## 2026-05-19 — augmented-force Villar calibration on hierarchical porosity

- Branch context: `dsm_native_hierarchical` (hierarchical porosity split active in RM DSM-native path).
- Repro command:
  `python3 Tests/Data/RichardsMechanics/ANCHORS_MS33_ModelI/run_villar_dense_dd_native_augmented_calibration.py --ogs-bin /Users/vinaykumar/git/build/release-omp-mfront/bin/ogs --lib-path /Users/vinaykumar/git/build/release-omp-mfront/lib`
- Result artefacts are written under
  `Tests/Data/RichardsMechanics/ANCHORS_MS33_ModelI/` with prefix
  `villar_dense_dd_native_augmented_lam1en06_*`.
- Run notes and outcomes are documented in:
  `Tests/Data/RichardsMechanics/ANCHORS_MS33_ModelI/VILLAR_AUGMENTED_HIERARCHICAL_RUN_2026-05-19.md`.

### Swelling formula fix (same date)

Pre-fix, `computeReferenceMicroPorositySwellingStressIncrement` used `delta_phi_m`
as the swelling driver. In the hierarchical split `phi_m = (1-phi_M)*n_l`, so at
typical density (phi_M ≈ phi0 ≈ 0.5) the effective slope was halved vs. the
calibrated value. Fix: changed driver to `delta_n_l` (the `MicroWaterContent` state
variable), consistent with the parameter name `micro_water_content_swelling_slope`
and backward-compatible with all pre-hierarchical slope calibrations.

Code location: `ProcessLib/RichardsMechanics/RichardsMechanicsFEM-impl.h`,
`computeReferenceMicroPorositySwellingStressIncrement`.

### Pressure-cap guardrail

The LinearElastic skeleton (E=52 MPa, ν=0.3 → K_bulk ≈ 43 MPa) caps swelling at
`P_sw_max = K_bulk × slope × phi0 ≈ 2 MPa` regardless of augmentation prefactor K.
All Villar targets above ≈1450 kg/m³ require pressures well above this cap.
High-density Villar matching requires a pressure-dependent constitutive model
(ModCamClay via MFront, as used in the separate mfront calibration script).
Treat residual large errors as a constitutive-model scope limit, not a solver failure.

## 2026-05-19 — MCC Villar calibration attempt and swelling-pressure density-dependence problem

### What was tried

Created `Tests/Data/RichardsMechanics/ANCHORS_MS33_ModelI/run_villar_dense_dd_native_augmented_mcc_calibration.py`
using `ModCamClay_semiExpl_constE` (E=52 MPa, ν=0.3, M=1.2, κ=6.6e-3, λ=7.7e-2,
p_c0=0.2 MPa) with the same augmented vdW exchange. Also tested `ModCamClay_semiExpl_absP`.

### Numerical stability boundary

MCC constE is stable only for augmentation K ≤ ~260 J/kg. For K > ~300 the MFront
semi-explicit integration fails with status -1 (rapid n_l exchange drives trial stress
far past p_c0 in a single day's step). The calibrated K for ρ_d=1400 is 964 J/kg —
already 3× the stability limit. Finer time stepping (6-hourly) and adaptive stepping
both also fail.

Stable MCC constE at K=250-290 gives P_sw ≈ 0.37 MPa for ρ_d=1400 via plastic
hardening (vs LinearElastic's 0.20 MPa at same K), but the Villar target is 1.5 MPa.

### Fundamental model problem: swelling pressure decreases with dry density

The swelling formula in `RichardsMechanicsFEM-impl.h`:
```cpp
delta_eps_sw = slope * delta_n_l;
delta_sigma_sw -= C_el * ((delta_eps_sw / 3.0) * identity2);
```
treats swelling as an **elastic eigenstrain** and generates
`P_sw_max = K_bulk × slope × phi0`. Because `phi0 = 1 − ρ_d/ρ_solid`
**decreases** with dry density, the cap **falls** as density rises:

| ρ_d (kg/m³) | phi0  | P_sw_max (MPa) | Villar target (MPa) |
|-------------|-------|----------------|---------------------|
| 1400        | 0.496 | 2.15           | 1.50  ✓             |
| 1475        | 0.469 | 2.03           | 2.50  ✗             |
| 1600        | 0.424 | 1.84           | 5.82  ✗             |
| 1800        | 0.353 | 1.53           | 22.6  ✗             |

**Root cause — wrong physical origin of swelling stress**: Swelling pressure originates
from the **disjoining pressure** Π — the reduction of combined micro-scale potentials
(μ_vdW, osmotic, hydration) as water content changes in the micro pores. It is of
**hydraulic origin**, not mechanical:

```
Π(ξ) = ρ_LR × K × exp(−ξ),  ξ = n_l / (λ × n_S × ρ_SR × Sa)
```

The correct macroscopic swelling stress is Π scaled by the **mineral surface area per REV**:

```
σ_sw = Π(ξ) × n_S × ρ_SR × Sa × f_geom
```

where `n_S = 1 − φ_M` is the aggregate fraction. Higher ρ_d → larger n_S → more
mineral surfaces per REV → also smaller ξ → stronger Π(ξ). Both effects raise P_sw
with density. Whether linear or exponential is material behaviour from Π(ξ).

The current eigenstrain formula has no n_S scaling, so the cap falls with density —
the opposite of physical reality.

### Suggested fix (not yet implemented)

**Primary fix — replace the swelling stress formula in
`computeReferenceMicroPorositySwellingStressIncrement`**:
Change from the elastic eigenstrain `K_bulk × slope × Δn_l` to the disjoining-pressure
coupling `Π(ξ) × n_S × ρ_SR × Sa × f_geom`. This removes the φ₀ cap entirely; P_sw
is then set by the exchange equilibrium `K × exp(−ξ_eq) = p_L_eq`, which naturally
grows with density through the n_S dependence of ξ_eq. No change to the constitutive
model is strictly required.

**Secondary concern — MFront numerical stability**: `semiExpl` integration fails
(status -1) for K > ~260 J/kg because rapid exchange drives the trial stress far past
p_c0 in a single day-step. This is an independent issue — a more implicit MFront
behaviour or an exchange-rate limiter is needed, regardless of which swelling-stress
formula is used.

**Documentation**: The problem and fix are described in the beamer presentation
`~/tex/cc2024/VK_B35_Pinion_May_2026/nagel_porosity_split.tex`, section
"Open problem: swelling-pressure density dependence" (last 3 slides before Questions).

### Agent implementation protocol (model correction)

- Scope: apply this protocol whenever editing
  `computeReferenceMicroPorositySwellingStressIncrement` or any RM DSM-native
  swelling-pressure term.
- Replace the swelling-stress driver from elastic eigenstrain form
  `K_bulk * slope * delta_n_l` to a disjoining-pressure-based coupling at REV scale:
  `sigma_sw = Pi(xi) * n_S * rho_SR * Sa * f_geom`.
- Keep exchange and swelling terms consistent in variables and scaling:
  use the same `n_l`, `n_S`, and REV conventions as in `PotentialExchange.h` and
  the hierarchical porosity split derivation.
- Verification requirement: include at least one density-trend check proving that
  predicted `P_sw` increases with `rho_d` for Villar points (1400 to 1800 kg/m^3).
- Numerical-stability guardrail: if MCC semi-explicit still fails for high K,
  treat that as a separate integration issue (implicit behaviour or exchange-rate
  limiting), not a reason to revert the physical swelling-stress correction.
- Commit hygiene: do not commit generated `.pvd/.vtu/.prj` artefacts from
  calibration sweeps unless explicitly requested for reference updates.
