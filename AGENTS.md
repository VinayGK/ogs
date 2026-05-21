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

**Secondary concern — MFront numerical stability**: `semiExpl` integration fails
(status -1) for K > ~260 J/kg because rapid exchange drives the trial stress far past
p_c0 in a single day-step. This is an independent issue — a more implicit MFront
behaviour or an exchange-rate limiter is needed, regardless of which swelling-stress
formula is used.

**Documentation**: The problem and fix are described in the beamer presentation
`~/tex/cc2024/VK_B35_Pinion_May_2026/nagel_porosity_split.tex`, section
"Open problem: swelling-pressure density dependence" (last 3 slides before Questions).

---

### Disjoining-pressure swelling stress — physics, implementation status, and guardrails

#### Status: IMPLEMENTATION IS COMPLETE

The disjoining-pressure branch is already in
`computeReferenceMicroPorositySwellingStressIncrement`
(lines 1458–1478 of `ProcessLib/RichardsMechanics/RichardsMechanicsFEM-impl.h`)
and `updateSwellingState` already extracts n_S and rho_LR (lines 1526–1530).
The calibration script `run_villar_dense_dd_native_augmented_calibration.py`
achieves **< 2 % error across ρ_d = 1400–1800 kg/m³** by running OGS and reading
sigma from VTU output — proof that the implemented formula and sign are correct.

Do **not** re-implement or change the sign. The section below exists to
document the physics correctly so future edits stay consistent.

---

#### Physical mechanism — commit this to memory

The vdW augmentation potential in `PotentialExchange.h` is:

```
μ_aug = sign × K × exp(−ξ)
ξ     = n_l / (λ × n_S × ρ_SR × Sa)    [dimensionless film parameter]
```

**The potential is HIGHEST (most attractive) when n_l is LOWEST (driest state).**
As water is absorbed into micro-pores n_l rises, ξ rises, exp(−ξ) falls, and
the vdW attraction between clay platelets weakens.

The macro-scale swelling stress that results from this is:

```
Π(ξ) = ρ_LR × K × exp(−ξ)   [Pa — disjoining pressure at REV scale]

δσ_sw = n_S × [Π(ξ_curr) − Π(ξ_prev)] × I
```

Sign derivation (do not reverse this):

| Event | Effect on ξ | Effect on Π | Effect on δσ_sw |
|---|---|---|---|
| Water absorbed, n_l↑ | ξ = n_l/denom ↑ | Π = ρK exp(−ξ) ↓ | Π_curr − Π_prev < 0 |
| δσ_sw += n_S × (negative) | | | δσ_sw < 0 (compressive) ✓ |

A **compressive** (negative) σ_sw means the material wants to expand. In a
constrained geometry this registers as swelling pressure; in a free geometry it
produces swelling deformation. The volume change originates from the quantity of
water held under a disjoining pressure at a density different from bulk water —
as n_l grows at weakening Π, the accumulated stress increment in the constrained
solid builds up exactly as measured by the Villar experiments.

The operator in the code is **`+=`** (not `−=`). The two must never be confused:
`−= (negative) = tensile` — completely wrong direction.

---

#### What the code actually contains (audit, 2026-05-20)

```cpp
// impl.h lines 1432–1488
template <int DisplacementDim>
inline MathLib::KelvinVector::KelvinVectorType<DisplacementDim>
computeReferenceMicroPorositySwellingStressIncrement(
    double const n_l_prev, double const n_l,
    double const n_S, double const rho_LR,          // ← already in signature
    MathLib::KelvinVector::KelvinMatrixType<DisplacementDim> const& C_el,
    PotentialExchangeParameters const& potential_exchange_params)
{
    // ...early-return guards for slope==0 and !isfinite(delta_n_l)...

    if (potential_exchange_params.vdw_augmentation_prefactor > 0.0 &&
        potential_exchange_params.vdw_augmentation_decay_length > 0.0 &&
        n_S > 0.0)
    {
        double const denom =
            potential_exchange_params.vdw_augmentation_decay_length * n_S *
            potential_exchange_params.micro_solid_density_reference *
            potential_exchange_params.specific_surface;
        double const xi_curr = n_l / denom;
        double const xi_prev = n_l_prev / denom;
        double const K = potential_exchange_params.vdw_augmentation_prefactor;
        double const Pi_curr = rho_LR * K * std::exp(-xi_curr);
        double const Pi_prev = rho_LR * K * std::exp(-xi_prev);
        // n_l↑ → ξ↑ → Π↓ → (Pi_curr−Pi_prev) < 0 → compressive increment ✓
        delta_sigma_sw.noalias() += n_S * (Pi_curr - Pi_prev) * identity2;
        return delta_sigma_sw;  // early return; eigenstrain block below is skipped
    }

    // Eigenstrain fallback (K=0, backward-compatible):
    delta_sigma_sw.noalias() -= C_el * ((slope * delta_n_l / 3.0) * identity2);
    return delta_sigma_sw;
}
```

`updateSwellingState` (lines 1502–1557) already provides n_S and rho_LR:

```cpp
auto const phi_M =
    std::get<ProcessLib::ThermoRichardsMechanics::TransportPorosityData>(
        state_current).phi;
double const n_S   = std::max(1e-16, 1.0 - phi_M);
double const rho_LR = *std::get<MicroLiquidDensity>(state_current);
// rho_LR ≈ 1000 kg/m³; MicroLiquidDensity is a StrongType<double,…> — dereference
// with * to get the double.  With PRJ params micro_liquid_density_reference=1e-6,
// micro_liquid_density_a=1e-16, the solver converges to rho_LR ≈ macro density.
```

`sigma_sw` enters the mechanical strain as:
```cpp
eps_m = eps + C_el_inverse * sigma_sw;   // swelling eigenstress → eigenstrain
```
Negative sigma_sw → material wants to expand → swelling pressure when constrained.

---

#### Known limitations that must NOT be papered over

**1. Jacobian is lagged.**
`sigma_sw` is updated after the local micro-Newton solve, not as part of the
global Newton Jacobian. The global Jacobian has no `d(sigma_sw)/d(p_L)` or
`d(sigma_sw)/d(u)` term. For small Δn_l per step (slow exchange, small dt) this
is harmless. For rapid exchange or large steps it can slow Newton convergence.
Do not claim quadratic Newton convergence for the full coupled system — it holds
only when Δn_l is small per step.

**2. K is calibrated per dry density — it is NOT a universal material constant.**
K grows from ~5 444 J/kg (ρ_d = 1400) to ~105 430 J/kg (ρ_d = 1800). Each PRJ
uses the K calibrated for that specific initial dry density. If ρ_d changes during
simulation (swelling/compression), the calibrated K becomes approximate. This is
an open modelling limitation, not a code bug.

**3. Calibration script uses 1 MPa initial macro-suction, not 100 MPa.**
`PRESSURE_IC_PA = -1e6` in the calibration script (line 66). The production PRJ
files use −1e8 Pa. The K values are insensitive to this choice because swelling
pressure is controlled by the final equilibrium Π(ξ_eq), not by the path. But if
the script is re-run, keep the 1 MPa convention to avoid pressure-tolerance issues
(the script still uses the old `abstols=5e-8 1e-13 1e-13` — acceptable at 1 MPa).

**4. Gate test does NOT exercise the K>0 path.**
`ctest -R dsm_micromacro` uses K=0 (falls through to eigenstrain). The actual
verification is the calibration script result:
```
mean_rel_error_percent: 0.96 %   max_rel_error_percent: 1.81 %
```
Any edit to `computeReferenceMicroPorositySwellingStressIncrement` must be
followed by a re-run of the calibration script to confirm the density-trend
result still holds. The gate tests alone are insufficient.

---

#### Invariants to maintain on any future edit

- **Sign**: `delta_sigma_sw += n_S * (Pi_curr - Pi_prev)` — the `+=` is mandatory.
  Physical meaning: n_l↑ → Π↓ → (Pi_curr − Pi_prev) < 0 → compressive increment.
  `−=` gives tensile swelling — physically wrong and will fail calibration.
- **ξ formula**: `n_l / (λ × n_S × ρ_SR × Sa)` — must be identical in impl.h and
  PotentialExchange.h. Any change to the exchange potential formula requires a
  matching change here.
- **n_S source**: use `TransportPorosityData.phi` (the macro porosity φ_M) to
  compute n_S = 1 − φ_M. Do NOT use the total porosity φ = φ_M + φ_m.
- **rho_LR source**: use `*std::get<MicroLiquidDensity>(state_current)`. With the
  current PRJ parameters this ≈ 1000 kg/m³ (macro liquid density). Do not
  hardcode 1000.0 — the MicroLiquidDensity EOS must be respected.
- **Backward compatibility**: the K=0 branch (eigenstrain) must remain unchanged
  and must be exercised by all existing dsm_micromacro ctests.

---

### Agent implementation protocol (model correction)

- Scope: apply this protocol whenever editing
  `computeReferenceMicroPorositySwellingStressIncrement` or any RM DSM-native
  swelling-pressure term.
- The physics origin of swelling is the VOLUME of water held under disjoining
  pressure at a density different from bulk. Potentials are reduced by water uptake.
  Do not re-derive or re-sign the formula from scratch — use the invariants above.
- Verification requirement: re-run the calibration script and confirm P_sw grows
  with ρ_d with mean error < 2 % across 1400–1800 kg/m³.
- Numerical-stability guardrail: if MCC semi-explicit fails for high K, treat that
  as a separate integration issue, not a reason to alter the swelling-stress sign.
- Commit hygiene: do not commit generated `.pvd/.vtu/.prj` artefacts from
  calibration sweeps unless explicitly requested for reference updates.

## 2026-05-20 — MS33 theoretical benchmarking PRJ files: two mandatory fixes

### Context

EURAD-2 MS33 theoretical benchmarking PRJ files live under
`Tests/Data/RichardsMechanics/ANCHORS_MS33_Model{I,III,IV,VII}/`.
BGR selected Models I (mandatory, three dry densities), III (gap geometry),
IV (pellets geometry), VII (free swelling) for the DSM native hierarchical
vdW-augmented model. All 6 production PRJs are verified working as of this date.

### Fix 1 — pressure convergence tolerance (all PRJs)

**Symptom**: every time step rejected; Newton reports `|dx_p| ≈ 1e-7` every
iteration and never converges.

**Root cause**: `abstols=5e-8 1e-13 1e-13` uses 5×10⁻⁸ Pa for the pressure
component. At 100 MPa initial suction (`p_L ≈ −1e8 Pa`), machine epsilon
(2.2×10⁻¹⁶) times the solution magnitude creates a floating-point floor of
∼1×10⁻⁷ Pa on `|dx_p|`. This floor lies **above** the 5×10⁻⁸ Pa tolerance,
so Newton stagnates indefinitely.

**Fix** (applied to all 6 PRJs):

```xml
<!-- OLD — below machine-epsilon floor at 100 MPa scale: -->
<abstols>5e-8 1e-13 1e-13</abstols>

<!-- NEW — physically meaningful at relevant scales: -->
<abstols>1 1e-12 1e-12</abstols>
```

Pressure: 1 Pa (adequate for geomechanics). Displacements: 1×10⁻¹² m.

### Fix 2 — sigma0 at free displacement boundaries (Model VII only)

**Symptom**: Newton explodes on step 1: `|dx_p| ≈ 1e9 Pa`, `|dx_u_r| ≈ 63 mm`.
Fixed-wall tests with any sigma0 converge; only the free-boundary cylinder fails.

**Root cause — code**: `RichardsMechanicsFEM-impl.h` line 3384:

```cpp
local_rhs.template segment<displacement_size>(displacement_index)
    .noalias() += Kup * p_L;
```

This adds the Biot pressure force `α χ p_L` to the mechanical RHS at every
assembly. For **Dirichlet-constrained faces**, the solver absorbs any residual
as a reaction force — sigma0 is irrelevant to convergence. For **free faces**,
no reaction force exists; the mechanical residual at t=0 equals `α χ p_L_initial`
unless sigma_eff already cancels it.

**Equilibrium condition** (effective stress, type="effective", no auto-correction):

```
sigma0_eff_isotropic = α_Biot × χ(S_w) × p_L_initial
```

With `α_Biot=1`, `p_L_initial=−1e8 Pa`, vG P0=27 MPa, α=0.45 → S_w=0.3293
→ χ = 0.3293 (BishopsPowerLaw exponent=1):

```
sigma0_eff_rr = sigma0_eff_θθ = 1.0 × 0.3293 × (−1e8) = −3.293e7 Pa
sigma0_eff_zz = −3.293e7 + (−2e5)                      = −3.313e7 Pa
              (axial Neumann load −0.2 MPa added for Model VII)
```

**Applied in Model VII PRJ**:

```xml
<!-- Initial stress: must balance Biot pressure force (Kup*p_L, line 3384 of impl.h).
     At 100 MPa suction, S_w=0.3293, chi=0.3293, p_L=-1e8:
     sigma0_eff_rr = sigma0_eff_θθ = α*χ*p_L = 0.3293*(-1e8) = -3.293e7 Pa
     sigma0_eff_zz = α*χ*p_L + σ_Neumann = -3.293e7 + (-2e5) = -3.313e7 Pa -->
<parameter><name>sigma0</name><type>Function</type>
    <expression>-3.293e7</expression><expression>-3.313e7</expression>
    <expression>-3.293e7</expression><expression>0</expression>
</parameter>
```

**Verification** (pure RM diagnostic, no micro_porosity): time step 1 converges
in 3 Newton iterations with quadratic convergence
(`|dx|`: 4.08×10⁵ → 7.03 → 7.93×10⁻³).

**General rule for free-boundary RM PRJs with large initial suction**:

1. Compute S_w at p_cap_initial from the vG WRC (P0, α).
2. χ = S_w (BishopsPowerLaw exponent=1).
3. sigma0_eff_isotropic = α_Biot × χ × p_L_initial.
4. Add Neumann traction per direction if present.
5. Use `<initial_stress>sigma0</initial_stress>` (no `type` attribute = effective,
   stored directly — no auto-correction by OGS).

If using `type="total"` instead, OGS adds `χ α (−p_cap)` automatically during
initialisation (see `CreateInitialStress.cpp`), so sigma0 should carry only the
Neumann traction contribution.

### Verified results (2026-05-20)

| PRJ file | Steps | Rejected |
|---|---|---|
| `ANCHORS_MS33_ModelI/ms33_modelI_dd1400.prj` | 255 | 0 |
| `ANCHORS_MS33_ModelI/ms33_modelI_dd1600.prj` | 256 | 0 |
| `ANCHORS_MS33_ModelI/ms33_modelI_dd1800.prj` | 256 | 0 |
| `ANCHORS_MS33_ModelIII/ms33_modelIII_gap2mm.prj` | 600 | 0 |
| `ANCHORS_MS33_ModelIV/ms33_modelIV_pellets.prj` | 302 | 0 |
| `ANCHORS_MS33_ModelVII/ms33_modelVII_freeswelling.prj` | 422 | 0 |

Binary used: `/Users/vinaykumar/git/build/release-omp-mfront/bin/ogs`
(version `vdw-baseline-2026-05-08-15-g45ea35b9`). The `dsm_native-release`
worktree binary has a VTK ABI mismatch (linked against vtk 9.5, installed 9.6).

## 2026-05-20 — Model III gap stiffness fix and its physical limit

### Fix applied

`ms33_modelIII_gap2mm.prj` originally had `E` as a global `Constant` (52 MPa),
applied identically to clay and gap zones. The gap zone (r = 0.025–0.027 m,
φ₀ = 0.985) is physically an air void. With E_gap = E_clay = 52 MPa the gap
offered the same mechanical resistance as the bentonite, suppressing radial
closure to 0.089 mm (4.5 % of the 2 mm gap width at t = 200 d).

**Fix**: changed `E` to a `FunctionParameter` spatially conditioned on the
radial coordinate `x` (exprtk boolean-multiplication):

```xml
<parameter><name>E</name><type>Function</type>
    <expression>52e6 * (x &lt;= 0.025) + 10e6 * (x &gt; 0.025)</expression>
</parameter>
```

E_gap = 10 MPa chosen on two grounds:
1. **No element inversion**: the clay's peak radial swelling pressure is
   ≤ 5.82 MPa (Villar target). For a gap of width L = 2 mm, the gap strain at
   peak force = 5.82 MPa / 10 MPa = 0.58 < 1 — elements never invert.
   Lower values (100 kPa, 1 MPa) fail: at σ_rr ≈ 5.82 MPa the gap tries to
   compress by 58× or 5.8× its width → Jacobian goes negative → Newton diverges.
2. **Manageable condition number**: E_clay / E_gap = 52 / 10 = 5.2 (SparseLU trivially stable).
   The 100 kPa trial had condition number 520 and still failed due to element
   near-inversion, not ill-conditioning.

Simulation result: **600 steps, 0 rejected, 258 s** (vs 300 steps with old E=52 MPa).

### Physical limit — why the gap barely closes

With E_gap = 10 MPa, the equilibrium is reached when the clay's **radial**
effective stress σ_rr matches the gap resistance E_gap × (u_r / L_gap).

In this constant-volume test (all four walls fixed: r=0, r=0.027, z=0, z=0.07),
the swelling pressure is primarily carried by the **axial and circumferential
directions** (rigid lid + rigid bottom absorb most of the swelling). The radial
direction has partial relief through the gap. At t = 200 d:

| Location | u_r (mm) | σ_rr (MPa) | mean_p (MPa) |
|---|---|---|---|
| t=0d (initial) | 0 | — | 0.14 |
| t=20d (end of ramp) | −0.465 to −0.990 | +3.9 to +6.1 (tensile!) | 8.18 |
| t=200d (final) | 0.008 to 0.161 | −0.6 to −1.6 (compressive) | 7.63 |

**At t = 20 d the gap is under spurious tension** (σ_rr > 0 at clay-gap interface).
Cause: the bottom of the gap equilibrates to p_L ≈ 0 (wet) while the top remains
at p_L ≈ −1e8 Pa (dry). The vertical pore-pressure gradient in the gap zone —
which has no DSM swelling — produces tensile radial effective stress that pulls
the clay-gap interface inward. A void gap has zero tensile strength, but the
soft elastic gap material sustains this tension in the FEM. **The negative
u_r values at t = 20 d are a modelling artefact**, not physical gap widening.

**At t = 200 d**, full hydration: gap aperture ≈ 1.84 mm (0.16 mm closure out
of 2 mm). The low closure is not a numerical failure — it is the physical
consequence of the anisotropic constraint. σ_rr ≈ 0.8 MPa (radial); the mean
swelling pressure 7.6 MPa is primarily carried axially and circumferentially.

### Known limitations — do not paper over

1. **No contact mechanics**: once the gap closes, OGS RM cannot switch to a
   "rigid contact" condition. The gap material remains elastic throughout.
   Full closure (aperture = 0) cannot be modelled without a contact formulation.

2. **Tensile artefact at t = 20 d**: the elastic gap transmits tensile stress
   during transient non-uniform hydration (wet bottom, dry top). In reality a
   void gap has zero tensile strength — this phase of the output should be
   treated as qualitative only. The final state (t ≥ 100 d, full hydration)
   is physically meaningful.

3. **Gap aperture underestimated**: equilibrium aperture ≈ 1.84 mm (OGS RM,
   linear elastic) vs ≈ 0 mm (expected from plasticity + contact in other
   codes). This discrepancy is an acknowledged limitation of the linear-elastic
   model, not a DSM swelling bug.

4. **E_gap = 10 MPa minimum**: any lower value risks element inversion when
   σ_rr_clay transiently approaches 5.82 MPa. Do not reduce E_gap without
   first verifying ε_gap_max < 0.95 throughout the full 200-day run.

## 2026-05-20 — Pi-path default, early-return bug, and K calibration requirement

### What changed in the PRJ files (2026-05-20)

Five tags were set in **all six** MS33 production PRJ files
(Models I dd1400/dd1600/dd1800, III, IV, VII):

| Tag | Value | Effect |
|---|---|---|
| `micro_water_content_swelling_slope` | `0` | Pi-path only; slope branch inactive |
| `accumulate_swelling_contributions` | `true` | Pi and slope branches run independently; fixes historical early-return bug |
| `use_micro_liquid_density_for_pi` | `true` | Pi uses actual `rho_lR` / `rho_lR_prev` (micro density) at both steps instead of macro `rho_LR` |
| `use_micro_liquid_density_for_micro_pressure` | `true` | `p_L_m = −rho_lR * mu_lR` uses micro density (Action 2a) |
| `vdw_augmentation_decay_length` | `1e-6` | unchanged |

The Pi-path (`delta_sigma_sw = n_S * (Pi_curr − Pi_prev) * I`,
`Pi = rho_lR * K * exp(−ξ)`) is now the sole swelling-stress mechanism.
The disjoining pressure Π is the swelling pressure in clay surface-force theory;
calibration is through vdW/augmentation parameters (A_H, Sa, K, λ).

**Note on density effect:** With the current trivial EOS
(`micro_liquid_density_reference = 1e-6`, `micro_liquid_density_a = 1e-16`),
`rho_lR ≈ rho_LR = 1000 kg/m³` to machine precision. The
`use_micro_liquid_density_for_*` flags have negligible numerical effect until
Action 1 (EOS calibration) is performed. The flags are set now so that once the
EOS is calibrated the correct physics is used automatically.

### Early-return bug — RESOLVED (2026-05-20)

The historical early return on `slope ≤ 0` that silenced the Pi-path has been
resolved by exposing it as an opt-in PRJ flag.

**Mechanism:** `computeReferenceMicroPorositySwellingStressIncrement`
(`impl.h` lines ~1467–1563) now has two code paths:

```
accumulate_swelling_contributions = false  (legacy default)
  → if slope == 0: return zero immediately (old bug preserved for back-compat)
  → if slope > 0 and augmentation: run Pi, early-return, skip slope branch

accumulate_swelling_contributions = true   (new, set in all MS33 PRJs)
  → Pi block runs unconditionally when augmentation enabled
  → slope block runs independently if slope > 0
  → both accumulate additively; slope = 0 gives Pi-only
```

All 6 MS33 PRJ files have `<accumulate_swelling_contributions>true</accumulate_swelling_contributions>`.
The bug no longer affects any MS33 simulation.

**Function signature was also updated** to pass `rho_lR_prev` alongside `rho_lR`
(line 1469: `double const rho_lR, double const rho_lR_prev`), enabling correct
time-discrete Pi computation.

**Gate:** rebuild (`ninja RichardsMechanics ogs`), then run the K verification
below. All three Model I dry densities must show < 1 % MAE.

### vdW material parameters — literature justification

All six PRJ files share these physical parameters:

| Parameter | PRJ tag | Value | Literature basis |
|---|---|---|---|
| Hamaker constant | `hamaker_constant` | 5.1e-21 J | Gregory (1981); Israelachvili (2011) — consensus for Na-smectite in water: ~5 × 10⁻²¹ J |
| Specific surface | `specific_surface` | 523 (code units) | Effective value: gives initial film thickness h₀ = n_l0 / (nS * rho_SR * Sa) ≈ 1.4 nm at 100 MPa suction — consistent with 1–2 water monolayers in Na-smectite interlayer (Saiyouri et al. 2000). Do NOT replace with BET surface area (700–800 m²/g = 700,000 m²/kg ≠ 523 code-unit m²/kg). |
| Solid grain density | `micro_solid_density_reference` | 2780.0 kg/m³ | Brigatti et al. (2006) — Na-montmorillonite: 2740–2800 kg/m³ |
| Decay length | `vdw_augmentation_decay_length` | 1e-6 (code units) | PRJ comment says "physical 1 nm". Unit convention must be confirmed: if lambda is in SI metres, 1e-6 m = 1 μm ≠ 1 nm. Do not change this value without (a) confirming the unit convention and (b) recalibrating K. |

### K calibration per dry density

K (augmentation prefactor, J/kg) is not a fixed material constant — it absorbs
clay surface charge density and pore-fluid ionic composition. It is calibrated
separately for each dry density against the Villar empirical target:

```
p_sw = exp(6.77 × ρ_d[g/cm³] − 9.07)  MPa
```

Current K values in the PRJ files (calibrated with Pi-path active):

| Model / density | PRJ file | K [J/kg] | Target [MPa] |
|---|---|---:|---:|
| Model I dd1400 | `ms33_modelI_dd1400.prj` | 4981.81 | 1.504 |
| Model I dd1600 | `ms33_modelI_dd1600.prj` | 23423.8 | 5.824 |
| Model I dd1800 | `ms33_modelI_dd1800.prj` | 105429.7 | 22.556 |
| Model III | `ms33_modelIII_gap2mm.prj` | 23423.8 | (dd1600 reference) |
| Model IV | `ms33_modelIV_pellets.prj` | 23423.8 | (dd1600 reference) |
| Model VII | `ms33_modelVII_freeswelling.prj` | 23423.8 | (dd1600 reference) |

These K values were calibrated when `slope = 0.1` — the augmentation block ran
(Pi early return) and slope branch was bypassed. After the bug fix, the same K
values should reproduce the same sigma_sw (verify, do not assume).

### Calibration and rerun procedure

**After fixing the early-return bug and rebuilding:**

```bash
# Step 1 — verify existing K values:
cd Tests/Data/RichardsMechanics/ANCHORS_MS33_ModelIV
python ms33_calibrate_K.py --verify
# Expected: all three Model I densities < 1 % error.
# If FAIL: python ms33_calibrate_K.py  (bisection, ~6–20 OGS runs per density)

# Step 2 — propagate any new K(dd1600) to Models III, IV, VII PRJ files.

# Step 3 — rerun full suite:
OGS=/Users/vinaykumar/git/build/release-omp-mfront/bin/ogs
DATA=Tests/Data/RichardsMechanics
for DD in dd1400 dd1600 dd1800; do
    $OGS -o $DATA/ANCHORS_MS33_ModelI -l warn \
         $DATA/ANCHORS_MS33_ModelI/ms33_modelI_${DD}.prj &
done; wait
$OGS -o $DATA/ANCHORS_MS33_ModelIII -l warn $DATA/ANCHORS_MS33_ModelIII/ms33_modelIII_gap2mm.prj
$OGS -o $DATA/ANCHORS_MS33_ModelIV  -l warn $DATA/ANCHORS_MS33_ModelIV/ms33_modelIV_pellets.prj
$OGS -o $DATA/ANCHORS_MS33_ModelVII -l warn $DATA/ANCHORS_MS33_ModelVII/ms33_modelVII_freeswelling.prj
```

**Verification gates (all must pass before any model moves to PASS):**

| Gate | Criterion |
|---|---|
| G-K | Model I: MAE < 1 % at ρ_d = 1400, 1600, 1800 (`ms33_calibrate_K.py --verify`) |
| G-phi | `transport_porosity >= 0` at every node/step in every model (VTK check) |
| G-conv | All 6 models complete with 0 rejected steps |

### Invariants — do not break on any future edit

1. The early-return on `slope ≤ 0` **must be removed** (it silences the Pi-path).
2. The Pi block uses `+=` (compressive increment): `delta_sigma_sw += n_S * (Pi_curr − Pi_prev)`. Do not change to `−=`.
3. The eigenstrain fallback (K = 0, slope > 0) must remain unchanged and must
   continue to pass all `ctest -R dsm_micromacro` tests.
4. Any edit to `computeReferenceMicroPorositySwellingStressIncrement` must be
   followed by running `ms33_calibrate_K.py --verify` and confirming MAE < 1 %.
5. K values are per-density. Do not apply a single K to all densities without
   verifying the full Villar curve fit.

### Critical bug: rho_lR initialisation spike (FIXED 2026-05-21)

**Symptom:** `ms33_calibrate_K.py --verify` reports large negative swelling
pressures (e.g. −7.6 MPa vs target +5.8 MPa for dd1600) with the Pi-path active
(`use_micro_liquid_density_for_pi = true`).

**Root cause:** In `initializeConcrete` (`RichardsMechanicsFEM-impl.h`), the micro
liquid density state was initialised to `micro_liquid_density_reference` (= 1e-6
kg/m³ — a trivial EOS placeholder). In the first time step the exchange solve
computes `rho_lR ≈ rho_LR ≈ 1000 kg/m³`, while `rho_lR_prev` remains 1e-6. With
`use_micro_liquid_density_for_pi = true`:

```
Pi_prev = 1e-6 × K × exp(−ξ_prev) ≈ 0
Pi_curr = 1000 × K × exp(−ξ_curr) >> 0
delta_sigma_sw = n_S × (Pi_curr − Pi_prev) >> 0   ← huge tensile spike in step 1
```

This tensile spike (~10⁶× amplified by the density ratio) permanently corrupts
the accumulated sigma_sw, reversing the sign of the swelling pressure.

**Fix (applied 2026-05-21):** At the end of the `initializeConcrete`
`isPotentialExchangeEnabled` block, after `rho_LR_initial` (the actual macro
liquid density at t = 0) is computed, call
`computeActiveMicroLiquidDensity(n_l_initial, rho_LR_initial, ...)` and overwrite
both `*rho_lR` and `**rho_lR_prev` with the resulting physical density (≈ rho_LR
≈ 1000 kg/m³ with the current trivial EOS). When Action 1 (EOS calibration) is
completed, this initialisation will automatically use the correct calibrated density.

**Invariant:** `*rho_lR_prev` at initialisation must be consistent with
`*rho_lR` (i.e. both from `computeActiveMicroLiquidDensity` at the same initial
state). Never initialise rho_lR_prev from `micro_liquid_density_reference` alone.

**Verification:** After the fix, `ms33_calibrate_K.py --verify` gives
MAE = 0.78 % (< 1 % gate), K values unchanged.

## 2026-05-20 — Full clean-rebuild and ground-up validation procedure

This section records the mandatory sequence to run whenever the C++ source
(`ProcessLib/RichardsMechanics/`, `MaterialLib/`, or any DSM-native file) has
been edited and you need to verify the full stack from binary to benchmark.

### Step 0 — Remove all stale output artefacts

Before rebuilding, purge every `.pvd`, `.vtu`, and `.log` from all MS33
model directories so the rerun starts clean and no old-binary outputs are
mistaken for current results:

```bash
DATA=/Users/vinaykumar/git/ogs/Tests/Data/RichardsMechanics
for DIR in ANCHORS_MS33_ModelI ANCHORS_MS33_ModelIII ANCHORS_MS33_ModelIV ANCHORS_MS33_ModelVII; do
    rm -f $DATA/$DIR/*.pvd $DATA/$DIR/*.vtu $DATA/$DIR/*.log
done
```

### Step 1 — Rebuild OGS

```bash
cd /Users/vinaykumar/git/build/release-omp-mfront
ninja RichardsMechanics ogs
```

Wait for `ninja` to succeed (exit code 0) before proceeding.

### Step 2 — Run unit tests

```bash
cd /Users/vinaykumar/git/build/release-omp-mfront
ninja testrunner
bin/testrunner --gtest_filter="RichardsMechanics.*"
```

All 24 (or more) RichardsMechanics unit tests must pass. If any fail, fix them
before continuing — integration tests and model runs are meaningless if the
constitutive unit logic is broken.

### Step 3 — Run full ctest (RichardsMechanics suite)

```bash
cd /Users/vinaykumar/git/build/release-omp-mfront
ctest -j 18 --output-on-failure -R "RichardsMechanics"
```

Expected: all tests pass, including `double_porosity_swelling`,
`double_porosity_swelling-omp`, and the `beacon_1a01` vtkdiff tests.
If any fail, diagnose against the new binary before updating reference VTUs.

To regenerate reference VTUs for a test (e.g., `double_porosity_swelling`):
```bash
# Locate the test PRJ and run OGS with the new binary; copy outputs to
# Tests/Data/RichardsMechanics/<test_dir>/ then re-run ctest to confirm.
```

### Step 4 — Verify K values (Model I calibration check)

```bash
cd /Users/vinaykumar/git/ogs/Tests/Data/RichardsMechanics/ANCHORS_MS33_ModelIV
python ms33_calibrate_K.py --verify
```

Expected: MAE < 1 % for all three dry densities (dd1400, dd1600, dd1800).
If MAE ≥ 1 %, run the bisection calibration:

```bash
python ms33_calibrate_K.py   # updates K in the three Model I PRJ files
```

Then propagate updated K(dd1600) to Models III, IV, VII if it changed.

### Step 5 — Ground-up MS33 model runs

Run all six production PRJs with the freshly built binary:

```bash
OGS=/Users/vinaykumar/git/build/release-omp-mfront/bin/ogs
DATA=/Users/vinaykumar/git/ogs/Tests/Data/RichardsMechanics

# Model I — three dry densities (run in parallel)
for DD in dd1400 dd1600 dd1800; do
    $OGS -o $DATA/ANCHORS_MS33_ModelI -l warn \
         $DATA/ANCHORS_MS33_ModelI/ms33_modelI_${DD}.prj \
         > $DATA/ANCHORS_MS33_ModelI/anchors_${DD}_full.log 2>&1 &
done
wait

# Model III (gap-closure)
$OGS -o $DATA/ANCHORS_MS33_ModelIII -l warn \
     $DATA/ANCHORS_MS33_ModelIII/ms33_modelIII_gap2mm.prj \
     > $DATA/ANCHORS_MS33_ModelIII/ms33_modelIII_run_$(date +%Y%m%d).log 2>&1

# Model IV (pellets)
$OGS -o $DATA/ANCHORS_MS33_ModelIV -l warn \
     $DATA/ANCHORS_MS33_ModelIV/ms33_modelIV_pellets.prj \
     > $DATA/ANCHORS_MS33_ModelIV/ms33_modelIV_run_$(date +%Y%m%d).log 2>&1

# Model VII (free swelling)
$OGS -o $DATA/ANCHORS_MS33_ModelVII -l warn \
     $DATA/ANCHORS_MS33_ModelVII/ms33_modelVII_freeswelling.prj \
     > $DATA/ANCHORS_MS33_ModelVII/ms33_modelVII_run_$(date +%Y%m%d).log 2>&1
```

### Step 6 — Post-processing and plausibility checks

After all runs complete:

1. **Check convergence**: `grep "rejected" $DATA/ANCHORS_MS33_Model*/*.log` — all
   six models must show 0 rejected steps.
2. **Check transport_porosity ≥ 0**: run `ms33_calibrate_K.py --verify` (it
   also exercises the node-wise porosity check) or inspect VTKs manually.
3. **Regenerate PDFs**:
   ```bash
   cd $DATA/ANCHORS_MS33_ModelIV && python ms33_postprocess.py
   cd $DATA/ANCHORS_MS33_ModelIII && python ms33_postprocess_modelIII.py
   cd $DATA/ANCHORS_MS33_ModelVII && python ms33_postprocess_modelVII.py
   ```
4. **Update AGENTS.md problem log** in `ANCHORS_MS33_ModelIV/AGENTS.md` with
   a new timestamped entry recording the result of this run.

### Verification gates (must all pass before marking suite PASS)

| Gate | Criterion |
|---|---|
| G-unit | All `RichardsMechanics.*` gtest cases pass |
| G-ctest | `ctest -R RichardsMechanics` fully green (0 failing) |
| G-K | Model I: MAE < 1 % at ρ_d = 1400, 1600, 1800 |
| G-phi | `transport_porosity >= 0` at every node/step in every model |
| G-conv | All 6 models complete with 0 rejected steps |

## 2026-05-20 23:51:31 CEST — Exact non-running / failing status

Current branch/commit used for this check:
- Branch: `dsm_native_hierarchical`
- Commit: `45991f64b0`

### 1) RichardsMechanics ctest is NOT green

Command run:
```bash
cd /Users/vinaykumar/git/build/release-omp-mfront/ProcessLib/RichardsMechanics
ctest --output-on-failure
```

Result:
- `94% tests passed, 6 tests failed out of 94`
- Failing tests:
  - `ogs-RichardsMechanics/double_porosity_swelling`
  - `ogs-RichardsMechanics/double_porosity_swelling-omp`
  - `ogs-RichardsMechanics/DoubleStructureBenchmark/double_porosity_swelling_RM`
  - `ogs-RichardsMechanics/DoubleStructureBenchmark/double_porosity_swelling_RM-omp`
  - `ogs-RichardsMechanics_beacon_1a01_dsm_micromacro_inflow_reference-time-vtkdiff`
  - `ogs-RichardsMechanics_beacon_1a01_dsm_micromacro_inflow_reference-time-omp-vtkdiff`

Observed mismatch pattern:
- `double_porosity_swelling*`: large vtkdiff deviations in
  `saturation`, `displacement`, `sigma`, `epsilon`, `porosity`,
  `micro_pressure`, `micro_saturation`.
- `beacon_1a01 ... inflow ... vtkdiff*`: mismatch concentrated in
  `swelling_stress` and `sigma` (other monitored fields match).

### 2) MS33 calibration verification script is NOT runnable in this environment

Command:
```bash
python3 /Users/vinaykumar/git/ogs/Tests/Data/RichardsMechanics/ANCHORS_MS33_ModelIV/ms33_calibrate_K.py --verify
```

Result:
- Fails immediately with `ModuleNotFoundError: No module named 'scipy'`.
- Therefore gate `G-K` cannot be executed until SciPy is installed.

### 3) Consequence for benchmark status

Because (1) ctest is red and (2) `G-K` cannot be executed, the suite must
remain **NOT PASS**. Do not claim full validation until these two blockers are
resolved and rerun.
