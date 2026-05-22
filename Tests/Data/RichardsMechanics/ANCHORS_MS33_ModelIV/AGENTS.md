# AGENTS.md — ANCHORS MS33 DSM Pi-path benchmarks

**Scope:** EURAD-2 MS33 theoretical benchmarking suite, BGR DSM native
hierarchical vdW-augmented model. Covers Models I / III / IV / VII and the
shared implementation file
`ProcessLib/RichardsMechanics/RichardsMechanicsFEM-impl.h`.

**OGS binary (release-omp-mfront):**
```
/Users/vinaykumar/git/build/release-omp-mfront/bin/ogs
```
Do **not** use the `dsm_native-release` worktree binary — VTK ABI mismatch
(links vtk 9.5, system has 9.6).

**Branch:** `dsm_native_hierarchical` (remote: origin)

---

## Roadmap summary

| Step | Description | Status |
|------|-------------|--------|
| 1 | REV-scale n_l mass balance (φ_m·ρ_lR, not n_l·ρ_lR); cross-coupling implicit in storage diff | ✅ done — commit `0d7a9edd64` (2026-05-20) |
| 2 | Thermodynamic swelling stress σ_sw = −φ_m·Π; K recalibration | ✅ done — commit `88d42c98fd` (2026-05-21) |
| 3 | EOS argument unification: decouple ρ_lR(EOS) from Π; Gibbs–Duhem consistency | ✅ done — see Action 6 completion (2026-05-21) |
| 4 | Comment cleanup: remove additive-approximation text from PRJ files and paper | ✅ done — 2026-05-21 |
| 5 | Fix vdW base dimensional error: /ρ_lR added, A=2.2e-20 J (literature), K recalibrated | ✅ done — 2026-05-22 |
| 6 | DSM consistency hardening: viscosity guards, micro-pressure density default, test updates, LE rerun verification | ✅ done — commit `66b782afa1` (2026-05-22) |

---

## 2026-05-22 14:10 CEST — Action 7: DSM consistency hardening + full LE rerun verification

**Commit:** `66b782afa1`

**What was done:**

1. **Potential-exchange derivative tests updated**  
   `Tests/ProcessLib/RichardsMechanics/PotentialExchange.cpp` now checks
   finite-difference consistency for `dmu_lR/drho_lR` (non-zero after vdW
   dimensional correction), instead of legacy zero-derivative expectations.

2. **Positive-viscosity guard applied across DSM exchange paths**  
   In `RichardsMechanicsFEM-impl.h`, added shared
   `requirePositiveViscosity(...)` and applied it to:
   - predictor/coupled micro mass-storage solves,
   - implicit `n_l` solve,
   - implicit `dn_l/dp_L` helper,
   - local exchange Jacobian assembly paths using `alpha_bar * rho_LR / mu`.

3. **Micro-pressure density default corrected**  
   `use_micro_liquid_density_for_micro_pressure` default switched to `true`
   (confined micro-liquid density basis) in:
   - `PotentialExchangeParameters.h`,
   - `CreateRichardsMechanicsProcess.cpp`.
   If explicitly set `false`, a runtime warning is issued.

4. **Build + rerun verification**  
   Rebuilt `release-omp-mfront` (`--parallel 18`) and reran canonical LE set:
   - Model I dd1400/dd1600/dd1800: 308/0 each
   - Model III: 331/0
   - Model IV: 313/0
   - Model V LE: 308/0
   - Model VII: 433/0
   Total LE accepted steps: 2309, rejected: 0.

5. **Artefacts committed**  
   Source changes, tests, and regenerated run outputs/logs were committed in
   `66b782afa1`.

---

## 2026-05-20 — Action 1–4: nS fix, n_l balance, Model VII IC, K calibration scaffold

**Commit:** `0d7a9edd64`

**What was done:**

1. **nS fix** (`computeActiveMicroSolidVolumeFraction`): corrected to return
   `1 − φ_M` in `MicroSolidVolumeFractionMode::Reference` mode instead of a
   wrong formula that gave values > 1.

2. **REV-scale n_l balance** (`solveReferenceMassStoragePredictorState` and
   `...CoupledState`): storage term changed from `n_l · ρ_lR` (aggregate) to
   `(1−φ_M)/(1−n_l) · n_l · ρ_lR` (REV-scale `φ_m · ρ_lR`). This is the
   hierarchical split `φ_m = (1−φ_M) · n_l`.

3. **Cross-coupling** (`−φ_m · φ̇_M / (1−φ_M)` in n_l balance): proved absent
   from code is correct — the backward-Euler storage difference
   `φ_m^{n+1}·ρ_lR^{n+1} − φ_m^n·ρ_lR^n` implicitly contains both the micro
   uptake term and the cross-coupling term when φ_M is encoded at both time
   levels via the hierarchical split. No separate code term needed. Documented
   in FIX-4 slide of `nagel_porosity_split.tex`.

4. **Model VII IC**: sigma0 corrected to `[-3.293e7, -3.313e7, -3.293e7, 0]` Pa
   to balance Biot pressure force at free boundary.

5. **All 6 PRJ files**: abstols changed from `5e-8 1e-13 1e-13` → `1 1e-12 1e-12`.

**Verification gate (step 1):** All 6 models run 0 rejected steps with K values
from this era. Model I Villar errors were then ~28–43% (expected — step 2 not
yet done).

---

## 2026-05-21 — Action 5: Step 2 — thermodynamic swelling stress formula and K recalibration

**Commit:** `88d42c98fd`

### Code change

**File:** `ProcessLib/RichardsMechanics/RichardsMechanicsFEM-impl.h`  
**Function:** `computeReferenceMicroPorositySwellingStressIncrement` (line ~1497)

**Old formula (both legacy and active paths):**
```cpp
delta_sigma_sw.noalias() += n_S * (Pi_curr - Pi_prev) * identity2;
```

**New formula (both paths):**
```cpp
// Thermodynamic form: sigma_sw = -phi_m * Pi = -n_S * n_l * Pi
// (compressive eigenstress, tension-positive convention).
// Sign: n_l*Pi increases during hydration (n_l << C = lambda*nS*rho_SR*Sa),
// so n_l_prev*Pi_prev - n_l*Pi_curr < 0 -> compressive increment.
delta_sigma_sw.noalias() +=
    n_S * (n_l_prev * Pi_prev - n_l * Pi_curr) * identity2;
```

**Why the sign works:** At equilibrium `n_l · Π = n_l · ρ_lR · K · exp(−ξ)`.
Since `n_l ≪ C = λ · nS · ρ_SR · Sa ≈ 0.628` for dd1600 parameters, `n_l · Π`
is an increasing function of `n_l`. Hydration increases `n_l`, so the increment
`n_l_prev · Π_prev − n_l · Π_curr < 0` → compressive, consistent with swelling.

### K recalibration

**Script:** `ms33_calibrate_K.py` (this directory)  
**Method:** `scipy.optimize.brentq` bisection on residual `p_sw − p_Villar`  
**Villar target:** `p_sw = exp(6.77 · ρ_d[g/cm³] − 9.07)` MPa

| ρ_d [kg/m³] | K_old [J/kg] | K_new [J/kg] | Ratio | Villar [MPa] | err [%] |
|------------|-------------|-------------|-------|-------------|---------|
| 1400 | 3622.66 | 6963.70 | ×1.92 | 1.50381 | −0.009 |
| 1600 | 16792.03 | 27371.90 | ×1.63 | 5.82407 | −0.011 |
| 1800 | 75853.21 | 107246.59 | ×1.41 | 22.55598 | +0.005 |

K_dd1600 = 27371.9007 propagated to ModelIII, IV, VII (same reference density).

**PRJ files updated (all 6):**
```
ANCHORS_MS33_ModelI/ms33_modelI_dd1400.prj   ← K=6963.6999
ANCHORS_MS33_ModelI/ms33_modelI_dd1600.prj   ← K=27371.901
ANCHORS_MS33_ModelI/ms33_modelI_dd1800.prj   ← K=107246.59
ANCHORS_MS33_ModelIII/ms33_modelIII_gap2mm.prj  ← K=27371.9007
ANCHORS_MS33_ModelIV/ms33_modelIV_pellets.prj   ← K=27371.9007
ANCHORS_MS33_ModelVII/ms33_modelVII_freeswelling.prj ← K=27371.9007
```

**Verification gate (step 2):**

| Model | p_sw [MPa] | Villar [MPa] | err [%] | Note |
|-------|-----------|-------------|---------|------|
| I dd1400 | 1.50368 | 1.50381 | −0.009 | calibration target |
| I dd1600 | 5.82345 | 5.82407 | −0.011 | calibration target |
| I dd1800 | 22.55698 | 22.55598 | +0.005 | calibration target |
| III gap2mm | 4.88266 | (5.82407) | −16.2 | expected: gap lowers mean |
| IV pellets | 4.29377 | (5.82407) | −26.3 | expected: pellet zone lower ρ_d |
| VII freeswelling | 2.23 | — | — | free boundary; no constraint |

All 6 models: 0 rejected time steps.

---

## 2026-05-21 — Action 6: Step 3 — EOS argument unification and Gibbs–Duhem fix (COMPLETED)

**Status:** Completed on 2026-05-21 (commit: this EOS-fix commit).  
**Scope:** Decouple swelling-pressure Pi path from `rho_lR(EOS)` and make Pi use
bulk density `rho_LR`, then recalibrate `K` and rerun full MS33 suite.

### Background

The micro physics has two empirical EOS objects that use different coordinate
variables for the same physical state n_l:

**Density EOS** (`computeReducedMicroLiquidDensity`,
`RichardsMechanicsFEM-impl.h` line ~440):
```
ω_l = n_l · ρ_lR / (n_S · ρ_SR)          [gravimetric, kg_water/kg_solid]
ρ_lR = ρ_LR + ρ_l0 · exp(−a · ω_l^b)     [implicit fixed-point, Newton]
```
`ρ_lR` appears in its own argument via ω_l → requires Newton iteration to
converge.

**Pi-path chemical potential** (`computeVanDerWaalsMicroPotential`,
`ConstitutiveRelations/PotentialExchange.h` line ~81):
```
ξ = n_l / (λ · n_S · ρ_SR · S_a)         [geometric, h/λ, no ρ_lR]
μ_aug = K · exp(−ξ)                       [J/kg]
```
Code explicitly: `out.dmu_lR_drho_lR = 0.0;` (line 176) — potential is blind
to ρ_lR.

**Gibbs–Duhem violation:** The swelling stress uses `Π = ρ_lR · μ_aug`.
For thermodynamic consistency: `dμ = (1/ρ_lR) · dΠ`. Substituting:
`(1/ρ_lR) · ∂Π/∂ρ_lR = K · exp(−ξ) ≠ 0 = ∂μ_aug/∂ρ_lR`.
The relation is violated in the ρ_lR direction.

**Parametrisation (MS33 PRJ files, before fix):**
```xml
<micro_liquid_density_reference>100.0</micro_liquid_density_reference>
<micro_liquid_density_a>1e-16</micro_liquid_density_a>
<micro_liquid_density_b>1.0</micro_liquid_density_b>
```
With `a = 1e-16 ≈ 0`, the density EOS reduces to `ρ_lR ≈ ρ_LR + 100 ≈ 1100 kg/m³`
(constant). So Π used a constant 1100 kg/m³ prefactor. The Gibbs–Duhem
violation was numerically small (~10%) but conceptually wrong.

### The fix

**Principle:** separate the mechanical pressure (Π [Pa]) from the specific
chemical potential (μ [J/kg]). Use bulk density ρ_LR (not the structured micro
density ρ_lR) in the Π ↔ μ conversion. This satisfies Gibbs–Duhem exactly:
`dμ = (1/ρ_LR) · dΠ` ✓ (since ρ_LR is constant in the MS33 setup).

The density EOS (ρ_lR) is retained for the **mass balance** storage term
`ρ_l = φ_m · ρ_lR`, where structured water density is physically meaningful.
It is removed from the **swelling stress** calculation.

**Code change (DONE — commit 2026-05-21):**

In `computeReferenceMicroPorositySwellingStressIncrement`
(`RichardsMechanicsFEM-impl.h`, lines 1531–1562), both code paths use:
```cpp
double const Pi_curr = rho_LR * K * std::exp(-xi_curr);
double const Pi_prev = rho_LR * K * std::exp(-xi_prev);
```
ρ_lR (micro EOS density) is no longer used in the Pi-path. Gibbs–Duhem is
satisfied exactly since ρ_LR is constant in the MS33 setup.

**`use_micro_liquid_density_for_pi` flag removed (2026-05-22):**
The flag was stored in `PotentialExchangeParameters` and parsed from PRJ, but
was never branched on in the implementation — dead code. Removed from:
- `PotentialExchangeParameters.h`
- `CreateRichardsMechanicsProcess.cpp`
- `DSMMicroMacroSingleIntegrationPoint.cpp` (unit test)
- All MS33 PRJ files that had `<use_micro_liquid_density_for_pi>true</use_micro_liquid_density_for_pi>`

### K recalibration after the fix

After the code change, the swelling stress prefactor drops from ρ_lR ≈ 1100
to ρ_LR ≈ 1000 kg/m³, reducing Π by ~9.1%. To restore the Villar targets,
K must increase by 1100/1000 = 1.10:

| ρ_d | K_current [J/kg] | K_expected [J/kg] |
|-----|-----------------|-------------------|
| 1400 | 6963.70 | ~7660 |
| 1600 | 27371.90 | ~30110 |
| 1800 | 107246.59 | ~117970 |

Run `ms33_calibrate_K.py` (this directory) to get exact values. The script
uses brentq bisection and converges to < 0.015% Villar error in ~7 OGS calls
per density.

### Optional Part B: EOS argument unification (paper clarity)

After the code fix above, ρ_lR enters only the mass balance. The argument
choice (ω_l vs ξ) no longer affects the swelling stress. However, for paper
consistency, switch the density EOS to use ξ (same variable as μ_aug):

```
ρ_lR = ρ_LR + A · exp(−B · ξ^C)    [explicit, no Newton iteration]
ξ = n_l / (λ · n_S · ρ_SR · S_a)
```

This change is in `computeReducedMicroLiquidDensity` (lines 440–533). Replace
the `ω_l`-based implicit solve with an explicit evaluation in ξ. Parameters
(A, B, C) must be refit to adsorption isotherm data under the new argument
(not interchangeable with the current ω_l-based (A=100, B=1e-16, C=1)).

Part B is lower priority — the current (A=100, B≈0) parametrisation already
yields near-constant ρ_lR ≈ 1100, so the ω_l vs ξ distinction is numerically
negligible with current parameters.

### Verification gate (step 3)

**Implemented and verified on 2026-05-21:**
1. Compile cleanly: `cmake --build /Users/vinaykumar/git/build/release-omp-mfront -j8` ✅
2. K recalibration complete (Model I):
   - dd1400: `K = 7656.5016`, `p_sw = 1.50378 MPa`, Villar err `-0.0020%`
   - dd1600: `K = 29999.2513`, `p_sw = 5.82306 MPa`, Villar err `-0.0174%`
   - dd1800: `K = 118585.8600`, `p_sw = 22.56225 MPa`, Villar err `+0.0278%`
   - mean absolute Villar error: `0.0158%` ✅
3. Full six-model MS33 suite rerun: all converged with `0` rejected steps ✅
   - I dd1400: accepted/rejected = `308/0`
   - I dd1600: accepted/rejected = `308/0`
   - I dd1800: accepted/rejected = `308/0`
   - III gap2mm: accepted/rejected = `299/0`
   - IV pellets: accepted/rejected = `279/0`
   - VII freeswelling: accepted/rejected = `458/0`
4. Final-time swelling pressure sign check (`p_sw > 0`) ✅
   - III: `5.493395 MPa`
   - IV: `4.294508 MPa`
   - VII: `4.746589 MPa`

---

## Action 7 (DONE): Step 4 — comment cleanup

**Date completed:** 2026-05-21

Updated all 6 PRJ files:
- Header K values (lines ~8): dd1400 6963.70→7656.50, dd1600 27371.90→29999.25,
  dd1800 107246.59→118585.86 J/kg
- Inline K comments: replaced "recalibrated with thermodynamic σ_sw / ρ_lR=1100"
  → "Π uses bulk ρ_LR (Gibbs–Duhem consistent); steps 2+3, 2026-05-21"
- ModelIV Pi formula comment: corrected to step-2/3 formula
  `nS*(n_l_prev*Pi_prev - n_l*Pi_curr)*I` with `Pi = rho_LR*K*exp(-xi)`
- Removed "Pending recalibration" language from ModelIV
Presentation `nagel_porosity_split.tex`:
- Step 3 marked done in FIX-5, bottom-line roadmap, FIX-3 alertblock, formula-status slide

---

## Quick-reference: PRJ parameter locations

All 6 PRJ files share the same `<potential_exchange>` block structure:

```
Tests/Data/RichardsMechanics/
  ANCHORS_MS33_ModelI/
    ms33_modelI_dd1400.prj    line ~52–75
    ms33_modelI_dd1600.prj    line ~52–75
    ms33_modelI_dd1800.prj    line ~52–75
  ANCHORS_MS33_ModelIII/
    ms33_modelIII_gap2mm.prj  line ~69–88
  ANCHORS_MS33_ModelIV/
    ms33_modelIV_pellets.prj  line ~73–92
  ANCHORS_MS33_ModelVII/
    ms33_modelVII_freeswelling.prj  line ~80–98
```

Key parameters and current values (as of Action 6 completion, 2026-05-21):

| Parameter | Value | Role |
|-----------|-------|------|
| `micro_liquid_density_reference` | 100.0 kg/m³ | ρ_l0 in density EOS (constant excess above bulk) |
| `micro_liquid_density_a` | 1e-16 | exp decay rate — near zero makes EOS constant |
| `micro_liquid_density_b` | 1.0 | exponent in ω_l^b |
| `hamaker_constant` | 2.2e-20 J | Hamaker constant A; Israelachvili & Adams 1978; **do NOT calibrate** |
| `vdw_augmentation_prefactor` | see table above | K [J/kg] |
| `vdw_augmentation_decay_length` | 7.5e-7 m | λ = 0.75 nm |
| `micro_solid_density_reference` | 2780.0 kg/m³ | ρ_SR |
| `specific_surface` | 523 m²/kg | S_a (BET surface area, MX-80 montmorillonite) |

---

## Calibration script

`ms33_calibrate_K.py` (this directory):

```bash
# calibrate all three densities:
python ms33_calibrate_K.py

# verify existing K values without changing PRJ:
python ms33_calibrate_K.py --verify

# calibrate one density only:
python ms33_calibrate_K.py --density 1600
```

Results are written to `ms33_calibrate_K_results.txt`.

After calibration, propagate K_dd1600 manually to Models III/IV/VII:
```bash
K=<new_K_dd1600>
sed -i '' "s|<vdw_augmentation_prefactor>[^<]*</vdw_augmentation_prefactor>|<vdw_augmentation_prefactor>${K}</vdw_augmentation_prefactor>|g" \
  ../ANCHORS_MS33_ModelIII/ms33_modelIII_gap2mm.prj \
  ../ANCHORS_MS33_ModelIV/ms33_modelIV_pellets.prj \
  ../ANCHORS_MS33_ModelVII/ms33_modelVII_freeswelling.prj
```

---

## STRICT INVARIANT: vdW base potential must never be removed or replaced

**File:** `ProcessLib/RichardsMechanics/ConstitutiveRelations/PotentialExchange.h`
**Function:** `computeVanDerWaalsMicroPotential`

### What the code does (correct, dimensionally consistent after step 5)

```cpp
// 1. Hamaker vdW base potential — SET first [J/kg] (step 5 fix: divides by rho_lR)
double const prefactor = hamaker_constant * Sa³ / (6π);
out.mu_lR = potential_sign_factor * prefactor * nS³ * rho_SR³ / (n_l³ * rho_lR);
out.dmu_lR_drho_lR = -out.mu_lR / rho_lR;  // non-zero after /rho_lR fix

// 2. Exponential augmentation — ADDED on top [J/kg]
if (vdw_augmentation_prefactor > 0.0) {
    double const mu_aug = potential_sign_factor * K * std::exp(-xi);
    out.mu_lR += mu_aug;   // ← += not =, MUST stay additive ALWAYS
    out.dmu_lR_dnl += -mu_aug * xi / n_l;
    // augmentation does not depend on rho_lR; dmu_lR_drho_lR unchanged
    out.dmu_lR_dnS += mu_aug * xi / nS;
    out.dmu_lR_drho_SR += mu_aug * xi / rho_SR;
}
```

### What must NEVER happen

**NEVER** change `out.mu_lR += mu_aug` to `out.mu_lR = mu_aug`.
Doing so would silently discard the Hamaker vdW base potential and make
`K` carry ALL surface physics — the augmentation would no longer augment
anything. This is undetectable from calibration results alone (K absorbs
the larger mismatch) but invalidates the physical interpretation of both A and K.

### Hamaker constant A: material constant from literature — do NOT calibrate

A is the Hamaker constant for clay–water–clay interaction [J].
**It must be set from literature data, NOT from swelling-pressure calibration.**
Only K (augmentation amplitude) and λ (augmentation decay length) are calibrated.

| System | A (J) | Source |
|--------|-------|--------|
| Mica–water–mica (direct SFA) | 2.2×10⁻²⁰ | Israelachvili & Adams (1978) |
| Montmorillonite–water–montm. (DLVO) | 1–5×10⁻²⁰ | Novich & Ring (1984) |
| Na-montmorillonite (FHH isotherm) | ~1.5×10⁻²⁰ | Cases et al. (1992) |
| Smectite (DFT/Lifshitz) | 1–3×10⁻²⁰ | Šolc et al. (2011) |
| Silica–water–silica | 0.8–1.0×10⁻²⁰ | Bergström (1997) |
| **Value used (all MS33 PRJs)** | **2.2×10⁻²⁰** | Israelachvili & Adams 1978 (mica proxy for smectite) |

**PRJ tag:** `<hamaker_constant>2.2e-20</hamaker_constant>`

### Dimensional consistency (all three potentials must be J/kg)

```
μ_LR     = p_L / ρ_LR                                      [J/kg]  macro (Young-Laplace)
μ_lR_vdW = A·Sa³·nS³·ρ_SR³ / (6π·n_l³·ρ_lR)              [J/kg]  micro vdW base
μ_lR_aug = K·exp(−ξ)                                        [J/kg]  micro augmentation
μ_lR     = μ_lR_vdW + μ_lR_aug                             [J/kg]  total micro
Exchange: ρ̂_l = α_M_eff·(μ_LR − μ_lR)                   [kg/(m³·s)]
```

Proof: `p_L_m = −ρ·μ_lR` (impl.h lines 276, 1044) → μ_lR MUST be [J/kg].

### Which ρ_lR enters the vdW denominator

The formula is the **specific free energy** of the water film:

```
μ_lR_vdW = (surface energy/area) × (solid surface area / REV volume)
                                  / (water mass / REV volume)
          = E [J/m²] × nS·ρ_SR·Sa [m⁻¹]  /  n_l·ρ_lR [kg/m³]
```

The ρ_lR in the denominator is the **density of the confined water in the micro
pore** — the micro EOS density (~1100 kg/m³), not the bulk liquid density (~1000 kg/m³).

**Active code path** (`computeCompatibilityMicroHydraulicOutput` with local_context,
all MS33 call sites): uses micro density correctly via `computeActiveMicroPotential`
→ `computeReducedMicroLiquidDensity` (for `ScalarReferenceMassStorage` mode).
`p_L_m = −ρ_lR·μ_lR` also uses micro density when
`use_micro_liquid_density_for_micro_pressure = true` (set in all MS33 PRJs).

A 3-argument overload without `local_context` previously existed at impl.h line ~257
and used bulk ρ_LR in the vdW formula — ~10% error. It was dead code (never called)
and was removed on 2026-05-22.

### Consequence for adding new physical potentials

K is now calibrated to Villar swelling pressure AFTER the vdW base is
physically correct. A new physical term (osmotic, double-layer, structural
hydration) can be added as:
```cpp
out.mu_lR += new_term;   // must be in J/kg
```
Recalibrate K to Villar after adding any new term, since K carries only
the residual beyond vdW+new_term.

---

## Action 8 (DONE): Step 5 — vdW base dimensional fix

**Date completed:** 2026-05-22

### What was done

1. **Code fix** in `PotentialExchange.h`, `computeVanDerWaalsMicroPotential`:
   - Added `/ rho_lR` to vdW base formula → units now J/kg (was Pa)
   - Updated `dmu_lR_drho_lR = -out.mu_lR / rho_lR` (was 0.0)
   - Updated comment block with dimensional derivation and A literature values

2. **A updated** in all 27 PRJ files:
   - FROM: `<hamaker_constant>5.1e-21</hamaker_constant>` (numerically tuned, wrong units)
   - TO:   `<hamaker_constant>2.2e-20</hamaker_constant>` (Israelachvili & Adams 1978)
   - A is NOT calibrated; it is a material constant.

3. **K recalibrated** via `ms33_calibrate_K.py` (Villar bisection):

   | ρ_d (kg/m³) | K (J/kg) | p_sw,sim (MPa) | Villar (MPa) | Error |
   |-------------|----------|----------------|--------------|-------|
   | 1400 | 7 654.9 | 1.50376 | 1.50381 | −0.003% |
   | 1600 | 29 984.9 | 5.82286 | 5.82407 | −0.021% |
   | 1800 | 118 582.6 | 22.56322 | 22.55598 | +0.032% |

   K values are within 0.05% of step 3 values — the vdW contribution (~620 J/kg
   with A=2.2e-20 at MS33 initial state) is ~0.6% of the macro driving force
   (−10⁵ J/kg at 100 MPa), so K is essentially unchanged.

4. **All 6 MS33 models** run successfully, 0 rejected steps.
   - Model III (gap 2 mm): p_sw = 5.47 MPa at t=200 d; transport_porosity fills gap.
   - Model IV (pellets): p_sw = 4.25 MPa; transport_porosity = 0 (fully expanded).
   - Model VII (free swelling): axial displacement 17.2 mm (24.5% linear strain).

### For future agents: dimensional consistency check procedure

When modifying `computeVanDerWaalsMicroPotential` or adding new potential terms:

1. **Verify units** of every term entering `mu_lR`:
   - Each term must be [J/kg]
   - Check: does dividing [term units] by [kg/m³] give [J/kg]? If not, fix.

2. **Check proof**: `p_L_m = −rho_density * mu_lR` in impl.h lines 276 and 1044.
   rho_density ≈ 1000 kg/m³ and p_L_m must be in Pa → mu_lR must be J/kg.

3. **Augmentation is additive** (`+= mu_aug`, NEVER `= mu_aug`).

4. **A is not calibrated**. Current value: 2.2e-20 J (Israelachvili & Adams 1978).
   If A needs updating, require a literature source; do not bisect to Villar.
   After any A change, re-run `ms33_calibrate_K.py` to recalibrate K.

5. **K calibration** after any physics change:
   ```bash
   cd Tests/Data/RichardsMechanics/ANCHORS_MS33_ModelIV
   python ms33_calibrate_K.py
   ```
   Then propagate dd1600 K to Models III, IV, VII (all 15 ModelVII PRJs).
   Verification gate: Villar errors < 0.1% for all three densities.

6. **Run all 6 models** and confirm 0 rejected steps:
   ```bash
   OGS=/Users/vinaykumar/git/build/release-omp-mfront/bin/ogs
   $OGS -o . -l warn ms33_modelI_dd1400.prj
   $OGS -o . -l warn ms33_modelI_dd1600.prj
   $OGS -o . -l warn ms33_modelI_dd1800.prj
   $OGS -o . -l warn ../ANCHORS_MS33_ModelIII/ms33_modelIII_gap2mm.prj
   $OGS -o . -l warn ../ANCHORS_MS33_ModelIV/ms33_modelIV_pellets.prj
   $OGS -o . -l warn ../ANCHORS_MS33_ModelVII/ms33_modelVII_freeswelling.prj
   ```

---

## Related memory files

- `~/.claude/projects/-Users-vinaykumar-git-ogs/memory/project_ms33_benchmark_prj_status.md`
- `~/.claude/projects/-Users-vinaykumar-git-ogs/memory/feedback_ogs_rm_sigma0_free_boundary.md`
- `~/.claude/projects/-Users-vinaykumar-git-ogs/memory/feedback_ogs_rm_pressure_tolerance.md`

## Related presentation

`/Users/vinaykumar/tex/cc2024/VK_B35_Pinion_May_2026/nagel_porosity_split.tex`

FIX slides: FIX-1 (n_l balance), FIX-2 (code changes), FIX-3 (impact),
FIX-4 (cross-coupling implicit), FIX-5 (EOS argument — step 3, 2026-05-21).
