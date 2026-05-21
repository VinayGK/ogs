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
| 3 | EOS argument unification: decouple ρ_lR(EOS) from Π; Gibbs–Duhem consistency | ⬜ open — see Action 6 below |
| 4 | Comment cleanup: remove additive-approximation text from PRJ files and paper | ⬜ open |

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

## Action 6 (OPEN): Step 3 — EOS argument unification and Gibbs–Duhem fix

**Date open:** 2026-05-21  
**Priority:** Medium — affects time-path accuracy and paper correctness; calibration  
currently absorbs the inconsistency at equilibrium.

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

**Current parametrisation (MS33 PRJ files):**
```xml
<micro_liquid_density_reference>100.0</micro_liquid_density_reference>
<micro_liquid_density_a>1e-16</micro_liquid_density_a>
<micro_liquid_density_b>1.0</micro_liquid_density_b>
<use_micro_liquid_density_for_pi>true</use_micro_liquid_density_for_pi>
```
With `a = 1e-16 ≈ 0`, the density EOS reduces to `ρ_lR ≈ ρ_LR + 100 ≈ 1100 kg/m³`
(constant). So Π currently uses a constant 1100 kg/m³ prefactor. The Gibbs–Duhem
violation is numerically small (~10%) but conceptually wrong.

### The fix

**Principle:** separate the mechanical pressure (Π [Pa]) from the specific
chemical potential (μ [J/kg]). Use bulk density ρ_LR (not the structured micro
density ρ_lR) in the Π ↔ μ conversion. This satisfies Gibbs–Duhem exactly:
`dμ = (1/ρ_LR) · dΠ` ✓ (since ρ_LR is constant in the MS33 setup).

The density EOS (ρ_lR) is retained for the **mass balance** storage term
`ρ_l = φ_m · ρ_lR`, where structured water density is physically meaningful.
It is removed from the **swelling stress** calculation.

**Code change required:**

In `computeReferenceMicroPorositySwellingStressIncrement`
(`RichardsMechanicsFEM-impl.h`, lines 1497–1603), the function signature and
body need one change:

1. **Add `rho_LR` (bulk macro density) as a parameter** (currently not passed):
   ```cpp
   // Old signature:
   computeReferenceMicroPorositySwellingStressIncrement(
       double const n_l_prev, double const n_l,
       double const n_S, double const rho_lR, double const rho_lR_prev, ...)

   // New signature — add rho_LR:
   computeReferenceMicroPorositySwellingStressIncrement(
       double const n_l_prev, double const n_l,
       double const n_S, double const rho_lR, double const rho_lR_prev,
       double const rho_LR, ...)   // ← add this
   ```

2. **Replace state-varying ρ_lR with ρ_LR in Π** (both code paths, lines ~1544
   and ~1579):
   ```cpp
   // Old (both legacy and active paths):
   double const rho_curr = rho_lR;
   double const rho_prev =
       potential_exchange_params.use_micro_liquid_density_for_pi
           ? rho_lR_prev : rho_lR;
   double const Pi_curr = rho_curr * K * std::exp(-xi_curr);
   double const Pi_prev = rho_prev * K * std::exp(-xi_prev);

   // New:
   double const Pi_curr = rho_LR * K * std::exp(-xi_curr);
   double const Pi_prev = rho_LR * K * std::exp(-xi_prev);
   ```
   The `use_micro_liquid_density_for_pi` flag becomes redundant for the Pi-path
   augmentation and can be left in place for the vdW term if needed.

3. **Update all call sites** of `computeReferenceMicroPorositySwellingStressIncrement`
   and the wrapper `computeSwellingStressIncrement` to pass `rho_LR`. The call
   site in the assembly loop already has `rho_LR` available (it is used in the
   storage term assembly on the same integration point).

4. **Update `updateSwellingState`** (line ~1620) which calls
   `computeSwellingStressIncrement` — pass `rho_LR` through there too.

**No PRJ changes needed for the code fix.** The `use_micro_liquid_density_for_pi`
flag is rendered inert for the Pi-path by the code change and can be cleaned up
separately (step 4).

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

After code change + K recalibration:
1. All 3 Model I densities: Villar error < 0.1% (< 0.015% expected)
2. All 6 models: 0 rejected time steps
3. Compile cleanly: `ninja -j8 ogs` in `/Users/vinaykumar/git/build/release-omp-mfront`
4. Confirm `p_sw > 0` for all 6 models at final timestep

---

## Action 7 (OPEN): Step 4 — comment cleanup

**Date open:** 2026-05-21

Remove all remaining references to the additive-approximation model from:
1. PRJ file header comments (search for "additive", "K_code", "old formula",
   "absorb")
2. Paper draft: `/Users/vinaykumar/tex/dsm-bgr-paper/draft/paper_DSM.tex` —
   update EOS argument notation; document that Π uses bulk density; cite
   Gibbs–Duhem argument
3. Presentation: `nagel_porosity_split.tex` — steps 3/4 are already marked
   open in the roadmap slides (FIX-5 added 2026-05-21)

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

Key parameters and current values (as of commit `88d42c98fd`):

| Parameter | Value | Role |
|-----------|-------|------|
| `micro_liquid_density_reference` | 100.0 kg/m³ | ρ_l0 in density EOS (constant excess above bulk) |
| `micro_liquid_density_a` | 1e-16 | exp decay rate — near zero makes EOS constant |
| `micro_liquid_density_b` | 1.0 | exponent in ω_l^b |
| `use_micro_liquid_density_for_pi` | true | controls ρ_prev in Π (rendered inert by step 3 fix) |
| `vdw_augmentation_prefactor` | see table above | K [J/kg] |
| `vdw_augmentation_decay_length` | 7.5e-7 m | λ = 0.75 nm |
| `micro_solid_density_reference` | 2780.0 kg/m³ | ρ_SR |
| `specific_surface` | 8e5 m²/kg | S_a (BET surface area) |

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

## Related memory files

- `~/.claude/projects/-Users-vinaykumar-git-ogs/memory/project_ms33_benchmark_prj_status.md`
- `~/.claude/projects/-Users-vinaykumar-git-ogs/memory/feedback_ogs_rm_sigma0_free_boundary.md`
- `~/.claude/projects/-Users-vinaykumar-git-ogs/memory/feedback_ogs_rm_pressure_tolerance.md`

## Related presentation

`/Users/vinaykumar/tex/cc2024/VK_B35_Pinion_May_2026/nagel_porosity_split.tex`

FIX slides: FIX-1 (n_l balance), FIX-2 (code changes), FIX-3 (impact),
FIX-4 (cross-coupling implicit), FIX-5 (EOS argument — step 3, 2026-05-21).
