# Agent instructions — ANCHORS MS33 Model IV (and MS33 suite)

This folder simulates **ANCHORS MS33 Model IV** (`ms33_modelIV_pellets.prj`): a
confined bentonite block + pellet column under hydration, run with the BGR DSM
native hierarchical (double-structure) model in `RichardsMechanics`. It also
carries the **suite-wide** notes for the MS33 benchmark (Models I, III, IV, VII).
This is the authoritative copy; `ANCHORS_MS33_Model{I,III,VII}/AGENTS.md` are
short pointers back here.

Suite run log and scorecard: `../ANCHORS_MS33_RUN_SUMMARY_2026-05-20_anchors_dd.md`.
Status deck assembled from these artefacts:
`materialmodels/src/TPM/VK_SB_EURAD_DSM/` in the tex project tree.

---

## TODAY'S PHYSICS AUDIT — 2026-05-20

Guiding principle used in this audit:

1. Potentials are reduced by water uptake.
2. The microscale potential is the sum of vdW and augmentation terms.
3. Under the selected MS33 convention (`negative_attractive`), this combined micro potential should be most negative at the lowest water content.
4. Reduction in attraction drives microscale water uptake; this hydraulic-side evolution must emerge as swelling pressure (constrained) or swelling deformation (unconstrained) in mechanics.

### Process understanding (qualitative, for reviewers)

1. The implemented chain is: macro potential `mu_LR` and micro potential `mu_lR` define exchange `rho_l_hat = alpha_M * (mu_LR - mu_lR)`; exchange updates micro water content `n_l`; `n_l` then updates split porosities and swelling stress state.
2. In the MS33 setup (`micro_potential_convention = negative_attractive`), lower `n_l` gives stronger attraction (more negative `mu_lR`); hydration increases `n_l` and makes `mu_lR` less negative.
3. The model therefore represents hydration as a relaxation of attractive micro potential and a corresponding redistribution of liquid mass toward the micro domain.

### Source findings from today's audit

1. **Micro potential definition matches the intended sum form.** In `computeVanDerWaalsMicroPotential`, total `mu_lR` is vdW plus optional augmentation (`K * exp(-h/lambda)` with sign factor), so the combined potential is evaluated as one quantity in exchange.
2. **For the active MS33 settings, the combined potential is most negative at low water content.** With `negative_attractive` and positive augmentation prefactor, both terms are negative and decrease in magnitude with increasing `n_l` (film thickening).
3. **Micro-pressure compatibility variable needs targeted investigation.** Current output mapping uses `p_L_m = -rho_LR * mu_lR` (macro liquid density), while a physically stricter micro interpretation may require using `rho_lR`. This should be checked in theory notes and downstream use of `micro_pressure`.
4. **Micro-density EOS exists but is not yet physically calibrated.** The EOS path is implemented (`computeReducedMicroLiquidDensity`), but calibration against literature-backed microscale average density targets is still required; current MS33 values are numerically close to a near-trivial correction.
5. **Swelling-stress path needs hydraulic-mechanical consistency review.** With augmentation enabled, swelling increment follows the augmentation-pressure branch (`Pi`) and bypasses the slope-based branch. We must verify that disjoining-pressure evolution implied by changing micro potential is the same pressure that emerges as swelling stress in mechanics.
6. **Test coverage is good for sign/exchange algebra and local solve behavior, but not yet explicit for the full disjoining-pressure-to-swelling transfer identity under augmentation.**

### Required follow-up from this audit

1. **Micro pressure definition check:** decide and document whether compatibility output should use `rho_LR` or `rho_lR` in `p_L_m`.
2. **EOS calibration task:** fit `micro_liquid_density_reference`, `micro_liquid_density_a`, and `micro_liquid_density_b` to literature-consistent microscale density data and expected suction range.
3. **Swelling consistency task:** derive and verify the mapping from `Delta mu_lR` / disjoining pressure to `Delta sigma_sw` so hydraulic origin and mechanical response are thermodynamically and dimensionally consistent.
4. **Regression tests:** add focused tests for (a) combined-potential monotonicity at active MS33 parameters, (b) micro-pressure density choice, and (c) augmentation-mode swelling-stress consistency.

### Source investigation results and required actions (2026-05-20)

Source files examined:
- `ProcessLib/RichardsMechanics/ConstitutiveRelations/PotentialExchange.h` — vdW potential
- `ProcessLib/RichardsMechanics/RichardsMechanicsFEM-impl.h` — compatibility output, EOS, swelling
- `Tests/Data/RichardsMechanics/ANCHORS_MS33_ModelIV/ms33_modelIV_pellets.prj` — active parameters

---

#### Action 1 — Calibrate micro liquid EOS: REQUIRED, not optional

**Source:** `computeReducedMicroLiquidDensity` (`impl.h` line ~416). Current PRJ:
```xml
<micro_liquid_density_reference>1e-6</micro_liquid_density_reference>
<micro_liquid_density_a>1e-16</micro_liquid_density_a>
<micro_liquid_density_b>1.0</micro_liquid_density_b>
```
This gives `rho_lR ≈ rho_LR + 1e-6 kg/m³` — numerically zero correction.
**The EOS must not be left in this state.** The micro liquid density is a load-bearing
physics quantity: it governs mass exchange rates, the pressure conversion `p_L_m`,
and — once action 2 is implemented — the swelling stress magnitude directly.

**Target:** Interlayer water in Na-smectite is ~15 % denser than bulk at the dry
end of the suction range (Marry et al. 2002; Skipper et al. 1991; MD consensus).
The EOS must reproduce `rho_lR / rho_LR ≈ 1.15` at low `omega_l` (dry, high
suction) and converge toward 1.0 as `omega_l → ∞` (bulk-like wet state).

**How to calibrate:** The EOS form is
`rho_lR = rho_LR + rho_l0 * exp(-a * omega_l^b)`
where `omega_l = n_l * rho_lR / (nS * rho_SR)` (micro liquid mass fraction).
To get 15 % excess at small omega_l: set `rho_l0 ≈ 0.15 * rho_LR ≈ 150 kg/m³`,
`a` to control the crossover water content (e.g. `a ≈ 5`), `b ≈ 1.0`. Iterate:
the EOS is implicit in `rho_lR` and solved by the Newton loop already in the code.
Verify by evaluating `rho_lR` at the initial `n_l` for each dry density in
Model I and checking it falls in the 1050–1150 kg/m³ range.

**Gate:** After fitting, rerun Model I (all three dry densities). Swelling pressure
targets must still be met (< 1 % MAE). If they drift, adjust the swelling slope
(action 3) first, then the EOS.

---

#### Action 2 — Make rho_lR consistent everywhere: IMPLEMENTED via PRJ flags

**Status: code is in place; numerically active only after Action 1 (EOS calibration).**

Both density-consistency fixes are now controlled by PRJ flags and set in all
MS33 PRJ files.

**2a — p_L_m compatibility output — DONE:**
```xml
<use_micro_liquid_density_for_micro_pressure>true</use_micro_liquid_density_for_micro_pressure>
```
When `true`, `computeCompatibilityMicroHydraulicOutput` uses `micro_liquid_density.rho_lR`
instead of `rho_LR` for computing `p_L_m = −rho_lR * mu_lR` (`impl.h` line ~1013).

**2b — Pi computation — DONE:**
```xml
<use_micro_liquid_density_for_pi>true</use_micro_liquid_density_for_pi>
```
When `true`, `computeReferenceMicroPorositySwellingStressIncrement` uses
`rho_lR` (current) and `rho_lR_prev` (previous step) for `Pi_curr` and `Pi_prev`.
The function signature was updated to accept both (`impl.h` line 1469:
`double const rho_lR, double const rho_lR_prev`).

**Current numerical effect:** With the trivial EOS (`rho_l0 = 1e-6`,
`a = 1e-16`), `rho_lR ≈ rho_LR = 1000 kg/m³` to machine precision. The flags
have no measurable effect today. They become load-bearing once Action 1 calibrates
the EOS to give `rho_lR / rho_LR ≈ 1.15` at low water content.

**After Action 1:** rerun Model I and verify Villar MAE < 1 % still holds. If
swelling pressure shifts by more than ~15 % (the expected density correction),
refit K via `ms33_calibrate_K.py`.

---

#### Action 3 — Two valid mechanical paths; expose both, make selectable: IMPLEMENTED

**Status: DONE via `accumulate_swelling_contributions = true` in all PRJ files.**

**Background — Pi is not only hydraulic:**

In clay surface-force physics the disjoining pressure Pi between clay platelets
IS the swelling pressure. Two constitutive paths are available:

- **Pi-path** (augmentation-based): `delta_sigma_sw = nS * (Pi_curr - Pi_prev) * I`,
  `Pi = rho_lR * K * exp(-xi)`. Calibration through vdW/augmentation parameters
  (K, lambda, A_H). Physically motivated from surface-force theory.

- **Slope-path** (phenomenological): `delta_sigma_sw = -C_el * (slope * delta_n_l / 3) * I`.
  Calibration through slope coefficient fitted to swelling-pressure data.

**The historical bug (now resolved):**
The legacy code path (`accumulate_swelling_contributions = false`) had an early
return on `slope ≤ 0` that prevented the Pi block from running when `slope = 0`.

**Resolution — `accumulate_swelling_contributions = true`:**
The new code path (`impl.h` lines ~1527–1562) runs Pi and slope as independent
branches with no early return between them. `slope = 0` gives Pi-only;
`slope > 0` gives Pi + slope additively.

```xml
<accumulate_swelling_contributions>true</accumulate_swelling_contributions>
<micro_water_content_swelling_slope>0</micro_water_content_swelling_slope>
```
→ Pi-path only. This is the active configuration in all MS33 PRJ files.

**Selecting a path via PRJ parameters:**

| Intent | PRJ settings |
|---|---|
| Pi-path only | `vdw_augmentation_prefactor > 0`, `slope = 0`, `accumulate_swelling_contributions = true` |
| Slope-path only | `vdw_augmentation_prefactor > 0` (for exchange), `slope > 0`, `accumulate_swelling_contributions = true` |
| Both additive | `slope > 0`, augmentation enabled, `accumulate_swelling_contributions = true` |
| Legacy (pre-fix) | `accumulate_swelling_contributions = false` (do not use for new work) |

---

#### Action 4 — Regression tests: ADD AFTER ACTIONS 1–3

Tests to add in `Tests/` once actions 1–3 are complete:

1. **EOS physical range:** For the calibrated `(rho_l0, a, b)`, verify
   `rho_lR / rho_LR ∈ [1.05, 1.20]` at `omega_l = 0.01` and `→ 1.0` at
   `omega_l = 10`. Pure unit test on `computeReducedMicroLiquidDensity`.

2. **Slope branch always active:** With augmentation enabled, verify that a
   unit increase in `n_l` produces a non-zero `delta_sigma_sw` from the slope
   branch (would have been zero before this fix).

3. **Combined-potential monotonicity:** For active MS33 parameters
   (Hamaker = 5.1e-21, Sa = 523, K = 23423.8, λ = 1e-6), verify `mu_lR(n_l)`
   strictly monotone increasing over `n_l ∈ [1e-4, 0.5]` with
   `negative_attractive` convention.

4. **Density consistency round-trip:** After action 2, verify
   `p_L_m / rho_lR = -mu_lR` to machine precision at a known `(n_l, rho_lR)`.

---

**Sequencing:** Do actions in order 3 → 5 → 1 → 2 → 4.
Action 3 (slope fix) is a pure source change, independent of the EOS.
Action 5 (K calibration) must run after Action 3 — with the early return
present, sigma_sw = 0 and the calibration loop cannot converge.
Actions 1 and 2 are coupled — calibrate the EOS first, then propagate
`rho_lR` through the density-dependent expressions.
Rerun Model I as the gate after each action before proceeding.

---

#### Action 5 — vdW literature parameters and Pi-path K calibration: REQUIRED after Action 3

##### Background: what parameters govern the disjoining pressure

From `PotentialExchange.h` (the code comment block above
`computeVanDerWaalsMicroPotential`):

```
mu_lR_vdW = (A_H * Sa³ / 6π) * (nS * rho_SR)³ / n_l³
h = n_l / (nS * rho_SR * Sa)           [mean water film thickness, m]
mu_lR_aug = K * exp(−h / λ)
Pi = rho_lR * K * exp(−h / λ)          [disjoining pressure, Pa]
```

Parameters:
| Symbol | PRJ tag | Role |
|---|---|---|
| A_H | `hamaker_constant` | vdW dispersion amplitude [J] — from mineralogy |
| Sa  | `specific_surface` | effective specific surface [m²/kg code units] — sets film thickness h |
| ρ_SR | `micro_solid_density_reference` | solid grain density [kg/m³] — from mineralogy |
| K   | `vdw_augmentation_prefactor` | augmentation amplitude [J/kg] — **calibrated to swelling data** |
| λ   | `vdw_augmentation_decay_length` | characteristic film thickness [m code units] — from surface-force physics |

##### Literature values for fixed (physical) parameters

**Hamaker constant A_H — clay-water-clay geometry:**

| Source | A_H |
|---|---|
| Gregory (1981), J. Colloid Interface Sci. 83:138 | 5 × 10⁻²¹ J |
| Israelachvili (2011), *Intermolecular and Surface Forces*, 3rd ed., Table 13.1 | 1–5 × 10⁻²¹ J (clay in water) |
| Novich & Ring (1984), Clays and Clay Minerals 32:400 | 1.4–7.2 × 10⁻²¹ J (smectite, method-dependent) |
| Pashley & Israelachvili (1984) | ~5 × 10⁻²¹ J (mica-water-mica) |

**Literature consensus: A_H ≈ 5 × 10⁻²¹ J for Na-smectite in water.**
Current PRJ value `hamaker_constant = 5.1e-21 J` — within consensus. **No change.**

**Solid grain density ρ_SR:**

Na-montmorillonite grain density: 2740–2800 kg/m³ (Brigatti et al. 2006,
*Handbook of Clay Science*). Current `micro_solid_density_reference = 2780.0 kg/m³`. **No change.**

**Effective specific surface Sa — unit convention and physical range:**

The code formula `h = n_l / (nS * rho_SR * Sa)` requires Sa in m²/kg so that h
has units of metres. Verification at the Model I dd1400 initial state:
```
h₀ = n_l0 / (nS * rho_SR * Sa)
   = 1.0237e-3 / (0.5036 * 2780 * 523)
   = 1.0237e-3 / 732917
   = 1.40 × 10⁻⁹ m = 1.40 nm
```
A water film of 1.4 nm at 100 MPa suction (initial state) is physically correct
for 1–2 monolayers of interlayer water in Na-smectite (Saiyouri et al. 2000,
*Clay Minerals* 35:301; Marry et al. 2002, *J. Mol. Liq.* 96:455).

Literature BET/crystallographic specific surface for Na-montmorillonite:
700–800 m²/g = 700,000–800,000 m²/kg. This is **NOT** what `specific_surface`
in the code represents. The DSM model lumps the complex interlayer geometry
(stacked platelets, interlayer galleries) into a single "mean water film" model;
the effective Sa in that model (523 m²/kg) is a *phenomenological* parameter
calibrated so the film thickness h matches physical interlayer spacing at known
water content. **Do not replace with BET surface area.** Instead, verify that
h is physically reasonable (as above) and leave Sa = 523 unchanged.

**Augmentation decay length λ — physical surface-force decay lengths:**

The augmentation `K * exp(−h/λ)` represents the exponential short-range
repulsive contribution (electrical double-layer, EDL, or hydration forces).

Physically relevant decay lengths in clay systems:
| Force type | λ_physical |
|---|---|
| Hydration forces (water-molecule diameter scale) | 0.25–0.35 nm |
| EDL in distilled water (~1 μM ionic strength) | κ⁻¹ ≈ 300 nm |
| EDL in 1 mM NaCl | κ⁻¹ ≈ 9.6 nm |
| EDL in 10 mM NaCl | κ⁻¹ ≈ 3.0 nm |
| EDL in 100 mM NaCl (typical bentonite pore water) | κ⁻¹ ≈ 0.96 nm |

For bentonite barrier pore water (Na⁺ dominated, moderate salinity),
λ_physical ≈ 1–3 nm.

The PRJ comment `<!-- lambda=1e-6 code units = physical 1 nm decay length -->`
is physically incorrect if λ is in SI metres (1e-6 m = 1 μm ≠ 1 nm).
**This discrepancy must be investigated before λ is changed.**

Checking the effect of current λ = 1e-6 (code units) at dd1400 initial state:
```
xi₀ = h₀ / λ = 1.40e-9 m / 1e-6 = 1.40e-3   (if λ in metres)
exp(−xi₀) = 0.99861  →  augmentation barely decays at all
```
If λ were 1e-9 m (1 nm physical):
```
xi₀ = 1.40e-9 / 1e-9 = 1.40
exp(−xi₀) = 0.247   →  strong initial decay, exponential goes to ~0 at saturation
```

**The two physically consistent interpretations:**

A. λ = 1e-6 in code represents a **dimensionless ratio** or is in non-SI units
   (e.g., if distances in the model are in mm, 1e-6 mm = 1 pm — implausibly small).
   With xi₀ ≈ 1.4e-3, the augmentation changes by <0.14% over the entire
   hydration range; the exponential is essentially constant and acts as a simple
   additive offset.  K then sets the magnitude of a nearly-constant force, and
   must be very large to produce MPa-scale swelling.

B. λ = 1e-6 is a typo/legacy value, and the intended physical value is 1e-9 m (1 nm).
   With xi₀ = 1.4 and xi_final → large (as n_l grows to ~0.5 at saturation),
   exp(−xi) decays from 0.247 to ~0, and K is smaller but physically motivated.

**Resolution before calibration:**
Confirm the λ unit convention by reading
`ProcessLib/RichardsMechanics/ConstitutiveRelations/PotentialExchange.h`
comment and cross-checking the film-thickness formula units at the calibrated
initial state.  Update this AGENTS.md section with the confirmed interpretation.

For now, **treat λ as a fixed empirical parameter** (do not change 1e-6 until
the unit is confirmed) and calibrate K at this fixed λ.

##### Why K is the calibration target

K (the augmentation prefactor, J/kg) is not a pure mineralogical constant — it
absorbs the clay surface charge density, pore-fluid ionic composition, and
specific geometry of the interlayer.  These depend on the particular bentonite
and saturation conditions and are not directly measurable from mineral data
alone.  K is therefore the appropriate fit parameter.

λ is more fundamental (Debye length) and should come from surface-force
physics; once the unit convention is confirmed, λ should be set to its
physical value and K refitted.

##### Villar calibration target

MS33 specification uses the empirical dry density–swelling pressure relation
from Villar et al. (FEBEX bentonite data, multiple sources):

```
p_sw = exp(6.77 × ρ_d[g/cm³] − 9.07)  MPa
```

| ρ_d [kg/m³] | ρ_d [g/cm³] | p_sw_target [MPa] |
|---:|---:|---:|
| 1400 | 1.40 | 1.504 |
| 1600 | 1.60 | 5.824 |
| 1800 | 1.80 | 22.556 |

Note: these three points span a 15× range.  A single universal K cannot match
all three simultaneously unless the simulation dynamics (n_l evolution, Pi decay
path) create exactly the right density-dependent ratio.  Calibrate K separately
for each dry density; a global fit with both K and λ free is possible but
requires ~20–100 OGS runs (use the calibration script).

##### Calibration workflow

```
Step 1 — Apply Action 3 (remove early return in impl.h).
          Rebuild: cd /Users/vinaykumar/git/build/release-omp-mfront && ninja RichardsMechanics ogs
          All PRJ files must have micro_water_content_swelling_slope = 0. ✓ (done 2026-05-20)

Step 2 — Verify that existing K values still give < 1% MAE after the bug fix.
          The K values in the PRJ files were calibrated with the Pi-path active
          (slope=0.1 with augmentation enabled caused Pi early return → K was
          fitted to Pi-path output). After fixing the bug, the same K should
          reproduce the same sigma_sw.

          cd ANCHORS_MS33_ModelIV
          python ms33_calibrate_K.py --verify

          Expected: all three densities < 1% error, mean < 0.8%.
          If PASSED: skip Steps 3–4.

Step 3 — If verification fails, calibrate K for each density (bisection):
          python ms33_calibrate_K.py
          (or python ms33_calibrate_K.py --density 1400  for a single density)
          Results written to ms33_calibrate_K_results.txt.

Step 4 — Propagate calibrated K to all PRJ files.
          Models III, IV, VII all use ρ_d = 1600 kg/m³ reference → K(1600).

          Files to update:
            ANCHORS_MS33_ModelI/ms33_modelI_dd1400.prj  → K(1400)
            ANCHORS_MS33_ModelI/ms33_modelI_dd1600.prj  → K(1600)
            ANCHORS_MS33_ModelI/ms33_modelI_dd1800.prj  → K(1800)
            ANCHORS_MS33_ModelIII/ms33_modelIII_gap2mm.prj → K(1600)
            ANCHORS_MS33_ModelIV/ms33_modelIV_pellets.prj  → K(1600)
            ANCHORS_MS33_ModelVII/ms33_modelVII_freeswelling.prj → K(1600)

Step 5 — Rerun all ANCHORS models and verify gates (see below).
```

##### Simulated swelling pressure extraction

In Model I (constant-volume isotropic), the swelling pressure at t = 200 d is
extracted as the spatially averaged mean compressive stress:

```python
# In the VTU, sigma is a Kelvin vector (2D axisymmetric):
# [sigma_xx, sigma_yy, sigma_zz, sigma_xy*sqrt(2)]
mean_stress = (sigma_xx + sigma_yy + sigma_zz) / 3
p_sw = -mean(mean_stress over all nodes)   # compressive → positive
```

The calibration script (`ms33_calibrate_K.py` in this directory) automates
this extraction using the vtk Python package.

##### Rerun commands for all ANCHORS models

After calibrating K and propagating to all PRJ files:

```bash
BUILD=/Users/vinaykumar/git/build/release-omp-mfront
OGS=$BUILD/bin/ogs
DATA=/Users/vinaykumar/git/ogs/Tests/Data/RichardsMechanics

# Rebuild (only if impl.h was changed for Action 3):
cd $BUILD && ninja RichardsMechanics ogs

# Model I — three dry densities
for DD in dd1400 dd1600 dd1800; do
    PRJ=$DATA/ANCHORS_MS33_ModelI/ms33_modelI_${DD}.prj
    OUTDIR=$DATA/ANCHORS_MS33_ModelI
    LOG=$OUTDIR/ms33_modelI_${DD}_run_$(date +%Y%m%d).log
    $OGS -o $OUTDIR -l warn $PRJ > $LOG 2>&1 &
done
wait

# Model III
PRJ=$DATA/ANCHORS_MS33_ModelIII/ms33_modelIII_gap2mm.prj
$OGS -o $DATA/ANCHORS_MS33_ModelIII -l warn $PRJ \
    > $DATA/ANCHORS_MS33_ModelIII/ms33_modelIII_run_$(date +%Y%m%d).log 2>&1

# Model IV
PRJ=$DATA/ANCHORS_MS33_ModelIV/ms33_modelIV_pellets.prj
$OGS -o $DATA/ANCHORS_MS33_ModelIV -l warn $PRJ \
    > $DATA/ANCHORS_MS33_ModelIV/ms33_modelIV_run_$(date +%Y%m%d).log 2>&1

# Model VII
PRJ=$DATA/ANCHORS_MS33_ModelVII/ms33_modelVII_freeswelling.prj
$OGS -o $DATA/ANCHORS_MS33_ModelVII -l warn $PRJ \
    > $DATA/ANCHORS_MS33_ModelVII/ms33_modelVII_run_$(date +%Y%m%d).log 2>&1
```

##### Verification gates after rerun

| Gate | Criterion | Check command |
|---|---|---|
| **G-K**: Model I swelling pressure | MAE < 1% at all three densities | `python ms33_calibrate_K.py --verify` |
| **G-phi**: transport_porosity ≥ 0 | No negative values at any node/step | `python check_transport_porosity_nonnegative.py` |
| **G-conv**: All models complete | 0 rejected steps in logs | `grep -i "rejected\|failed" *.log` |
| **G-phi2**: micro_porosity ≤ porosity | 0 violations node-wise | VTK check (see Model IV defect section) |

Gate G-K is the primary calibration gate. All four must pass before models
move from PARTIAL to PASS and before PDFs are regenerated.

##### PRJ documentation standard for K

After calibration, annotate the K value in every PRJ with the target it was
fitted to and the calibration method:

```xml
<!-- K calibrated to Villar swelling pressure at rho_d=1400 kg/m³ via Pi-path.
     Villar target: exp(6.77*1.4 - 9.07) = 1.504 MPa.
     Calibration script: ANCHORS_MS33_ModelIV/ms33_calibrate_K.py
     A_H=5.1e-21 J (literature), Sa=523 (code units), lambda=1e-6 (code units, see AGENTS.md).
     Recalibrate if lambda unit convention is resolved and lambda is updated. -->
<vdw_augmentation_prefactor>4981.81</vdw_augmentation_prefactor>
```

---

## PRIMARY OPEN DEFECT — negative macro porosity in Model IV

### Problem statement

In the DSM the pore space is split into two structural levels. OGS enforces, at
every output node, exactly:

```
porosity (total)  =  micro_porosity  +  transport_porosity
```

where `transport_porosity` is the **macro** (inter-aggregate, advective)
porosity = `phi_macro`, and `micro_porosity` is the intra-aggregate porosity.

In the Model IV final state (`ms33_modelIV_pellets_ts_302_t_17280000.000000.vtu`,
t = 200 d) `transport_porosity` is **negative** in the most-confined nodes:

| Field (final step) | min | max |
|---|---:|---:|
| `porosity` (total)           | 0.290 | 0.719 |
| `micro_porosity`             | 0.001 | 0.719 |
| `transport_porosity` (macro) | **-0.0134** | 0.294 |

Verified: `porosity - micro_porosity - transport_porosity = 0` exactly at all
90 nodes. Negative macro porosity occurs at 6 of 90 nodes; at those nodes
`micro_porosity = 0.407` while total `porosity` has fallen to ~0.394, so the
macro level is forced to `0.394 - 0.407 = -0.013`. Porosity must lie in [0, 1];
a negative macropore volume is non-physical.

### Root cause

1. **Physical driver.** Under confinement the clay aggregates swell on
   hydration — `micro_porosity` (and `micro_water_content`) rises — faster than
   the bulk pore volume can yield. The macropore network closes completely and
   the model then keeps driving it past closure. This is worst in Model IV
   because the pellet/block density contrast is large
   (`dry_density_solid` spans ~782..1973 kg/m3), so the loose pellet zone swells
   hardest into the densest, most-confined macro voids.

2. **Missing constitutive coupling.** The micro-porosity / micro-macro exchange
   update (`micro_exchange_source`) and the total-porosity update are not
   mutually constrained, so nothing stops `micro_porosity` from exceeding total
   `porosity`. When `phi_macro -> 0` the micro structure should no longer be
   able to draw volume from a macropore network that no longer exists.

3. **The clamp is not governing the reported field.** `ms33_modelIV_pellets.prj`
   declares `<minimal_porosity>0</minimal_porosity>` on the
   `TransportPorosityFromMassBalance` property for both materials (clay:
   lines ~122-127; pellet: lines ~197-202). The reported `transport_porosity`
   still goes negative. Because the decomposition above holds *exactly*, the
   output field is the derived difference `porosity - micro_porosity`, i.e. the
   `TransportPorosityFromMassBalance` clamp is bypassed or applied before the
   micro-exchange contribution. Confirm this in the OGS source.

---

### What needs to be done

Work the three layers in order; each has a verification gate.

1. **OGS source — locate and fix the production of `transport_porosity`.**
   In `ProcessLib/RichardsMechanics` (the DSM / micro-porosity path), find where
   `transport_porosity` is assembled for output. Either (a) make
   `TransportPorosityFromMassBalance` with `minimal_porosity` actually bound the
   reported field, or (b) if it is the derived difference, enforce
   `phi_macro = max(porosity - micro_porosity, 0)` consistently for the field,
   the permeability input (`KozenyCarman` uses it), and any transport term.
   *Gate:* a rebuilt OGS rerun of `ms33_modelIV_pellets.prj` yields
   `transport_porosity >= 0` at every node and step.

#### Investigation log — OGS source fix (completed 2026-05-20)

**Key source file:** `ProcessLib/RichardsMechanics/RichardsMechanicsFEM-impl.h`

**Root cause identified:**

`updateSwellingStressAndVolumetricStrain` (line ~1560) is the injection point.
When `potential_exchange_enabled && saturation_micro` is present, it does:
```cpp
phi_M.phi = phi_M_prev->phi;   // copies PREVIOUS step value — NEGATIVE if prev was negative
```
`phi_M` is a direct reference to `current_states_[ip].TransportPorosityData`.
The negative propagates: once one step produces a negative `phi_M`, that value is
copied into `prev_states_` via `postTimestepConcrete`, and the next step injects
it again via this line — creating a self-sustaining negative that never clears.

**Fixes applied (committed in `179559bbd4`; plus output-layer backstop in current tree):**

1. **Source fix — stop injection at origin** (approved, line 1644):
   ```cpp
   // RichardsMechanicsFEM-impl.h, inside updateSwellingStressAndVolumetricStrain
   phi_M.phi = phi_M_prev->phi;
   // Prevent propagation of a non-physical negative macro porosity
   // from previous-step state into current assembly/output.
   phi_M.phi = std::max(0.0, phi_M.phi);
   ```
   This stops the negative from entering `current_states_[ip]` at all. Safe
   because `phi_M >= 0` is a physical invariant; the negative in `phi_M_prev`
   was itself an artefact.

2. **Gate 1 clamp in `assembleWithJacobianEvalConstitutiveSetting`** (line ~2770):
   After all porosity-split updates, clamps `phi_m <= phi_total` and sets
   `phi_M = phi_total - phi_m >= 0` in `state_current`. Prevents any residual
   negative from reaching the Jacobian assembly.

3. **Gate 1 clamp in `computeSecondaryVariableConcrete`** (line ~3681):
   Same clamp applied to `current_states_[ip]` in the secondary-variable path,
   ensuring the value written via the extrapolation is non-negative.

4. **Output-level clamp in `AddProcessDataToMesh.cpp`**:
   After nodal extrapolation, any `transport_porosity` value < 0 is clamped to 0
   before being written to the VTU. This is a backstop against extrapolation
   overshoot at material interfaces and is **not a substitute** for fixes 1–3:
   fixes 1–3 are needed to keep `current_states_[ip]` clean so the negative does
   not affect assembly physics or propagate into `prev_states_`.

**Why fixes 2–4 alone were insufficient (historical note):**

Fixes 2 and 3 were applied and tested first; VTU still showed negatives. The
reason: `updateSwellingStressAndVolumetricStrain` re-injects the negative from
`phi_M_prev` every timestep, and because it is called *before* the Gate 1 clamps
in both assembly and secondary-variable paths, the negative was being set, then
the downstream overwrites (`updateMicroscaleHydraulicState`,
`updatePorositySplitState`) were somehow not overriding it in the output path
(root cause of that failure not fully traced). Fix 1 at line 1644 is the correct
and sufficient intervention: it prevents injection so no downstream cleanup is
needed. Fixes 2–4 remain as defence-in-depth.

**GATE 1 PASSED — 2026-05-20 17:02.** No `NEG` lines. All 4 timesteps clean:

| VTU | transport_porosity min | max |
|---|---:|---:|
| ts_0   t=0 d   | 0.424108 | 0.676061 |
| ts_76  t=20 d  | 0.340831 | 0.707953 |
| ts_178 t=100 d | 0.000000 | 0.316757 |
| ts_302 t=200 d | 0.000000 | 0.294324 |

Min reaches 0.000 (macropores fully closed) but never negative. Fix confirmed.

**GATE 2 PASSED (output invariants) — 2026-05-20.**
After rerun with output-field reconciliation in
`ProcessLib/Output/AddProcessDataToMesh.cpp`, all 4 VTUs satisfy
`micro_porosity <= porosity` node-wise (0 violations) and keep
`transport_porosity >= 0` (0 negative nodes).

Remaining for full benchmark PASS: implement constitutive-level saturation
(`micro_porosity <= porosity` by state evolution, not only output reconciliation),
regenerate PDFs, update run summary.

```bash
cd /Users/vinaykumar/git/build/release-omp-mfront
ninja RichardsMechanics ogs

OUTDIR=/Users/vinaykumar/git/ogs/Tests/Data/RichardsMechanics/ANCHORS_MS33_ModelIV
PRJ=$OUTDIR/ms33_modelIV_pellets.prj
OGS=/Users/vinaykumar/git/build/release-omp-mfront/bin/ogs
$OGS -o $OUTDIR -l warn $PRJ > $OUTDIR/rerun_modelIV.log 2>&1

python3 -c "
import vtk, glob
files = sorted(glob.glob('$OUTDIR/ms33_modelIV_pellets_ts_*_t_*.vtu'))
for f in files:
    r = vtk.vtkXMLUnstructuredGridReader()
    r.SetFileName(f); r.Update()
    a = r.GetOutput().GetPointData().GetArray('transport_porosity')
    if a:
        vals = [a.GetValue(i) for i in range(a.GetNumberOfTuples())]
        if min(vals) < 0:
            print('NEG', f[-60:], 'min=%.6f' % min(vals))
print('done')
"
```
*Gate:* no `NEG` lines — `transport_porosity >= 0` at every node and timestep.

---

2. **Constitutive — cap micro swelling against available macro volume.**
   Bound the micro-macro exchange / micro-porosity growth so
   `micro_porosity <= porosity` always holds: as `phi_macro -> 0` the
   micro->macro water exchange and micro swelling must saturate. This is the
   physically correct fix; the clamp in step 1 only prevents the symptom.
   *Gate:* in the rerun, `micro_porosity <= porosity` node-wise, and the
   macro-porosity trajectory reaches but does not cross zero.

3. **Postprocess — `ms33_postprocess_modelIV.py`.**
   - The dry-density fallback is wrong: the docstring/`zone_stats` use
     `rho_d = RHO_S * (1 - phi_M)` with `phi_M = transport_porosity` (macro
     only). Dry density is solid mass per bulk volume and must use **total**
     porosity: `rho_d = RHO_S * (1 - porosity)`. The script currently prefers
     the native `dry_density_solid` field (correct) and only falls back to this
     formula — fix the fallback so it is right if ever used.
   - Add a diagnostic: count and print nodes/steps where
     `transport_porosity < 0`, and clamp negatives before plotting so a clamp
     artefact is never silently presented as a result.
   - **Status (2026-05-20): DONE.** Both fixes applied to
     `ms33_postprocess_modelIV.py`. Awaiting clean VTU from source fix before
     regenerating PDFs.

After all three gates pass: rerun the suite, regenerate
`ms33_modelIV_dry_density.pdf` / `ms33_modelIV_mean_stress.pdf`, update
`../ANCHORS_MS33_RUN_SUMMARY_2026-05-20_anchors_dd.md`, and only then may
Model IV move from PARTIAL to PASS.

---

## Binding rule — record output and result paths

Whenever a figure, table, or reported number is produced from a simulation,
record the **full provenance chain** so reviewers and the status deck can trace
it: the result artefact, the postprocess script, and the `.prj` / `.pvd` /
`.vtu` / `.log` it came from. Paths are repo-relative, rooted at
`Tests/Data/RichardsMechanics/`.

- Put the chain in the run summary
  (`../ANCHORS_MS33_RUN_SUMMARY_2026-05-20_anchors_dd.md`) and in the docstring
  of each `ms33_postprocess_model*.py`.
- The status deck has the matching rule: **every slide that plots, displays, or
  pulls data from simulation results must carry a visible source line** naming
  these artefacts. When you regenerate a figure or value, the deck's source line
  must be updated in the same change. See
  `materialmodels/src/TPM/VK_SB_EURAD_DSM/AGENTS.md`.

Model IV artefacts: `ms33_modelIV_pellets.prj` -> `ms33_modelIV_pellets.pvd`
(+ `_ts_*.vtu`); postprocess `ms33_postprocess_modelIV.py` ->
`ms33_modelIV_mean_stress.pdf`, `ms33_modelIV_dry_density.pdf`; run logs
`anchors_modelIV_full.log`, `rerun_modelIV.log`.

---

## MS33 suite — other open items

- **Model V** — no `.prj` exists anywhere in the tree. The benchmark suite is
  incomplete until a Model V case is added.
- **Model VII** (`ANCHORS_MS33_ModelVII/`) — separate open defect: simulated
  void ratio during loading/unloading is `e ~ 2.81..2.96`, far above the
  benchmark reference band `0.4..1.2`. The stress-path interpretation is already
  fixed (total stress `sigma_total = sigma_effective + chi*p_L`, `chi = S`,
  matches targets within ~0.01 MPa); the residual is genuine model behaviour,
  not postprocessing. See `ANCHORS_MS33_ModelVII/AGENTS.md`.
- **Models I and III** — PASS. Model I dry-density calibration mean abs. error
  0.785 %; Model III outputs finite and bounded. Keep them as regression
  baselines: do not let the macroporosity fix change Model I/III reference VTUs.

---

## Timestamped problem log

### 2026-05-20 23:05:03 CEST — open problems after latest OGS fixes

1. **RichardsMechanics ctest gate still has 6 failing checks** (build tree:
   `build/release-omp-mfront/ProcessLib/RichardsMechanics`):
   - `ogs-RichardsMechanics/double_porosity_swelling`
   - `ogs-RichardsMechanics/double_porosity_swelling-omp`
   - `ogs-RichardsMechanics/DoubleStructureBenchmark/double_porosity_swelling_RM`
   - `ogs-RichardsMechanics/DoubleStructureBenchmark/double_porosity_swelling_RM-omp`
   - `ogs-RichardsMechanics_beacon_1a01_dsm_micromacro_inflow_reference-time-vtkdiff`
   - `ogs-RichardsMechanics_beacon_1a01_dsm_micromacro_inflow_reference-time-omp-vtkdiff`
   **Cause:** `ogs` binary is pre-hierarchical-split (built 17:19); these tests
   compare against VTU reference files. Requires OGS rebuild + VTU regeneration.

2. **DSM unit baseline drift — FIXED (2026-05-20)**:
   The following four unit tests were failing because `computeTransportPorosityUpdate`
   was changed to enforce the hierarchical split `phi_M = (phi − n_l)/(1 − n_l)`,
   while test expected values and CSV baselines still used the old flat split
   (`phi_m ≈ n_l`). Fixed by updating expected values and regenerating CSVs:
   - `RichardsMechanics.DSMMicroMacroTransportPorositySplitRecomposesTotalPorosity` ✓
   - `RichardsMechanics.DSMMicroMacroAdditiveMacroPorosityRateUpdate` ✓
   - `RichardsMechanics.DSMMicroMacroOverlapTransferBaselineHistory` ✓
   - `RichardsMechanics.DSMMicroMacroStrainCoupledOverlapBaselineHistory` ✓
   All 24 RichardsMechanics unit tests pass as of `bin/testrunner` rebuild at 22:29.

2c. **rho_lR initialisation bug — FIXED (2026-05-21):**
   Root cause: `rho_lR_initial` was set to `micro_liquid_density_reference = 1e-6`
   (a trivial EOS placeholder) and written to both `rho_lR` and `rho_lR_prev` at
   init time (impl.h ~line 2015). In step 1 the exchange solve updates `rho_lR`
   to ≈ `rho_LR` ≈ 1000 kg/m³, while `rho_lR_prev` stays at 1e-6. With
   `use_micro_liquid_density_for_pi = true`, this gives:
   `Pi_prev = 1e-6 * K * exp(-ξ_prev) ≈ 0` vs `Pi_curr = 1000 * K * exp(-ξ_curr)`
   → ~10⁶× tensile sigma_sw spike in step 1 that dominates the accumulation.
   Result: p_sw = −0.99, −7.6, −45.9 MPa (NEGATIVE) vs targets 1.50, 5.82, 22.6 MPa.
   **Fix:** after `rho_LR_initial` (actual macro density) is known during
   `initializeConcrete`, call `computeActiveMicroLiquidDensity(n_l_initial, rho_LR_initial, ...)`
   and overwrite `*rho_lR` and `**rho_lR_prev` with the corrected value (≈ rho_LR ≈ 1000).
   Location: `RichardsMechanicsFEM-impl.h`, added after `computeTransportPorosityUpdate`
   block inside `if (isPotentialExchangeEnabled(...))` in `initializeConcrete`.
   **Post-fix verification:** `ms33_calibrate_K.py --verify` → MAE = 0.78 %,
   all three densities GATE PASSED < 1 % (K values unchanged: 4981.81, 23423.8, 105429.7).

2b. **RichardsMechanics ctest gate — FIXED (2026-05-21):**
   After rebuilding `ogs` and regenerating reference VTUs with the new binary,
   all 6 previously failing tests now pass:
   - `ogs-RichardsMechanics/double_porosity_swelling` ✓
   - `ogs-RichardsMechanics/double_porosity_swelling-omp` ✓
   - `ogs-RichardsMechanics/DoubleStructureBenchmark/double_porosity_swelling_RM` ✓
   - `ogs-RichardsMechanics/DoubleStructureBenchmark/double_porosity_swelling_RM-omp` ✓
   - `ogs-RichardsMechanics_beacon_1a01_dsm_micromacro_inflow_reference-time-vtkdiff` ✓
   - `ogs-RichardsMechanics_beacon_1a01_dsm_micromacro_inflow_reference-time-omp-vtkdiff` ✓
   Cause: reference VTUs were produced before the hierarchical split and
   `computeCompatibilityMicroHydraulicOutput` (3→4 param) changes. Updated:
   `Tests/Data/RichardsMechanics/double_porosity_swelling_t_{50000,100000}.vtu`,
   all `DoubleStructureBenchmark/double_porosity_swelling_t_*.vtu`,
   `beacon_1a01_dsm_micromacro_inflow_reference_t_100000.vtu`.
   Full suite: 208/209 pass. Remaining failure:
   `ThermoRichardsMechanics/.../bentonite_column-LARGE-omp` — pre-existing OMP
   floating-point non-determinism at picoscale tolerance (5.7e-12), unrelated to
   DSM changes; serial version passes.

3. **Model IV submission readiness remains blocked** by constitutive-level
   saturation coupling (`micro_porosity <= porosity` must hold by evolution,
   not only by output reconciliation) and by full rerun/regenerated PDFs under
   the updated physics switches.

### 2026-05-20 — Full clean-rebuild and ground-up validation procedure

Run this sequence whenever C++ source has been edited. Also required now because
the `ogs` binary predates the hierarchical-split and Pi-path changes.

**Step 0 — Remove stale output artefacts:**
```bash
DATA=/Users/vinaykumar/git/ogs/Tests/Data/RichardsMechanics
for DIR in ANCHORS_MS33_ModelI ANCHORS_MS33_ModelIII ANCHORS_MS33_ModelIV ANCHORS_MS33_ModelVII; do
    rm -f $DATA/$DIR/*.pvd $DATA/$DIR/*.vtu $DATA/$DIR/*.log
done
```

**Step 1 — Rebuild OGS binary:**
```bash
cd /Users/vinaykumar/git/build/release-omp-mfront
ninja RichardsMechanics ogs
```

**Step 2 — Rebuild testrunner and run unit tests:**
```bash
ninja testrunner
bin/testrunner --gtest_filter="RichardsMechanics.*"
```
All unit tests must pass before proceeding. Fix any failures first.

**Step 3 — Full RichardsMechanics ctest:**
```bash
ctest -j 18 --output-on-failure -R "RichardsMechanics"
```
If vtkdiff tests fail, check whether reference VTUs need regeneration with the
new binary (the old VTUs were produced with the pre-hierarchical-split binary).

**Step 4 — Verify K calibration:**
```bash
cd /Users/vinaykumar/git/ogs/Tests/Data/RichardsMechanics/ANCHORS_MS33_ModelIV
python ms33_calibrate_K.py --verify
# Expected: MAE < 1 % for dd1400, dd1600, dd1800
# If FAIL: python ms33_calibrate_K.py   (bisection; propagate new K(dd1600) to
#           Models III, IV, VII PRJs if the value changes)
```

**Step 5 — Ground-up MS33 model runs:**
```bash
OGS=/Users/vinaykumar/git/build/release-omp-mfront/bin/ogs
DATA=/Users/vinaykumar/git/ogs/Tests/Data/RichardsMechanics

# Model I — three dry densities (parallel)
for DD in dd1400 dd1600 dd1800; do
    $OGS -o $DATA/ANCHORS_MS33_ModelI -l warn \
         $DATA/ANCHORS_MS33_ModelI/ms33_modelI_${DD}.prj \
         > $DATA/ANCHORS_MS33_ModelI/anchors_${DD}_full.log 2>&1 &
done; wait

# Model III
$OGS -o $DATA/ANCHORS_MS33_ModelIII -l warn \
     $DATA/ANCHORS_MS33_ModelIII/ms33_modelIII_gap2mm.prj \
     > $DATA/ANCHORS_MS33_ModelIII/ms33_modelIII_run_$(date +%Y%m%d).log 2>&1

# Model IV
$OGS -o $DATA/ANCHORS_MS33_ModelIV -l warn \
     $DATA/ANCHORS_MS33_ModelIV/ms33_modelIV_pellets.prj \
     > $DATA/ANCHORS_MS33_ModelIV/ms33_modelIV_run_$(date +%Y%m%d).log 2>&1

# Model VII
$OGS -o $DATA/ANCHORS_MS33_ModelVII -l warn \
     $DATA/ANCHORS_MS33_ModelVII/ms33_modelVII_freeswelling.prj \
     > $DATA/ANCHORS_MS33_ModelVII/ms33_modelVII_run_$(date +%Y%m%d).log 2>&1
```

**Step 6 — Post-run checks:**
- `grep "rejected" $DATA/ANCHORS_MS33_Model*/*.log` → must show 0 rejected steps for all 6 runs.
- `transport_porosity >= 0` at every node/step (VTK inspection or `ms33_calibrate_K.py --verify`).
- Regenerate PDFs: `python ms33_postprocess.py` (Model IV), `ms33_postprocess_modelIII.py` (III), `ms33_postprocess_modelVII.py` (VII).
- Add a new timestamped entry to this problem log recording the outcome.
