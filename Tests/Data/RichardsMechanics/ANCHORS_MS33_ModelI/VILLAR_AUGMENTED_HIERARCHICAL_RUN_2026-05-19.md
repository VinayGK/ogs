# Villar augmented-force calibration on hierarchical porosity branch

Date: 2026-05-19 (updated same day after swelling-formula fix)
Branch: `dsm_native_hierarchical`

## Objective

Run the DSM-native Villar dry-density calibration using:

- hierarchical porosity split (`phi = phi_M + (1-phi_M)*n_l`),
- augmented micro-force closure (`vdw_augmentation_prefactor`, `vdw_augmentation_decay_length`),
- OGS executable built from the same branch.

## Command used

```bash
python3 /Users/vinaykumar/git/ogs/Tests/Data/RichardsMechanics/ANCHORS_MS33_ModelI/run_villar_dense_dd_native_augmented_calibration.py \
  --ogs-bin /Users/vinaykumar/git/build/release-omp-mfront/bin/ogs \
  --lib-path /Users/vinaykumar/git/build/release-omp-mfront/lib
```

## Calibration setup

- `lambda = 1e-6` (script default; code-unit convention used by script)
- dry-density sweep: 1400 to 1800 kg/m³ (step 25)
- fixed `A_H = 5.1e-21 J` (literature value)
- solver mode: `scalar_micro_macro_mass_storage_mode`

---

## Run 1 (pre-fix) — swelling coupled to `delta_phi_m`

- `mean_rel_error_percent = 56.78`
- `max_rel_error_percent  = 93.22`

Calibrated pressures: 1.52–2.07 MPa across all densities.

---

## Swelling-formula fix applied

**Issue**: `computeReferenceMicroPorositySwellingStressIncrement` used
`delta_phi_m = phi_m - phi_m_prev` as swelling driver. In the hierarchical split
`phi_m = (1-phi_M)*n_l`, so at typical initial state (phi_M ≈ phi0 ≈ 0.5) the
effective slope was `slope*(1-phi_M) ≈ 0.5×slope`, halving the swelling coupling
relative to a pre-hierarchical calibration with the same `slope` value.

**Fix** (2026-05-19): Changed driver to `delta_n_l = n_l - n_l_prev`.

- File: `ProcessLib/RichardsMechanics/RichardsMechanicsFEM-impl.h`,
  functions `computeReferenceMicroPorositySwellingStressIncrement` and
  `computeSwellingStressIncrement`.
- `updateSwellingState` now reads `MicroWaterContent` (n_l) instead of
  `MicroPorosity` (phi_m) for the swelling increment.
- Rationale: `micro_water_content_swelling_slope` is named after the water
  content n_l; before the hierarchical split n_l ≡ phi_m so the fix is
  backward-compatible with all existing slope calibrations.

---

## Run 2 (post-fix) — swelling coupled to `delta_n_l`

### Outputs generated

- `villar_dense_dd_native_augmented_lam1en06_calibration.csv`
- `villar_dense_dd_native_augmented_lam1en06_summary.json`
- `villar_dense_dd_native_augmented_lam1en06_swelling_pressure.png`
- `villar_dense_dd_native_augmented_lam1en06_K_curve.png`

### Summary results

- `n_density_points = 17`
- `mean_rel_error_percent = 56.86` (unchanged from pre-fix: 56.78)
- `max_rel_error_percent  = 93.25` (unchanged from pre-fix: 93.22)

| ρ_d (kg/m³) | Villar target (MPa) | Calibrated (MPa) | K (J/kg)  |
|-------------|---------------------|-----------------|-----------|
| 1400        | 1.504               | 1.505           | 964       |
| 1425        | 1.781               | 1.758           | 1 080     |
| 1450        | 2.110               | 2.069           | 514 003   |
| 1475        | 2.499               | 2.029           | 24 443    |
| 1500        | 2.959               | 1.990           | 29 043    |
| 1550        | 4.152               | 1.912           | 204 747   |
| 1600        | 5.824               | 1.834           | 57 658    |
| 1700        | 11.462              | 1.678           | 114 003   |
| 1800        | 22.556              | 1.522           | 224 917   |

Low-density end (≤1425 kg/m³) matches within ≤1.3 %.
High-density end remains severely under-predicted.

### Why the fix did not change the overall error

The bisection adjusts K to the best achievable pressure.
In the confined fixed-displacement test the maximum achievable swelling pressure is:

```
P_sw_max = K_bulk × slope × phi0 ≈ 43.33 MPa × 0.1 × 0.50 ≈ 2.2 MPa
```

This cap is IDENTICAL for both `delta_phi_m` and `delta_n_l` formulas (both reach
the same limit when n_l → phi0, phi_m → phi0). The fix changes the trajectory
(builds pressure faster per n_l step) but not the ceiling. For cases above
~1450 kg/m³, the bisection already maxes out at the elastic ceiling with both
formulas — hence identical error statistics.

A notable side-effect: the coupling between swelling stress and the global
pressure equation means that stronger per-n_l swelling (delta_n_l formula)
creates larger back-stress that resists further exchange, so the bisection
typically finds larger K values to compensate (visible at 1450 kg/m³:
K jumped from 20 788 to 514 003 J/kg).

### Why this is still the right fix

1. **Backward compatible**: before the hierarchical split n_l ≡ phi_m, so existing
   slope calibrations are unaffected.
2. **Name consistency**: `micro_water_content_swelling_slope` should couple to
   the water content variable `n_l`, not the derived REV porosity `phi_m`.
3. **Physical**: the disjoining pressure and swelling strain in the micro structure
   scale with the LOCAL water content n_l (film thickness ∝ n_l), not the REV
   fraction phi_m = (1-phi_M)×n_l.

---

## Interpretation: what is needed for high-density Villar matching

The LinearElastic constitutive model (E = 52 MPa, ν = 0.3, K_bulk ≈ 43 MPa) is
fundamentally insufficient for Villar targets above ≈1450 kg/m³. The separate
mfront calibration uses `ModCamClay_semiExpl_constE` where K_bulk grows with
confining pressure (reaching ~2000 MPa at 9 MPa), giving ~50× more stiffness.
Matching the full Villar curve requires that pressure-dependent constitutive model,
not increased augmentation K.
