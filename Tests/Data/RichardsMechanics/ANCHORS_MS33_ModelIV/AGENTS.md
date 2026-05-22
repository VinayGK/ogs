# AGENTS.md — ANCHORS MS33 DSM Pi-path benchmarks

**Scope:** EURAD-2 MS33 theoretical benchmarking suite, BGR DSM native
hierarchical vdW-augmented model. Models I / III / IV / VII.  
**Branch:** `dsm_native_hierarchical`  
**Binary:** `/Users/vinaykumar/git/build/release-omp-mfront/bin/ogs`  
(Do not use the `dsm_native-release` worktree binary — VTK ABI mismatch.)

---

## Roadmap

| Step | Description | Status |
|------|-------------|--------|
| 1 | REV-scale n_l mass balance: `φ_m·ρ_lR` storage, cross-coupling implicit in backward-Euler diff | ✅ `0d7a9edd64` (2026-05-20) |
| 2 | Thermodynamic swelling stress `σ_sw = −φ_m·Π`; K recalibration | ✅ `88d42c98fd` (2026-05-21) |
| 3 | Gibbs–Duhem fix: Pi-path uses `ρ_LR` (not `ρ_lR`); dead `use_micro_liquid_density_for_pi` flag removed | ✅ `c4888b6db4` / `ce9178fa96` (2026-05-21–22) |
| 4 | PRJ comment cleanup; presentation step marked done | ✅ 2026-05-21 |
| 5 | vdW base dimensional fix: `÷ρ_lR` added → J/kg; A = 2.2e-20 J (literature); K recalibrated | ✅ `0d579e8aeb` (2026-05-22) |
| 6 | DSM consistency hardening: viscosity guards, `use_micro_liquid_density_for_micro_pressure` default `true`, LE rerun | ✅ `66b782afa1` (2026-05-22) |
| 7 | Dead code removal: 3-arg `computeCompatibilityMicroHydraulicOutput` overload + flag | ✅ `4d47efff55` / `ce9178fa96` (2026-05-22) |
| 8 | Fix 5 failing DSM unit tests after step-5 vdW fix (reference functions + CSV baselines) | ✅ `3ac6b7de1f` (2026-05-22) — 13/13 pass |

---

## Physics summary

### Hierarchical porosity split

Total porosity `φ = φ_M + φ_m`. Micro porosity: `φ_m = (1−φ_M) · n_l` where
`n_l ∈ [0,1]` is the aggregate saturation of the micro pore space.

**Active solid fraction** (`computeActiveMicroSolidVolumeFraction`,
`MicroSolidVolumeFractionMode::CurrentPorositySplit`):
`n_S = 1 − φ_M` (micro-aggregate fraction, NOT pure solid).  
Commit `0d7a9edd64` corrected this from a formula that exceeded 1.

### REV-scale n_l mass balance

Storage term: `φ_m · ρ_lR = (1−φ_M)/(1−n_l) · n_l · ρ_lR` (REV scale, not
the aggregate `n_l · ρ_lR`). Cross-coupling term `−φ_m · φ̇_M / (1−φ_M)` is
implicit in the backward-Euler storage difference when φ_M is encoded at both
time levels — no separate code term needed.

### Swelling stress

```
σ_sw = −φ_m · Π  →  Δσ_sw = n_S · (n_l_prev·Π_prev − n_l·Π_curr) · I
```

Sign: n_l·Π increases during hydration → increment is negative → compressive,
consistent with swelling (tension-positive convention).

`Π = ρ_LR · K · exp(−ξ)` where `ξ = n_l / (λ · n_S · ρ_SR · S_a)`.

**Gibbs–Duhem consistency:** `Π` uses bulk `ρ_LR` (constant in MS33 setup),
not the micro EOS density `ρ_lR`. This ensures `dμ = (1/ρ_LR) · dΠ`.
The micro EOS density `ρ_lR` is retained in the mass-balance storage term only.

### vdW micro-potential (`computeVanDerWaalsMicroPotential`, `PotentialExchange.h`)

```
μ_lR_vdW = sign · A · Sa³ · n_S³ · ρ_SR³ / (6π · n_l³ · ρ_lR)   [J/kg]
μ_lR_aug = sign · K · exp(−ξ)                                      [J/kg]
μ_lR     = μ_lR_vdW + μ_lR_aug                                    [J/kg]
```

The `÷ρ_lR` in the vdW base was added in step 5 (`0d579e8aeb`); before that
it was in Pa (wrong). `dmu_lR_drho_lR = −μ_lR_vdW / ρ_lR` (non-zero after fix).

Exchange rate: `ρ̂_l = α_M_eff · (μ_LR − μ_lR)` where `μ_LR = p_L / ρ_LR`.  
Pressure check: `p_L_m = −ρ_lR · μ_lR` (impl.h lines 276, 1044) → μ_lR MUST
be in J/kg.

---

## Strict invariants

### 1. Augmentation is always additive — NEVER replace with assignment

```cpp
// CORRECT — augmentation adds on top of vdW base:
out.mu_lR += mu_aug;

// FORBIDDEN — discards vdW base entirely:
out.mu_lR = mu_aug;
```

Replacing `+=` with `=` silently discards the Hamaker base and makes `K` carry
all surface physics. This is undetectable from calibration alone (K absorbs the
mismatch) but invalidates the physical interpretation of both A and K.

### 2. A is a material constant — do NOT calibrate to swelling pressure

A is the Hamaker constant for clay–water–clay interaction.
**Set from literature only.** Current value: `2.2e-20 J` (Israelachvili & Adams
1978, mica–water–mica via SFA). Only K and λ are calibrated.

| System | A (J) | Source |
|--------|-------|--------|
| Mica–water–mica (SFA) | 2.2×10⁻²⁰ | Israelachvili & Adams (1978) |
| Montmorillonite–water (DLVO) | 1–5×10⁻²⁰ | Novich & Ring (1984) |
| Na-montmorillonite (FHH isotherm) | ~1.5×10⁻²⁰ | Cases et al. (1992) |
| Smectite (DFT/Lifshitz) | 1–3×10⁻²⁰ | Šolc et al. (2011) |

### 3. Adding new physical potential terms

Any new term (osmotic, double-layer, structural hydration) must:
1. Be in J/kg — check that `[term] / [kg/m³] = J/kg`
2. Be added as `out.mu_lR += new_term`
3. Be followed by K recalibration to Villar (K carries only the residual
   beyond vdW base + all additive terms)

---

## Parameters

### Current K values (after step-5 recalibration, all MS33 PRJ files)

| ρ_d (kg/m³) | K (J/kg) | p_sw,sim (MPa) | Villar (MPa) | Error |
|-------------|----------|----------------|--------------|-------|
| 1400 | 7 654.9 | 1.50376 | 1.50381 | −0.003% |
| 1600 | 29 984.9 | 5.82286 | 5.82407 | −0.021% |
| 1800 | 118 582.6 | 22.56322 | 22.55598 | +0.032% |

Villar target: `p_sw = exp(6.77 · ρ_d[g/cm³] − 9.07)` MPa.  
K_dd1600 = 29 984.9 J/kg propagated to Models III, IV, VII.

### Shared parameter values

| Parameter | Value | Role |
|-----------|-------|------|
| `hamaker_constant` | 2.2e-20 J | A — literature, do NOT calibrate |
| `vdw_augmentation_prefactor` | see table above | K [J/kg] — calibrated |
| `vdw_augmentation_decay_length` | 7.5e-7 m | λ = 0.75 nm |
| `micro_solid_density_reference` | 2780.0 kg/m³ | ρ_SR |
| `specific_surface` | 523 m²/kg | S_a (BET, MX-80 montmorillonite) |
| `micro_liquid_density_reference` | 100.0 kg/m³ | ρ_l0 (excess above bulk, near-constant EOS) |
| `micro_liquid_density_a` | 1e-16 | near-zero → ρ_lR ≈ ρ_LR + 100 ≈ 1100 kg/m³ |
| `use_micro_liquid_density_for_micro_pressure` | true | micro-pressure uses ρ_lR (default since step 6) |

### PRJ file locations (`potential_exchange` block)

```
Tests/Data/RichardsMechanics/
  ANCHORS_MS33_ModelI/
    ms33_modelI_dd1400.prj    ~line 52
    ms33_modelI_dd1600.prj    ~line 52
    ms33_modelI_dd1800.prj    ~line 52
  ANCHORS_MS33_ModelIII/
    ms33_modelIII_gap2mm.prj  ~line 69
  ANCHORS_MS33_ModelIV/
    ms33_modelIV_pellets.prj  ~line 73
  ANCHORS_MS33_ModelVII/
    ms33_modelVII_freeswelling.prj  ~line 80
```

---

## Calibration

```bash
# Calibrate all three densities (writes ms33_calibrate_K_results.txt):
python ms33_calibrate_K.py

# Verify existing K values without changing PRJ:
python ms33_calibrate_K.py --verify

# Single density:
python ms33_calibrate_K.py --density 1600
```

After calibration, propagate K_dd1600 to Models III/IV/VII:
```bash
K=<new_K_dd1600>
sed -i '' \
  "s|<vdw_augmentation_prefactor>[^<]*</vdw_augmentation_prefactor>|<vdw_augmentation_prefactor>${K}</vdw_augmentation_prefactor>|g" \
  ../ANCHORS_MS33_ModelIII/ms33_modelIII_gap2mm.prj \
  ../ANCHORS_MS33_ModelIV/ms33_modelIV_pellets.prj \
  ../ANCHORS_MS33_ModelVII/ms33_modelVII_freeswelling.prj
```

**Verification gate after any physics or K change:**
1. Villar errors < 0.1% for all three densities
2. All 6 MS33 models run with 0 rejected steps:
```bash
OGS=/Users/vinaykumar/git/build/release-omp-mfront/bin/ogs
$OGS -o . -l warn ms33_modelI_dd1400.prj
$OGS -o . -l warn ms33_modelI_dd1600.prj
$OGS -o . -l warn ms33_modelI_dd1800.prj
$OGS -o . -l warn ../ANCHORS_MS33_ModelIII/ms33_modelIII_gap2mm.prj
$OGS -o . -l warn ../ANCHORS_MS33_ModelIV/ms33_modelIV_pellets.prj
$OGS -o . -l warn ../ANCHORS_MS33_ModelVII/ms33_modelVII_freeswelling.prj
```
3. All 13 DSM unit tests pass:
```bash
cd /Users/vinaykumar/git/build/release-omp-mfront
./bin/testrunner --gtest_filter='*DSMMicroMacro*'
```

---

## Related files

- Memory: `~/.claude/projects/-Users-vinaykumar-git-ogs/memory/`
  - `project_ms33_benchmark_prj_status.md`
  - `feedback_ogs_rm_sigma0_free_boundary.md`
  - `feedback_ogs_rm_pressure_tolerance.md`
- Presentation: `/Users/vinaykumar/tex/cc2024/VK_B35_Pinion_May_2026/nagel_porosity_split.tex`
  (FIX slides 1–5 correspond to roadmap steps 1–3, 4, 5)
