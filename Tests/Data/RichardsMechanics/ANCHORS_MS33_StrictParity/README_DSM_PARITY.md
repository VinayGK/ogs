# ANCHORS MS33 — DSM Parity Benchmark

Two pairs of PRJ files for verifying that the native C++ and MFront bridge
DSM implementations produce identical physics.

## File pairs

### Pair 1: Non-DSM baseline (structural parity only)
| File | Branch binary | Purpose |
|---|---|---|
| `ms33_strict_parity_native.prj` | either | Standard RM with SaturationDependentSwelling — confirms solver/mesh are identical |
| `ms33_strict_parity_mfront.prj` | either | Identical to above (same content, different output prefix) |

### Pair 2: DSM parity (vdW potential exchange)
| File | Branch binary | Purpose |
|---|---|---|
| `ms33_dsm_parity_native.prj` | `dsm_native_hierarchical` | Native `<potential_exchange>` with hierarchical split |
| `ms33_dsm_parity_mfront.prj` | `dsm_mfront_hierarchical` | `RichardsMechanicsDSMMicroMacroBridge` with identical physics |

## How to run the DSM parity test

**Preferred — use the parity runner script** (runs both sims + reports MAE automatically):

```bash
# From repo root:
python3 scripts/run_dsm_parity.py

# Single suite only:
python3 scripts/run_dsm_parity.py --suite dsm_ms33

# Re-report without re-running (uses existing output in /tmp/dsm_parity_runs):
python3 scripts/run_dsm_parity.py --no-run
```

Output lands in `/tmp/dsm_parity_runs/dsm_ms33/{native,mfront}/`.
The script prints an MAE table split by early timesteps (ts≤50, 100→50 MPa suction)
and late timesteps (ts≥60, 50→0 MPa). Registered suites and field maps live in
`PARITY_SUITES` at the bottom of `scripts/run_dsm_parity.py`.

**Manual run (fallback):**

```bash
OGS_NATIVE=/path/to/native-release-omp-sharedcache/bin/ogs
OGS_MFRONT=/path/to/mfront-release-omp-sharedcache/bin/ogs
$OGS_NATIVE ms33_dsm_parity_native.prj -o /tmp/out/native
$OGS_MFRONT ms33_dsm_parity_mfront.prj -o /tmp/out/mfront
```

## Physics agreement (post-fix)

Both implementations use the same model after these fixes:

| Quantity | Native (`dsm_native_hierarchical`) | MFront bridge (`dsm_mfront_hierarchical`) |
|---|---|---|
| Porosity split | `phi_M = (phi-n_l)/(1-n_l)` (algebraic_split) | `phi_M = (phi0-n_l)/(1-n_l)` (phi≈phi0 small-strain) |
| Active nS | `1 - phi_M` (current_porosity_split mode) | `(1-phi0)/(1-n_l)` (from hierarchical split) |
| vdW potential | `A·Sa³·nS³·ρ_SR³ / (6π·n_l³·ρ_lR)` [J/kg] | same formula via omega formulation [J/kg] |
| Exchange | `α·(μ_LR - μ_lR)` where α in kg·s/m⁵ | same |
| Swelling | `ε̇_sw = β_sw · Φ̇_Micro` | same via SwellingSlope |

## Parameter values used

| Parameter | Value | Units | Notes |
|---|---|---|---|
| `phi0` | 0.41 | — | Initial total porosity |
| `hamaker_constant` / `HamakerConstant` | 2.2e-20 | J | Israelachvili & Adams 1978 |
| `specific_surface` / `SpecificSurface` | 5.23e5 | m²/kg | 523 m²/g FEBEX montmorillonite |
| `mass_exchange_coefficient` / `MassExchangeCoefficient` | 1e-10 | kg·s/m⁵ | J/kg domain |
| `micro_solid_density_reference` / `ReferenceDensitySolid` | 2780 | kg/m³ | |
| `n_l0` | 0.001 | — | Initial micro water content at 100 MPa suction |
| `micro_water_content_swelling_slope` / `SwellingSlope` | 0.1 | — | β_sw |

## Known residual differences

The mfront bridge uses the exponential micro-liquid EOS
(`rho_lR = rho_l0·exp(-A_eos·ω^B) + rho_LR_ref`) while the native
simple path uses `rho_LR_m = rho_LR` (constant macro density).
The parity PRJ minimises this difference by setting `rho_l0 = 1 kg/m³`
and `density_a = 50` so `rho_lR ≈ rho_LR_ref = 1000 kg/m³` in the bridge.
Residual difference from this source is < 0.1%.

## History of divergences (now fixed)

| Divergence | Before fix | After fix |
|---|---|---|
| Porosity split | mfront naive: `phi_m = n_l`, `phi_M = phi0 - n_l` | Both hierarchical |
| Active nS | mfront constant: `n_s = 1-phi0` | Both evolving |
| vdW units | mfront in Pa (missing `/rho_lR`) | Both J/kg |
