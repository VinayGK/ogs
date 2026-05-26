# DSM MFront Hierarchical — Patch Recipe and Engineering Log

## Objective
`dsm_mfront_hierarchical` implements the hierarchical DSM porosity split and J/kg
vdW potential via `RichardsMechanicsDSMMicroMacroBridge.mfront`, so that the physics
match the native `dsm_native_hierarchical` C++ implementation.

## Branch set
- `master`
- `dsm_native` — naive/additive porosity, naive vdW (Pa units)
- `dsm_native_hierarchical` — hierarchical split, J/kg vdW (commit 0d579e8aeb)
- `dsm_mfront` — mfront bridge, naive porosity, Pa-units vdW
- `dsm_mfront_hierarchical` — **this branch**: hierarchical split, J/kg vdW in bridge

## PATCH_RECIPE MAINTENANCE RULE
Any change to code, tests, or benchmarks that affects DSM constitutive physics
must be reflected here as a new section with: what changed, why, diff summary,
and validation outcome.

---

## Physics divergences identified and fixed (2026-05-25, commit de589d8fc0)

### Root-cause analysis
Comparison of `RichardsMechanicsDSMMicroMacroBridge.mfront` (this branch) against
`ProcessLib/RichardsMechanics/ConstitutiveRelations/PotentialExchange.h`
(dsm_native_hierarchical) revealed three independent divergences:

| # | Location | Was (WRONG) | Is (CORRECT) |
|---|---|---|---|
| 1 | `phi_m_prev`, `phi_m_trial`, `phi_M_trial`, `@UpdateAuxiliaryStateVariables` | `phi_m = n_l`, `phi_M = phi0 - n_l` (naive/additive) | `phi_M = (phi0-n_l)/(1-n_l)`, `phi_m = n_l*(1-phi_M)` (hierarchical) |
| 2 | `omega_from_state` lambda and `@UpdateAuxiliaryStateVariables` | `n_s = 1.0 - phi0` (constant initial solid fraction) | `active_nS = (1-phi0)/(1-n_l)` (evolving from hierarchical split) |
| 3 | `mu_micro_from_state` lambda and `@UpdateAuxiliaryStateVariables` | `(A * rho_lR^3 / 6pi) * (Sa/omega)^3` — **Pa** (wrong units) | `A / (6pi * rho_lR) * (Sa/omega)^3` — **J/kg** (correct, matches native commit 0d579e8aeb) |

Divergence 3 is the dimensional fix. The formula in omega notation:
- wrong:   `mu = A * rho_lR^3 / (6pi) * (Sa/omega)^3 = A*Sa^3*nS^3*rho_SR^3/(6pi*n_l^3)` [Pa]
- correct: `mu = A / (6pi * rho_lR) * (Sa/omega)^3 = A*Sa^3*nS^3*rho_SR^3/(6pi*n_l^3*rho_lR)` [J/kg]

### What was changed in `RichardsMechanicsDSMMicroMacroBridge.mfront`

**@Integrator block — porosity helpers (new, replaces `n_s = 1.0 - phi0`):**
```mfront
const auto phi_M_from_nl = [&](double nl) {
  const auto nl_safe = std::max(minimum_n_l, nl);
  return std::clamp((phi0 - nl_safe) / (1.0 - nl_safe), 0.0, phi0);
};
const auto phi_m_from_nl = [&](double nl) {
  return nl * (1.0 - phi_M_from_nl(nl));
};
const auto active_nS_from_nl = [&](double nl) {
  return 1.0 - phi_M_from_nl(nl);   // = (1-phi0)/(1-nl)
};
const auto phi_m_prev = phi_m_from_nl(n_l_prev);   // was: n_l_prev
```

**omega_from_state lambda (uses evolving nS):**
```mfront
// was: n_l_arg * rho_lR_arg / (n_s * rho_SR)
const auto nS = active_nS_from_nl(n_l_arg);
return n_l_arg * rho_lR_arg / (nS * rho_SR);
```

**mu_micro_from_state lambda (J/kg fix):**
```mfront
// was: (hamaker * rho_lR_arg^3 / (6pi)) * (Sa/omega)^3   [Pa]
// now: (hamaker / (6pi * rho_lR_arg)) * (Sa/omega)^3     [J/kg]
return (hamaker_constant / (6.0 * pi * rho_lR_safe)) *
       std::pow(specific_surface / omega, 3.0);
```

**Post-solver phi_m/phi_M assignment:**
```mfront
// was: const auto phi_m_trial = n_l_trial;
//      const auto phi_M_trial = phi0 - phi_m_trial;
const auto phi_m_trial = phi_m_from_nl(n_l_trial);
const auto phi_M_trial = phi_M_from_nl(n_l_trial);
```

**@UpdateAuxiliaryStateVariables (same three fixes applied):**
```mfront
const auto phi_M_updated = std::clamp(
    (phi0 - std::max(minimum_n_l, n_l)) / (1.0 - std::max(minimum_n_l, n_l)),
    0.0, phi0);
const auto active_nS_updated = 1.0 - phi_M_updated;
phi_m = n_l * active_nS_updated;          // was: n_l
phi_M = phi_M_updated;                    // was: phi0 - n_l
const auto omega = std::max(n_l * rho_lR / (active_nS_updated * rho_SR), minimum_n_l);
mu_lR_value = (hamaker_constant / (6.0 * pi * rho_lR_safe)) *
              std::pow(specific_surface / omega, 3.0);   // [J/kg]
```

### Exchange coefficient unit convention after fix

| Domain | mu units | mass_exchange_coefficient units | Relation |
|---|---|---|---|
| Before fix (mfront) | Pa | s/m | `rho_hat [kg/m³/s] = alpha [s/m] * delta_p [Pa]` |
| After fix (mfront, native) | J/kg = m²/s² | kg·s/m⁵ | `rho_hat [kg/m³/s] = alpha [kg·s/m⁵] * delta_mu [J/kg]` |

Existing calibrated `MassExchangeCoefficient` values (e.g. `1e-13`) were in old
Pa-domain units. In J/kg domain the equivalent value is `old * rho_lR ≈ old * 1000`.
The parity PRJ files use `1e-10 kg·s/m⁵` (new domain).

---

## Parity PRJ files added (2026-05-25, commit de589d8fc0)

Location: `Tests/Data/RichardsMechanics/ANCHORS_MS33_StrictParity/`

| File | Runs on | Physics path |
|---|---|---|
| `ms33_dsm_parity_native.prj` | `dsm_native_hierarchical` binary | `<potential_exchange>` with `algebraic_split` + `current_porosity_split` |
| `ms33_dsm_parity_mfront.prj` | `dsm_mfront_hierarchical` binary | `RichardsMechanicsDSMMicroMacroBridge` |

Both files use identical physical parameters. See `README_DSM_PARITY.md` in that
directory for the parameter table, run instructions, and residual differences.

---

## Earlier partial port (2026-05-25, commit e279020d03) — NOW SUPERSEDED

The commit `e279020d03` applied a hierarchical porosity mapping to
`RichardsMechanicsDSMMicroMacroBridge_MCC.mfront` (the MCC variant) but
**not** to `RichardsMechanicsDSMMicroMacroBridge.mfront` (the main bridge).
It also did not fix the J/kg dimensional error or the evolving-nS issue.

This pass (commit de589d8fc0) fixes the main bridge comprehensively. The MCC
variant (`_MCC.mfront`) should receive the same three fixes in a follow-up.

---

## Replay instructions from clean `master`

```bash
git checkout master
git checkout -b dsm_mfront_hierarchical
# Replay all commits:
git format-patch dsm_mfront..dsm_mfront_hierarchical -o /tmp/dsm_mfront_hier_patches
git am /tmp/dsm_mfront_hier_patches/*.patch
```

Or cherry-pick specific physics commits:
- `e279020d03` — MCC variant partial hierarchical port
- `de589d8fc0` — main bridge: all three divergences fixed + parity PRJ files

---

## Additional physics fixes (2026-05-25, commit TBD)

Four more structural divergences discovered and fixed by tracing
`RichardsMechanicsFEM-impl.h` and running the parity simulations:

| # | Was (WRONG) | Is (CORRECT) |
|---|---|---|
| 4 | `mu_micro = A/(6π·ρ_lR)·(Sa/ω)³` [wrong ρ_lR power via ω] | `mu_micro = −A·Sa³·nS³·ρ_SR³/(6π·n_l³·ρ_lR)` [J/kg, negative] |
| 5 | Missing factor ρ_LR/μ in exchange: `rho_l_hat = α·(μ_LR − μ_lR)` | `alpha_eff = α·ρ_LR_ref/μ_fluid; rho_l_hat = alpha_eff·(μ_LR − μ_lR)` |
| 6 | Aggregate-scale mass balance: `n_l·ρ_lR` | REV-scale: `φ_m·ρ_lR` (matches native `rho_l_prev = phi_m_prev * rho_lR_prev`) |
| 7 | Newton line search allowed n_l > φ₀ (through clamped phi_M formula) | Clamp `n_l_candidate` to `[minimum_n_l, φ₀]` in line search; saturation fallback when equilibrium demands n_l > φ₀ |

Additional fix: `@UpdateAuxiliaryStateVariables` was recomputing `rho_l_hat_value`
using `p_LR` (beginning-of-step value in OGS-MFront coupling) instead of the
correctly computed `rho_l_hat_out` already set in `@Integrator`. Removed the
redundant recomputation.

`FluidViscosity` material property added to bridge (was missing) and to
`ms33_dsm_parity_mfront.prj`.

`ms33_dsm_parity_native.prj` fixed: `local_nonlinear_solve_mode` changed from
`full_potential` to `scalar_micro_macro_mass_storage_mode` (only valid mode);
removed nonexistent `<mode>` child element; added micro liquid EOS parameters.

---

## Validation status (2026-05-25, second pass)

| Check | Status | Notes |
|---|---|---|
| Build | ✓ DONE | `ninja build_mfront` rebuilds `libOgsMFrontBehaviour.dylib` |
| ctest | NOT RUN | Parity PRJ files not yet registered in Tests.cmake |
| Parity run (native vs mfront) | ✓ DONE | Both complete 120 steps without error |
| n_l parity | ✓ machine-epsilon at ts≥60, 1.2e-6 at ts=10 | n_l is the primary physics state variable |
| phi_m parity | ✓ machine-epsilon at ts≥60, 1.7e-6 at ts=10 | |
| rho_l_hat parity | ✓ machine-epsilon at ts=120, 2.8e-6 at ts=10 | After removing aux-block recomputation |
| sigma parity | ✓ 0.1% at ts≥60; 20% at ts=10 | Early-time difference from micro-EOS approximation (rho_l0=1, density_a=50) |
| Mathematical correctness | ✓ VERIFIED | vdW formula, alpha_eff, REV-scale mass balance all match native |
| PRJ parameter consistency | ✓ VERIFIED | Both parity PRJs use same numerical values |

---

---

## Physics fixes ported to MCC variant (2026-05-25, commit TBD)

Applied the same four physics fixes from the main bridge to
`RichardsMechanicsDSMMicroMacroBridge_MCC.mfront`:

| Fix | Location in _MCC.mfront | Change |
|---|---|---|
| 1. vdW nS-form + `/rho_lR` | `mu_micro_from_state` lambda | Replaced omega-form (`A·rho_lR³/(6π)·(Sa/ω)³`) with nS-form: `sign·A·Sa³·nS³·rho_SR³/(6π·n_l³·rho_lR)`. Keeps existing `micro_potential_sign` for convention switching. |
| 2. alpha_eff reference density | `alpha_M_effective_from_density` | Changed numerator from dynamic `rho_macro_safe` to fixed `rho_LR_ref`. `macro_viscosity` already served as FluidViscosity — no new material property needed. |
| 3. REV-scale mass balance | `evaluate`, `evaluate_reduced` | **Already correct** — MCC file was already using `phi_m·rho_lR` (not `n_l·rho_lR`). No change needed. |
| 4. Newton ceiling clamp + saturation fallback | Coupled solver line search; `solve_microstate_bracketed`; post-loop fallback | Changed clamp upper bound from `n_l_upper_bound` to `phi0` in line search and bracketed solver. Added saturation fallback (clamp to phi0 when `residual_at_phi0 ≤ 0`) in `!accepted` branch and final fallback. |
| 5. No rho_l_hat recomputation in @UpdateAuxiliaryStateVariables | `@UpdateAuxiliaryStateVariables` | **Already correct** — MCC file copies `rho_l_hat_trial_value` (set by `solve_microstate` in `@InitializeLocalVariables` using end-of-step pressure). No recomputation from beginning-of-step `p_LR`. |

Build: `ninja build_mfront` — passes clean (all 55 targets, `libBehaviour.dylib` linked).

---

## DSM verification ramps (2026-05-25) — single-element replacement for Jupyter GP notebook

Location: `Tests/Data/RichardsMechanics/ANCHORS_MS33_DSM_VerificationRamps/`

Replaces the Gauss-point Jupyter notebook
`/Users/vinaykumar/tex/cc2024/examples/DSM/Potentials_penalty.ipynb`
(saturation + desaturation ramp on a single Gauss point) by running the
equivalent ramp as an OGS single-element test, at the three Model-I
calibrated dry densities (1400, 1600, 1800 kg/m³).

- **Template PRJ**: `ANCHORS_MS33_StrictParity/ms33_dsm_parity_native.prj`
  (DSM hierarchical path active; micro_* secondary variables wired up).
- **Per-density override**: `phi0` and `IntrinsicPermeability0` substituted
  from `ANCHORS_MS33_ModelI/ms33_model_i_dd{1400,1600,1800}.prj`.
- **Ramp**: triangle in pressure BC; suction goes
  `-100 MPa → +0.1 MPa → -100 MPa` over 120 days (dt = 1 day).
  Shape matches the notebook (saturation half then desaturation half);
  pressure range covers the suction regime over which the Model-I
  calibration is active.
- **Mesh**: reuses the single-element axisymmetric quad
  `../square_1x1_quad_1e0.vtu` from the sibling MS33 directories.
- **Binary**: native DSM hierarchical
  `/Users/vinaykumar/git/build/native-release-omp-sharedcache/bin/ogs`
  (worktree `dsm_native_hierarchical`).

**Run script**:
```bash
cd Tests/Data/RichardsMechanics/ANCHORS_MS33_DSM_VerificationRamps
python3 run_verification_ramps.py            # generate PRJs, run all three, plot
python3 run_verification_ramps.py --only-gen # regenerate PRJs only
python3 run_verification_ramps.py --skip-run # skip OGS, just post-process
```

The script does three things in sequence:
1. Generates `dsm_ramp_dd{1400,1600,1800}.prj` from the native parity
   template with per-density phi0/perm substitutions and the verification
   ramp BCs.
2. Runs each PRJ with the native binary; logs go to `dsm_ramp_dd{dd}.log`.
3. Extracts cell/point-averaged scalars to `dsm_ramp_dd{dd}.csv` and writes
   four PNGs per density (water-content vs time, water-content vs suction,
   macro saturation hysteresis, mean stress vs time).

**Verified outcome (2026-05-25)**: All three sims complete 120 steps without
rejection. dd1400 example: at peak saturation (t=55 d, p=-8.2 MPa) the micro
fully saturates (`n_l = phi0 = 0.4964`), `phi_M → 0` (rigid hierarchical
split closing the macro channel), and the swelling-stress trace reaches its
peak ~−2.15 MPa; on desaturation it returns close to the initial state with
visible hysteresis, as expected from the path-dependent DSM update.

The four PNGs per density live in the same directory and are intended as
drop-in replacements for the Jupyter-notebook figures currently in the
manuscript. They are NOT yet swapped into `paper_DSM.tex`.

---

## Open items

- [x] Apply same fixes to `RichardsMechanicsDSMMicroMacroBridge_MCC.mfront` — DONE 2026-05-25
- [x] Add parity runner script — DONE 2026-05-25; see `scripts/run_dsm_parity.py`
- [x] Replace Gauss-point Jupyter verification ramps with single-element OGS
      sims at varying dd — DONE 2026-05-25; see
      `Tests/Data/RichardsMechanics/ANCHORS_MS33_DSM_VerificationRamps/run_verification_ramps.py`
- [ ] Swap notebook figures in `paper_DSM.tex` for the new PNGs (user review pending)
- [ ] Register parity PRJ files in `Tests/Data/RichardsMechanics/Tests.cmake`
- [ ] Recalibrate `MassExchangeCoefficient` and `HamakerConstant` in existing
      BEACON PRJ files (those using old Pa-domain alpha values) to the J/kg domain
- [ ] Investigate early-time sigma discrepancy (20% at ts=10) — may be micro-EOS
      or Biot coupling initialization difference

---

## Parity runner script (2026-05-25, commit d2941d2aa7)

`scripts/run_dsm_parity.py` — canonical way to run native-vs-mfront comparisons.

```bash
python3 scripts/run_dsm_parity.py              # run all suites
python3 scripts/run_dsm_parity.py --suite dsm_ms33
python3 scripts/run_dsm_parity.py --no-run     # re-report without re-running
```

**To add a new model pair:** append one dict to `PARITY_SUITES` at the bottom of the
script. Fields: `name`, `description`, `native_prj`, `mfront_prj`, `field_map`, and
optionally `native_bin`, `mfront_bin`, `notes`. No other code changes required.

**Updated docs** (same commit):
- `Tests/Data/RichardsMechanics/SCRIPT_INVENTORY.md` — new entry in workflows table
- `Tests/Data/RichardsMechanics/ANCHORS_MS33_StrictParity/README_DSM_PARITY.md` — "How to run" section now points to script
- `~/.claude/projects/.../memory/reference_dsm_parity_script.md` — agent memory

---

## OPEN PHYSICS QUESTIONS (2026-05-25, mandatory investigation before publication)

### Q1 — Exchange coefficient density scaling

Current implementation:
```
alpha_eff = MassExchangeCoefficient * rho_LR_ref / fluid_viscosity
rho_hat   = alpha_eff * (mu_LR - mu_lR)                          [kg/m³_REV/s]
```
The effective exchange coefficient scales with the **macro** fluid density
`rho_LR_ref`. This matches `alpha_M_effective = alpha_bar * rho_LR / mu` in
`RichardsMechanicsFEM-impl.h`.

**MUST RESOLVE**: Is `rho_LR` the correct density here, or should the
exchange flux involve the respective micro/macro densities, or their average?
The question is non-trivial when `rho_lR ≠ rho_LR` (compressible micro fluid).
The manuscript must state the derivation and justify the choice.

Tracked in: `tex/dsm-bgr-paper/draft/paper_DSM.tex` `\VKtodo` near exchange equation;
            `tex/dsm-bgr-paper/AGENTS.md` section "2026-05-25 open physics questions".

### Q2 — Micro-saturation ceiling and large-alpha timestep behaviour

With `alpha_eff = 1e-10 * 1000 / 1e-3 = 1e-4`, exchange is so fast that
within one timestep the equilibrium `n_l` can exceed the current total
porosity `phi`. The native C++ path clamps `n_l <= phi` (current, evolving).

The MFront bridge was clamping at the **initial** `phi_0` and throwing an
exception when no root existed below that ceiling. Fix applied (commit
`ca360695b3`): ceiling clamp in Newton line search + saturation fallback that
returns `n_l = phi_0` when residual at ceiling is negative.

**MUST RESOLVE**:
1. Should the ceiling be `phi` (current, evolving with swelling) or `phi_0`
   (fixed, small-strain approximation)? Using `phi_0` is a known approximation.
2. Does the saturation-clamped state (n_l held at ceiling, exchange not fully
   equilibrated in one step) introduce a spurious artificial source term in the
   macro balance over multiple steps?
3. Is there an implication for time-step sensitivity — should a maximum
   timestep be imposed when `alpha_eff` is large?

Tracked in: `tex/dsm-bgr-paper/AGENTS.md` section "2026-05-25 open physics questions".

---

## DECISION (2026-05-26): mfront must use current evolving `phi`, not `phi_0`

### Decision

The MFront bridge SHALL match the native implementation in using the
**current evolving total porosity** `phi(t)` everywhere the hierarchical
algebraic split and the saturation ceiling appear. `phi_0` (initial
porosity) is no longer the constitutive reference; it is only the initial
condition of `phi(t)`. This aligns the paper's stated convention with both
implementations.

### Justification (mass-conservation consistency)

The paper now contains a derivation (paper_DSM.tex §2.2, eq.
`eq:micro_balance_cancellation` and `eq:micro_balance_after_cancellation`)
showing that with `phi_m = n_l * (1 - phi(t)) / (1 - n_l)` from the
hierarchical split using the EVOLVING phi(t), the explicit REV
deformation term `phi_m * rho_lR * div(u_dot)` in the micro mass balance
and the implicit `dot(phi)` contribution that propagates through the
algebraic relation **cancel exactly**. The cancellation is

    dot(phi_m * rho_lR) + phi_m * rho_lR * div(u_dot)
      = [dot(n_l) * rho_lR * (1 - phi) / (1 - n_l)^2]
        - phi_m * rho_lR * div(u_dot)         <-- from dot(phi)=(1-phi)*div(u_dot)
        + phi_m * dot(rho_lR)
        + phi_m * rho_lR * div(u_dot)         <-- explicit REV term
      = dot(n_l) * rho_lR * (1 - phi) / (1 - n_l)^2 + phi_m * dot(rho_lR)
      = rho_l_hat.

The aggregate-referenced micro liquid content n_l is, by construction of
the hierarchical split, blind to REV-scale mechanical dilation. The two
appearances of `div(u_dot)` in the unexpanded balance are NOT the same
effect counted twice; they precisely cancel. Mass conservation is exact
under any mechanical loading. Using `phi_0` in the algebraic split is a
small-strain approximation of this exact form; using `phi(t)` is the
exact form.

### Implementation changes required in mfront

The MFront bridge file
`MaterialLib/SolidModels/MFront/RichardsMechanicsDSMMicroMacroBridge.mfront`
currently uses `phi0` in several places. Replace each with current `phi`
read from the OGS side. Specifically:

1. **`phi_M_from_nl` lambda in `@Integrator`**: change
       (phi0 - nl_safe) / (1.0 - nl_safe)
   to use current phi (passed as an external state variable or computed
   from PorosityFromMassBalance result):
       (phi - nl_safe) / (1.0 - nl_safe)
   with the corresponding clamp upper bound updated from `phi0` to `phi`.

2. **Newton line-search ceiling** in `solve_microstate`: change
       std::clamp(n_l_candidate, minimum_n_l, phi0)
   to
       std::clamp(n_l_candidate, minimum_n_l, phi)
   so the ceiling moves with mechanical dilation.

3. **Saturation fallback** (both occurrences — line-search fallback and
   post-Newton fallback): the ceiling probe currently evaluates at
   `phi0`; change to `phi`. The fallback sets `n_l_out = phi` when the
   residual at ceiling is non-positive.

4. **`@UpdateAuxiliaryStateVariables`**: any `phi_M_updated` computed
   with `phi0` must use current `phi` from the same external-state input.

5. **External state variable wiring**: add a `@ExternalStateVariable`
   for `TotalPorosity` (named `phi_total` in mfront, mapped to OGS's
   current `phi` from `PorosityFromMassBalance`). The OGS side already
   evolves `phi(t)` from the solid balance and PorosityFromMassBalance.

### What this does NOT change

- The vdW and augmentation closures (`mu_micro_from_state`,
  `mu_aug_from_state`) — these depend on `n_l`, `n_S`, `rho_SR`, `s^a`
  via the existing formulas; they do not involve `phi_0` directly.
- The exchange-coefficient form `alpha_eff = alpha_bar * rho_LR / mu`.
- The macro/micro liquid balance forms in the OGS C++ side.
- The native binary is already correct in using current `phi(t)`; only
  the MFront bridge needs the change.

### Validation steps after the change

1. Re-run the strict-parity benchmark
   `Tests/Data/RichardsMechanics/ANCHORS_MS33_StrictParity/`
   (`ms33_dsm_parity_native.prj` vs `ms33_dsm_parity_mfront.prj`).
   In the parity test case (zero deformation, all-Dirichlet pressure
   BCs), the change should be a no-op since `phi(t) = phi_0` pointwise.
   Confirm machine-epsilon parity is preserved.
2. Re-run a coupled HM benchmark with non-zero `div(u_dot)` (any of
   the BEACON or MS33 coupled cases) and confirm that native and mfront
   results converge as the timestep shrinks (with the cancellation,
   they should agree exactly; previously they diverged in regions with
   `div(u_dot) != 0`).
3. Re-run the DSM verification ramps in
   `ANCHORS_MS33_DSM_VerificationRamps/`. These are confined
   single-element runs so the change is a no-op; outputs should be
   bit-identical to the current committed CSVs.

### Cross-references

- Paper derivation: `paper_DSM.tex` §2.2, eq.
  `eq:micro_balance_cancellation` and `eq:micro_balance_after_cancellation`
  (in red).
- Open-questions tracking: the Q2 entry above (micro-saturation ceiling)
  is resolved by this decision; `phi(t)` is the chosen convention.
- The kinematic-partition limitation (the rigid split forcing $\phi_M$
  and $\phi_m$ to share $\dot\phi$ in a fixed ratio) is a separate
  issue, addressed in paper_DSM.tex `\subsection{Discussion of results}`
  and `\subsection{Outlook}`. This decision does not change that
  limitation; it only fixes the small-strain approximation in the
  saturation ceiling.
