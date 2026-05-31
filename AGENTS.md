# AGENTS worklog — dsm_mfront_hierarchical worktree

Worklog for the pure-MFront DSM bridge worktree. Accretes over time;
old entries stay (marked historical), never deleted (CLAUDE.md §6/§11).

---

## 2026-05-31 — Full MS33 suite attempt across mfront branches (Vinay request)

Request: get Model I/III/IV/V → 200 d, VII → 240 d running in
(A) this worktree = pure MFront DSM bridge, LE + MCC; and
(B) `dsm_native_hierarchical_wt` MCC_NATIVE = native DSM + MFront
ModifiedCamClay. Task #38.

### Current matrix (verified 2026-05-31)

| model | native+mfront-MCC (B) | mfront LE (A) | mfront MCC (A) |
| :---- | :-------------------- | :------------ | :------------- |
| I ×3  | 200 d ✅              | 200 d ✅      | 200 d ✅       |
| III   | 200 d ✅              | ❌ blocker ①  | ❌ blocker ①   |
| IV    | 200 d ✅              | ❌ blocker ①  | ❌ blocker ①   |
| V     | 200 d ✅              | no prj in wt  | no prj in wt   |
| VII   | ❌ ~6.5 d (blocker ②) | ❌ blocker ①  | ❌ ①, then ②   |

Model I (and V) work in (A) only because they pin pressure on BOTH
ends → zero free interior pressure DOFs. III/IV/VII have free interior
pressure nodes → blocker ① triggers.

### Blocker ① — mfront DSM bridge: missing pressure-Jacobian (FIXABLE, design ready)

Symptom: III/IV/VII die at time step #1, Newton iter #2, "Failed during
Eigen linear solver initialization" (SparseLU). Root cause is iter #1:
the liquid pressure update blows up (`|dx|≈7.87e15`, component 0) while
displacement components are sane.

Root cause (confirmed by code read, native is the known-good control):
the bridge assembles the DSM liquid mass-exchange source into the
pressure **residual but not the pressure Jacobian**.
- Defect: `ProcessLib/RichardsMechanics/RichardsMechanicsFEM-impl.h`
  lines 1606–1617 — `local_rhs += N_p^T·liquid_mass_exchange_source·w`
  with an inline comment explicitly deferring the Jacobian term. No
  `local_Jac` contribution.
- Known-good native: `dsm_native_hierarchical_wt/.../RichardsMechanicsFEM-impl.h`
  3559–3569 adds BOTH the residual AND
  `local_Jac[p,p] -= N_p^T·drho_L_hat_dpL_direct·N_p·w`. Native
  derivative in `computePotentialExchangeUpdate` (native 130–228).

At a free pressure node the only pressure-row Jacobian terms are storage
(tiny at high suction) + Laplacian (tiny: k=1e-18, k_rel floored 1e-2).
The DSM source adds a finite residual forcing with ZERO diagonal
sensitivity → Newton dx ≈ residual / (≈0 diagonal) → 1e15 → iter-2
tangent unfactorable. The converged SOLUTION is defined by the
residual, which is already correct; the missing tangent only breaks
the Newton path → a consistent-tangent fix changes convergence, NOT
physics.

DECISION (Vinay, 2026-05-31): implement the **analytic** tangent
(not FD). VII-MCC: investigate apex first (blocker ②), no physics change.

#### Exact fix (analytic dominant term — NO new literal)

The bridge `.mfront` (`RichardsMechanicsDSMMicroMacroBridge.mfront`
269–272) computes
`rho_l_hat = alpha_eff·(mu_LR − mu_lR)`, `alpha_eff = alpha·rho_LR/mu`,
`mu_LR = mu_macro_from_pressure(p, rho)`.
Dominant sensitivity (the macro chemical-potential channel):
`∂rho_l_hat/∂p_L = alpha_eff·∂mu_LR/∂p = alpha_eff·(1/rho_LR) = alpha/mu`.
Source = `−rho_l_hat`; DOF is `p_cap = −p_L`; following this file's own
`saturation_micro` assembly pattern (rhs 1626–1628 ↔ Jac 1630–1633),
the contribution is a stabilizing positive diagonal:

```
local_Jac[pressure,pressure] += N_p.transpose() * (alpha/mu) * N_p * w;  // 1/(Pa·s) · ... → pressure-block units, matches saturation_micro term
```

`alpha` (bridge mass-exchange coefficient) is an MFront material
property — NOT available C++-side; only `mu` (LiquidViscosityData, line
1518) and `rho_LR` (1444) are. So route = **MFront exports the
derivative**:

1. `.mfront`: export a new aux/internal scalar
   `drho_l_hat_dpLR = alpha_eff * dmu_macro_from_pressure_dp`
   (`mu_macro_from_pressure` is Young–Laplace, analytic derivative
   available in-file). No new physical constant — reuses existing
   `alpha_eff`, `mu_LR` machinery.
2. `MFrontRichardsMechanics.h` (~85–98): extract via
   `getInternalScalarByName(*state, "drho_l_hat_dpLR")`, add field to
   `PressureCoupledResponse`.
3. `ConstitutiveRelations/PressureCoupledSolidData.h` (struct line 35 +
   `makePressureCoupledSolidData` 43–52): add
   `d_liquid_mass_exchange_source_dp_cap`, set with the **p_cap = −p_L
   sign flip** exactly like `dS_L_dp_cap` (line 47) — i.e. negate the
   MFront `d/dp_L` to get `d/dp_cap`. Carry the `−rho_l_hat` source
   sign through.
4. `RichardsMechanicsFEM-impl.h` 1606–1617: add the
   `local_Jac[p,p] += N_p^T·(d_source_dp_cap)·N_p·w` term beside the
   existing residual line. §4.1: `+=` (accumulator). §4.2: annotate units.

Sign safety net: derived sign = **positive diagonal** (`+alpha/mu`,
matching `saturation_micro`). If III still fails to converge after the
fix, FLIP the sign and rebuild — a wrong-sign tangent only worsens
convergence; it does NOT change the converged solution (residual
untouched). A clean III run to 200 d self-verifies the sign.

#### Build + run (next session — long op, do as one uninterrupted block)

```
# regen MFront behaviour + rebuild ogs (verify MFront target name):
cmake --build /Users/vinaykumar/git/build/dsm_mfront_hier_wt-release -t ogs -j
# binary: /Users/vinaykumar/git/build/dsm_mfront_hier_wt-release/bin/ogs
# run (from each prj dir, OMP_NUM_THREADS=1, -o <dir>/out_works):
#   ANCHORS_MS33_ModelIII/ms33_modelIII_gap2mm.prj           (LE)
#   ANCHORS_MS33_ModelIV/ms33_modelIV_pellets.prj            (LE)
#   ANCHORS_MS33_ModelVII/ms33_modelVII_freeswelling.prj     (LE, 240 d)
#   ANCHORS_MS33_MCC/ModelIII|ModelIV|ModelVII/*.prj         (MCC)
# Regression guard (§12.3): re-run ModelI (single elem) and confirm
# UNCHANGED (the new term is inactive when pressure is fully pinned, and
# only affects convergence elsewhere). MCC VII will still hit blocker ②.
```

EXECUTED + VERIFIED 2026-05-31 (analytic-direct term): implemented as
the direct macro channel `d(rho_l_hat)/d(p_L) = alpha_eff/rho_macro_safe
= mass_exchange_coefficient/viscosity`, exported as MFront aux var
`drho_l_hat_dpLR` from both .mfront files (@UpdateAuxiliaryStateVariables),
threaded through PressureCoupledResponse → PressureCoupledSolidData
(p_cap=-p_L flip) → RM assembler (`local_Jac[p,p] += N_p^T·(alpha/mu)·N_p·w`,
positive diagonal). Built clean (ogs relinked, both -generic.cxx carry the
new var).

RESULT — PARTIAL, INSUFFICIENT ALONE:
- Model I dd1600 LE: 308 steps → 200 d, UNCHANGED (regression-safe; pinned
  pressure → new diagonal inert). Tree is in a clean functional state.
- Model III LE: blow-up reduced 5 orders (pressure |dx| 7.87e15 → ~4e10)
  but STILL DIVERGES (|dx| stuck ~1e10, Newton oscillates, time stepper
  bottoms out → exit 1).
- ROOT CAUSE of insufficiency (architectural): the bridge re-solves the
  FULL microstate (n_l, rho_lR, mu_lR) implicitly inside MFront at every
  macro pressure, so the true d(rho_l_hat)/d(p) is dominated by the
  IMPLICIT microstate (vdW) response — ~5 orders larger than the direct
  macro channel. Native converges with direct-only ONLY because native
  LAGS its microstate (separate state var), making the direct channel
  consistent there. The bridge's tighter implicit coupling needs the FULL
  tangent.
- The direct term is correct (right sign, regression-safe) but incomplete;
  it is LEFT IN PLACE (harmless, improves conditioning).

NEXT (Vinay's call, 2026-05-31): direct-analytic was Vinay's earlier
choice; it is now shown insufficient for the bridge architecture. Options
re-offered: (A) FD tangent — perturb macro p, re-solve microstate, FD
rho_l_hat; robust, full sensitivity, matches native FD mode; needs the
scale-aware perturbation h=1e-8·max(1,|p|) (§1.2 scale-derived; native's
own approved default, §1.1-item-3). (B) full analytic microstate-chain
tangent (implicit-function-theorem on the 2×2 micro Jacobian → dn_l/dp,
drho_lR/dp → dmu_lR/dp); no literal but substantial algebra/risk.
(C) stop, document III/IV/VII mfront as still-blocked.

### Blocker ② — VII MCC cam-clay tension apex (PHYSICS, Vinay's call)

native+mfront-MCC VII dies ~6.3–6.5 d (`t≈5.0e5–5.6e5 s`), MFront
ModifiedCamClay status −1. Furthest outputs:
`dsm_native_hierarchical_wt/.../ANCHORS_MS33_MCC_NATIVE/ModelVII/_tens_out`
(t≈561726), `_tenscut_out`/`_tensdiag_out` (≈503254). This is the
documented MCC dry-side limit (VII is LE-only;
[[project_dsm_mcc_bishop_cutoff]]).
DECISION (Vinay): investigate the apex (stress path, which IP, distance
to cap, S_L/suction at failure) before any physics change. Read-only
agent was launched 2026-05-31 but hit the account session limit before
reporting — **re-run next session** (task #40). No physics literal to
be chosen without Vinay (§9/§1.1).

### LAG ATTEMPT RESULT (2026-05-31, Vinay-approved) — also insufficient; root cause is the SOURCE MAGNITUDE, not the tangent

Implemented the native-matching lag in the base .mfront: macro source
`rho_l_hat = alpha_eff*(mu_LR(p_current) - mu_lR_lagged)` with mu_lR frozen at
the begin-of-step micro potential (`mu_micro_from_state(n_l_prev,rho_lR_prev)`),
microstate solver UNCHANGED, consistent direct tangent +alpha_eff/rho_macro_safe.
Clear inline comment block added at the @Integrator store. Built clean.

Model III LE STILL diverges, but the signature changed: pressure |dx|=|x|=
1.6173e17, ratio exactly 1.0, IDENTICAL on every time-step retry (no more
feedback oscillation — the lag DID remove the feedback). The constant huge
value = the macro pressure trying to equilibrate to the frozen micro potential
mu_lR_lagged, which at the dry IC is absurdly negative (vdW ~ 1/n_l^3). So the
true blocker is the SOURCE MAGNITUDE: the bridge dumps the dry clay's full water
demand onto the macro pressure equation as a sink, and with NO separate
micro-pressure field to absorb it (native HAS one), free interior pressure DOFs
(III/IV/VII) can only balance it with an unphysical ~1e17 Pa excursion. Model I
survives only because its pinned pressure is supplied at the boundary.

Approaches tried + result (all on Model III LE, step-1 pressure |dx|):
  none 7.87e15 | analytic-direct 4e10 | FD-full 1.7e11 | lag 1.6e17. All diverge.
CONCLUSION: not fixable by tangent or lag. Pure-mfront multi-element needs a
fundamental architecture change (separate micro-pressure carrier field as in
native, OR a macro-pressure cavitation cap, OR an equilibrated/less-extreme
micro IC) — physics/architecture decisions, large effort, low marginal value
since native already covers it. Recommended: revert bridge to pristine, document
the limitation, rely on native LE (all models) + native+mfront-MCC (I/III/IV/V).

DECISION (Vinay, 2026-05-31): REVERT & rely on native. All 6 bridge source
files (the two .mfront, MFrontRichardsMechanics.h, MechanicsBase.h,
PressureCoupledSolidData.h, RichardsMechanicsFEM-impl.h) reverted to pristine
via `git checkout --` and ogs rebuilt to match. PRJs and this AGENTS.md left
intact. The pure-MFront DSM bridge is documented as supporting single-element
(all-pressure-pinned) problems only — Model I LE+MCC — and NOT multi-element
axisym (III/IV/VII), which require either a separate micro-pressure carrier
field (native's architecture) or a macro-pressure cavitation cap. Deliverable
stands on: native LE (I/III/IV/V→200d, VII→240d) + native+mfront-MCC
(I/III/IV/V→200d). VII-MCC apex remains the separate open physics item (#40,
"investigate first" — the investigation agent hit the account limit; re-run
pending).

### Notes / provenance flags carried from the 3-branch run (#34, done)
- mfront Model I dd1400 & dd1800 prjs carry `VdwAugmentationPrefactor=0`
  (uncalibrated; native K=7654.9 / 118582.6 J/kg not ported). Their
  reported swelling pressures are NOT calibrated values. §12 gap.
- native MCC dd1800 = −6.09 MPa < dd1600 = −14.66 MPa: non-monotonic
  ordering flagged predicted-not-verified (possible stale `_farm_out`
  config). Vinay to verify.
- §12.3 gating LE ctest: reltols 1e-6 verified bit-for-bit identical vs
  pre-reltol references on all 5 gating models (native).

## PENDING / TODO (noted at the 2026-05-31 "commit all" checkpoint)

Committed under a Vinay-approved §0.2 exemption (the modified ANCHORS_MS33
PRJs carry material-param changes — VdwAugmentationPrefactor K=71900 J/kg,
VdwAugmentationDecayLength λ=7.5e-7 m, min_relative_permeability 1e-12→1e-2 —
and several lack the §12.2 provenance header). OUTSTANDING:

1. **§12.2 provenance backfill (BLOCKING for §12.3 compliance).** Add the
   mandatory provenance header to every modified DSM PRJ that lacks it
   (Model III, MCC III/IV/VII at minimum): calibration anchor (Dixon/Villar)
   + dataset row + target sigma + fitted K, and per-group material-param
   sources. The committed K=71900 / λ=7.5e-7 / relperm-floor values currently
   have NO cited source on record — Vinay to supply.
2. **mfront Model I dd1400 & dd1800: VdwAugmentationPrefactor=0** (uncalibrated;
   native K=7654.9 / 118582.6 J/kg not ported). Swelling pressures are not the
   calibrated values. Decide: port native K, or document as intentionally bare.
3. **native MCC dd1800 ordering anomaly** (−6.09 < dd1600 −14.66 MPa): verify
   (possible stale `_farm_out` config).
4. **VII-MCC tension apex** (task #40): re-run the read-only apex investigation
   (prior agent hit the account session limit), then decide the physics.
5. **Housekeeping:** run-output dirs (`out_works/`, `_spec_*/`, `_t/`, `_probe/`,
   `*_run.log`) are NOT gitignored and were deliberately EXCLUDED from the
   commit. Add a `.gitignore` for these patterns so they stop showing as
   uncommitted and can't be accidentally committed.
