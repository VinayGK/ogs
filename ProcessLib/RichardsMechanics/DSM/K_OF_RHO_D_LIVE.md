# Live K(rho_d) — augmentation prefactor at the EVOLVING dry density

Status: first-cut trial implementation (2026-06-11), per Vinay's order
"K(rho_d) try it" (2026-06-10). Branch `dsm_native_h_of_eps`.

READ FIRST: the LIVE path interpolates LOG-LINEARLY since commit 1bb414ac05
(Vinay's interpolation-scheme decision, 2026-08-26). The normative statement
of the shipped value and tangent is the closing section "Interpolation
scheme: LOG-LINEAR" at the end of this file; passages above it that describe
K-linear interpolation between knots or a per-segment-CONSTANT tangent are
annotated `[HISTORICAL ...]` in place and describe the pre-1bb414ac05 scheme.

READ FIRST (2): two further changes sit UNCOMMITTED in this worktree as of
2026-08-31 (code-review rounds 1 and 2 on top of HEAD 7ec39ecf4c; worklog
`AGENTS.md` entries 32-38). Node preservation AT AN INTERIOR KNOT is now
enforced by an explicit branch rather than by the algebra, and every
`<prefactors>` entry must be STRICTLY POSITIVE at parse time. Both are
specified below — in "PRJ interface" and in the closing section — and
neither has been exercised by an MS33 run or by ctest.

## What it is

The existing parse-time K(rho_d) table feature (memory:
`project_dsm_k_of_dry_density`, branch `dsm_native_pdisj_maxwell_kofdd`;
ported to this branch) resolves the augmentation prefactor
K = `potential_augmentation_prefactor` ONCE at parse time from a
piecewise-linear table `potential_augmentation_prefactor_vs_dry_density`
evaluated at the PRJ-supplied initial/target `<dry_density>`. K is then a
per-material constant in time.

This LIVE variant keeps the table un-frozen and re-evaluates

    K = K_table(rho_d),    rho_d = rho_SR * (1 - phi)   [kg/m^3]

at run time, where `rho_SR` = `micro_solid_density_reference` and `phi` is
the CURRENT total porosity at the evaluation site. As the column swells
(phi up, rho_d down) the prefactor tracks the calibrated K(rho_d) curve.

## PRJ interface

```xml
<potential_exchange>
  ...
  <potential_augmentation_prefactor_vs_dry_density>
    <dry_densities>...</dry_densities>
    <prefactors>...</prefactors>
  </potential_augmentation_prefactor_vs_dry_density>
  <dry_density>1400</dry_density>  <!-- optional in live mode: fallback K -->
  <potential_augmentation_prefactor_live_dry_density>true</potential_augmentation_prefactor_live_dry_density>
</potential_exchange>
```

- `potential_augmentation_prefactor_live_dry_density` default `false` →
  parse-time freeze, bit-for-bit the existing behavior (verified: dd1400
  off-mode regression below).
- `true` REQUIRES the table (OGS_FATAL otherwise) and is mutually exclusive
  with the scalar `potential_augmentation_prefactor` (unchanged rule).
- In live mode `<dry_density>` is optional; when given, K(dry_density) is
  stored as the FALLBACK scalar used at evaluation sites that have no
  porosity in scope (see below).
- **Every `<prefactors>` entry must be STRICTLY POSITIVE (> 0)** — an
  authoring rule new on 2026-08-31, uncommitted in the code-review batch
  (`AGENTS.md` entry 33). `parsePotentialExchangeParameters` in
  `CreateRichardsMechanicsProcess.cpp` walks the parsed `<prefactors>` list
  before the table is constructed and `OGS_FATAL`s on the first entry that
  is not `> 0`, naming that entry's index in the PRJ list and its value.
  Reason: the LIVE interpolant evaluates `ln(K_r/K_l)` (closing section
  below), which a zero or negative knot poisons — in one of the cases
  silently, through the Jacobian only. The check is deliberately NOT
  restricted to live mode, so it WIDENS parser rejection: a frozen-K deck
  carrying a zero knot used to parse and no longer does. No shipped deck is
  affected — all seven `<prefactors>` tables under
  `Tests/Data/RichardsMechanics/` are strictly positive (re-grepped
  2026-08-31: the MS33 III/IV/VII gapswitch / pellets / freeswelling decks,
  plus the four `*_kofdd` / `*_livek` decks).
- `<dry_densities>` needs no monotonicity rule and none was added: the
  MathLib base constructor sorts the support points (carrying their values
  with them) and `OGS_FATAL`s on duplicate support points
  (`MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.cpp`,
  re-read 2026-08-31).

## Implementation

One inline helper, `effectiveAugmentationPrefactor(params, phi)` in
`PotentialExchangeParameters.h`:

- live mode && table && `std::isfinite(phi)` →
  `table->getValueLogLinear(rho_SR * (1 - phi))` [LOG-LINEAR since commit
  1bb414ac05, Vinay's interpolation-scheme decision 2026-08-26; the first
  cut called the K-linear `table->getValue(...)`];
- otherwise → the parse-time scalar `potential_augmentation_prefactor`.

`AugmentationPrefactorTable::getValueLogLinear` HOLDS the endpoint values
outside `[rho_d_min, rho_d_max]`: it mirrors
`MathLib::PiecewiseLinearInterpolation::getValue`'s `<=`/`>=` clamp branches
and its `lower_bound` interval selection exactly (both read in
`ProcessLib/RichardsMechanics/PotentialExchangeParameters.h` and
`MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.cpp`), so K is
still clamped at the table range ends — no extrapolation. The 2026-08-26
switch changed the INTERIOR chord only; extrapolation is untouched.

[HISTORICAL — the first-cut (2026-06-11) wording of the bullet and of the
clamp paragraph above, restored here VERBATIM on 2026-08-31. It had been
rewritten in place by the 2026-08-31 doc-correction pass (`AGENTS.md`
entry 31) rather than annotated, so that entry's claim that all five stale
passages were annotated `[HISTORICAL ...]` with "nothing deleted" did not
hold for this one; §6.3 wants the superseded sentence kept BESIDE the new
one, which is what this block restores. The `getValue` it describes is
still shipped — it is now the parse-time frozen-K path, no longer the LIVE
path:

--- begin restored first-cut text (byte-for-byte as at HEAD 7ec39ecf4c) ---

- live mode && table && `std::isfinite(phi)` →
  `table->getValue(rho_SR * (1 - phi))`;
- otherwise → the parse-time scalar `potential_augmentation_prefactor`.

`MathLib::PiecewiseLinearInterpolation::getValue` HOLDS the endpoint values
outside `[rho_d_min, rho_d_max]` (verified in
`MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.cpp`), so K is
clamped at the table range ends — no extrapolation.
--- end restored first-cut text ---

The "verified in" is a verification record and is restored as such: the
2026-08-31 rewrite had softened it to "read in". It still holds, re-checked
in that file on 2026-08-31 — `getValue`'s two clamp branches are
`pnt_to_interpolate <= supp_pnts_.front()` returning
`values_at_supp_pnts_[0]` and `supp_pnts_.back() <= pnt_to_interpolate`
returning the last value, unchanged.]

Threaded through every site in `RichardsMechanicsFEM-impl.h` that feeds
`potential_augmentation_prefactor` into `computeVanDerWaalsMicroPotential` /
`computeStrainedFilmState`:

- micro-potential fold point `applyFilmPressureMicroPotential` and
  `computeActiveMicroPotential`, plus the predictor/coupled
  ScalarReferenceMassStorage local solves — phi from
  `PotentialExchangeLocalSolveContext::phi` (infinity sentinel for
  default-constructed contexts → scalar fallback);
- swelling-eigenstress increment
  `computeReferenceMicroPorositySwellingStressIncrement` /
  `computeSwellingStressIncrement` — new trailing defaulted argument
  `total_porosity` (NaN sentinel → scalar), passed from
  `updateSwellingState` as the stateful total `PorosityData.phi`; one K for
  both the prev and curr Pi evaluations of an increment (mirrors the
  held-fixed p_conf telescoping convention);
- the assembleWithJacobian Maxwell p-u tangent and swelling-Jacobian blocks
  — the local total-porosity `phi` at the integration point.

Sites without any porosity in scope (the GP eigenstress-difference driver,
unit tests with default contexts) fall back to the scalar — in OFF mode this
makes every site bit-for-bit identical to before.

## Provisional choices (first cut — flagged)

- [HISTORICAL — superseded 2026-08-26 by Vinay's interpolation-scheme
  decision (commit 1bb414ac05): the LIVE path is log-linear, see
  "Interpolation scheme: LOG-LINEAR" at the end of this file. The K-linear
  `getValue()` this bullet describes survives unchanged, but only on the
  parse-time frozen-K path.] LINEAR interpolation between the calibrated
  knots (the PiecewiseLinearInterpolation shape); the interpolation shape
  between calibration anchors is UNDECIDED (memory:
  `project_dsm_k_of_dry_density`) — linear is the provisional placeholder,
  Vinay's call pending.
- [HISTORICAL — superseded 2026-06-12, see "Analytic tangent completion"
  below] The analytic dK/drho_d (= dK/dphi → strain chain) tangent
  contribution was OMITTED in the first cut: the residual saw the live K but
  the Jacobian did not carry its derivative.
- NO double counting: only the K value fed to the UNCHANGED augmentation
  law `mu_aug = sign*K*exp(-h/lambda)` varies with rho_d; no new physical
  term, no second rho_d-dependence is introduced anywhere.

## Verification (2026-06-11)

- Unit tests `RichardsMechanicsLiveKOfRhoD.*` in
  `Tests/ProcessLib/RichardsMechanics/StrainedFilmPotential.cpp`
  (structural in-test knots, not physical values): off-mode returns the
  scalar; live mode reproduces the table at knots and interior points
  (linear-interpolation identity, derived in-file); endpoint clamping.
  31/31 RichardsMechanics unit tests pass (28 prior + 3 new).
  [HISTORICAL — that live-mode linear-interpolation assertion no longer
  covers the live path. Commit 1bb414ac05 (2026-08-26) retargeted
  `LiveModeEvaluatesTable` at the table's K-linear `getValue()` and renamed
  it `TableLinearGetValueAnchorsFrozenKParsePath`, where it now pins the
  frozen-K parse path; the live path is covered by the new
  `LiveModeEvaluatesTableLogLinear`. `OffModeReturnsScalar` and
  `ClampsAtTableRangeEnds` are unchanged. The 2026-06-11 counts above are
  that day's record and are not re-asserted here.]
- Off-mode regression: `ANCHORS_MS33_ModelI/ms33_modelI_dd1400.prj` (no
  live flag) with this build → final mean sigma_zz = -4.9218 MPa, matching
  the recorded baseline (eurad-anchors run
  `2026-06-10_0841_dsm_native_h_of_eps_successful`, README_ms33.md §2).
- A live-mode production run (e.g. free-swelling Model VII with a Dixon
  K(rho_d) table) has NOT been performed yet — behavior under live K is
  PREDICTED only, pending Vinay's run plan.

## Analytic tangent completion (2026-06-12, Vinay approved)

The first-cut omission above is closed. Implemented as specified by Vinay
(2026-06-12): **Jacobian-only — the residual is untouched.**

- [HISTORICAL — this bullet specifies the tangent of the PRE-1bb414ac05
  K-LINEAR scheme; superseded 2026-08-26 by Vinay's interpolation-scheme
  decision (commit 1bb414ac05). The live tangent is no longer a per-segment
  CONSTANT: it is `getSegmentSlopeLogLinear`, proportional to the local
  K(rho_d) — see "Interpolation scheme: LOG-LINEAR" at the end of this
  file. `getSegmentSlope` itself is unchanged and remains the correct
  companion of the K-linear `getValue()` on the parse-time frozen-K path.]
  Chain: `rho_d = rho_SR*(1-phi)` => `drho_d/dphi = -rho_SR` =>
  `dK/dphi = -rho_SR * (dK/drho_d)` with `dK/drho_d` the EXACT per-segment
  slope of the piecewise-linear table — slope 0 at/outside the clamped
  edges, LEFT-segment slope at interior knots (mirrors `getValue`'s
  lower_bound interval selection). Exposed by `AugmentationPrefactorTable::
  getSegmentSlope` (`PotentialExchangeParameters.h`); NOT MathLib's
  `getDerivative`, which quadratically blends adjacent segment slopes and is
  not the derivative of the clamped value the residual actually uses.
- `effectiveAugmentationPrefactorPhiDerivative(params, phi)` returns
  `dK/dphi` [J/kg per unit phi]; 0 in every case where
  `effectiveAugmentationPrefactor` falls back to the parse-time scalar
  (mode off, no table, sentinel/NaN phi).
- mu-level K-partials: `mu_aug = sign*K*exp(-h/lambda)` is LINEAR in K, so
  `VanDerWaalsMicroPotentialData` carries exact `dmu_lR_dK` and
  `ddmu_lR_dnl_dK` (0 when the disjoining floor clamps).
- Wired into the live p-u augmentation Jacobian block in
  `RichardsMechanicsFEM-impl.h` (`assembleWithJacobian`): an additional
  `dmu_lR_dK_tot * dK/dphi * dphi/deps_v` contribution to d mu_lR/d eps_v,
  with `dphi/deps_v = (alpha - phi)/(1 + w)` derived analytically from
  `PorosityFromMassBalance` (the MS33 porosity law); other porosity laws are
  treated as strain-independent (chain 0, exact for Constant). The
  default-OFF `enable_dsm_swelling_up_jacobian` block carries a NOTE and is
  left without the chain (dead code).
- Off mode / frozen table / clamped edge: `dK/dphi == 0` -> the new block is
  skipped entirely -> bit-for-bit identical Jacobian (off-mode verified
  bitwise, see below).

### Verification (2026-06-12, measured)

- New unit tests `RichardsMechanicsLiveKOfRhoD.AnalyticPhiTangentMatchesFD
  InsideSegment` and `...ClampedEdgesAndKnots` (FD-vs-analytic, derived
  identity; scale-derived tolerances; structural in-test knots): PASS.
  Full RichardsMechanics suite: 41 tests, 39 PASS + 2 designed skips
  (build `Pi_fofnlev_20260611`, incremental).
  [HISTORICAL — `AnalyticPhiTangentMatchesFDInsideSegment` no longer exists
  under that name: commit 1bb414ac05 (2026-08-26) retargeted it at the
  table's K-linear `getSegmentSlope()` and renamed it
  `TableLinearSegmentSlopeAnchorsFrozenKParsePath`, and added
  `AnalyticPhiTangentMatchesFDLogLinear` for the live log-linear tangent.
  `AnalyticPhiTangentClampedEdgesAndKnots` is unchanged. The 2026-06-12
  counts above are that day's record and are not re-asserted here.]
- Off-mode regression: `ANCHORS_MS33_ModelI/ms33_modelI_dd1400.prj` run
  with the reference binary (`h_of_eps_20260609`, pre-tangent) vs this
  build: all 12 output VTUs bitwise-identical field data (points + every
  point/cell array).
- Live-K sanity (2-step truncated `1a_robin_A_Kl` copy,
  `task42_case1_2026-06-12/_diagnostics_1bKl/1a_robin_A_Kl_trunc2_diag.prj`):
  converges; measured Newton iterations 17/2/2 for steps 1-3, EQUAL to the
  pre-tangent full-run counts (no regression, no measured gain on this
  confined-path case).
- CURE TEST (the motivating 1b `*_Kl` step-1 singularity,
  `task42_case1_2026-06-12/_diagnostics_1bKl/README_DIAG.md`): NOT CURED.
  Both `1b_A_Kl` and `1b_B_Kl` still fail in time step #1 ("time stepper
  cannot reduce the time step size further"). Measured CHANGE in the
  Newton trajectory, though: with the tangent, |dx|_uz CONTRACTS for 7
  iterations (14.6 -> 11.4 -> 17.7 -> 14.0 -> 5.9 -> 2.4 -> 1.09) before
  blowing up at iteration 8 (168 -> 4.3e3 -> 5.1e4 -> ...), vs the
  pre-tangent monotonic divergence (14.6 -> 26 -> 44 -> 64 -> inf). The
  missing-tangent hypothesis is therefore PARTIALLY supported (basin
  improved) but the 1b compliant-top failure has an additional, still
  unidentified mechanism — evidence in `out_1b_{A,B}_Kl/run.log`
  (2026-06-12 17:04/17:05 runs); no further patching without Vinay's call.

## Interpolation scheme: LOG-LINEAR (Vinay 2026-08-26, commit 1bb414ac05)

CURRENT SPECIFICATION of the LIVE path. It supersedes every statement above
that describes K-linear interpolation between knots or a per-segment-CONSTANT
tangent; each of those is annotated `[HISTORICAL ...]` in place and describes
the pre-1bb414ac05 scheme.

**Value.** `effectiveAugmentationPrefactor(params, phi)` evaluates the table
through `AugmentationPrefactorTable::getValueLogLinear`:

    x    = rho_d = rho_SR*(1 - phi)                 [kg/m^3]
    t    = (x - x_l)/(x_r - x_l)                    [-]
    K(x) = K_l * exp(t * ln(K_r/K_l))               [J/kg]

with `(x_l,K_l)`, `(x_r,K_r)` the bracketing knots: ln(K) is linear in rho_d
between them. Node-preserving by construction (t=0 -> K_l, t=1 -> K_r), so
the calibrated knot K values are not altered by the scheme (the
floating-point knot round-trip is audited in `AGENTS.md` entry 24). Outside
`[rho_d_min, rho_d_max]` the flat endpoint hold of `getValue` is mirrored
exactly — the decision did not change extrapolation.

[CORRECTED 2026-08-31 — "node-preserving BY CONSTRUCTION" is a statement
about the ALGEBRA, not about the floating-point evaluation, and since the
uncommitted code-review batch of 2026-08-31 it is no longer where node
preservation comes from. The 2026-08-26 sentence above is kept as that
date's wording; the shipped guarantee now reads:

Node preservation at a knot is enforced by an explicit BRANCH.
`AugmentationPrefactorTable::logLinearValueOnSegment()` returns the STORED
knot value at `t == 1` (and at `t == 0`) instead of evaluating
`K_l*exp(t*ln(K_r/K_l))` there. Without that branch an INTERIOR knot — which
`lower_bound` places in the LEFT segment, at `t == 1` — goes through the
log/exp round-trip and can land ~1 ULP off `K_r`; a BOUNDARY knot is
unaffected either way, because it returns through the endpoint-hold clamp,
which contains no logarithm. The miss is not hypothetical: on the SUPERSEDED
dd900 knot pair (`K(900) = 4367.227700212952` J/kg → `K(1400) = 46000.0`
J/kg; the superseded value is kept on record in the provenance block of
`Tests/Data/RichardsMechanics/ANCHORS_MS33_ModelI/ms33_modelI_dd900.prj` and
in commit `642a8f867a`) the 1400 knot comes back as `45999.999999999993`,
rel `-1.58e-16` — measured this session, 2026-08-31, and the same miss the
2026-08-26 landing audit recorded in `AGENTS.md` entry 24 as "1400 low by
1 ULP, rel 1.6e-16". On the table shipped at `7ec39ecf4c` all four knots
round-trip bit-exactly (also measured 2026-08-31), so what the branch
removes is a LATENT defect, not a live one, and no shipped result moves.
Entry 24's ULP sentence therefore describes the PRE-cascade table and stands
as that record. Covered by the unit test
`InteriorKnotBitExactWhereRoundTripMisses`
(`Tests/ProcessLib/RichardsMechanics/StrainedFilmPotential.cpp`), added in
the same batch on a knot pair whose round-trip actually misses — the two
tables' own structural knots do not. See `AGENTS.md` entries 32 and 34.]

**Tangent.** `effectiveAugmentationPrefactorPhiDerivative(params, phi)`
evaluates `AugmentationPrefactorTable::getSegmentSlopeLogLinear`, the exact
derivative of the value above:

    dK/drho_d = K(rho_d) * ln(K_r/K_l) / (x_r - x_l)   [(J/kg)/(kg/m^3)]
    dK/dphi   = -rho_SR * dK/drho_d                    [J/kg per unit phi]

so the live tangent is PROPORTIONAL TO THE LOCAL K(rho_d), not the
per-segment constant `getSegmentSlope` returns. The one-sided/clamp
convention is unchanged: 0 at and outside the boundary knots, LEFT-segment
value at an interior knot (`idx = lower_bound - 1`). Value and tangent
therefore come from the same interpolant — which is why the tangent had to
move together with the value.

**Wiring.** Nothing in the Jacobian assembly changed shape: both live-K
blocks in `RichardsMechanicsFEM-impl.h` call
`effectiveAugmentationPrefactorPhiDerivative` and pick up the new tangent
automatically — the p-u augmentation block (`dK_dphi_pu`) and the
2026-06-14 swelling-eigenstress block (`dK_dphi_sw`, not mentioned in
"Analytic tangent completion" above, which predates it). Both still
short-circuit on `dK/dphi == 0` (mode off, no table, sentinel phi, clamped
edge).

**Unchanged by this decision.** `getValue()` / `getSegmentSlope()` (the
K-linear pair) are deliberately retained: the parse-time frozen-K resolution
in `CreateRichardsMechanicsProcess.cpp` (`table->getValue(*dry_density)`,
which is also live mode's phi-less fallback scalar) still calls them, so
non-live PRJs take the K-linear value exactly as before. And
`mu_aug = sign*K*exp(-h/lambda)` is still LINEAR IN K, so the `dmu_lR_dK` /
`ddmu_lR_dnl_dK` partials described in "Analytic tangent completion" above
stand as written — only the K(rho_d) chord and its dK/drho_d changed.

[CORRECTED 2026-08-31 — "still calls them" is right about `getValue()` and
wrong about `getSegmentSlope()`; the sentence above is kept as the
2026-08-26 wording. Call-site census, re-grepped over `ProcessLib/` and
`Tests/ProcessLib/` on 2026-08-31 (the same census is recorded in the
`AugmentationPrefactorTable` class comment): `getValue()` on this table has
exactly ONE production caller, the parse-time frozen-K resolution in
`CreateRichardsMechanicsProcess.cpp` — which is also live mode's phi-less
fallback, so the sentence above is right that non-live PRJs are untouched;
`getSegmentSlope()` has NONE. That is structural rather than an oversight:
the frozen-K path resolves K to a scalar at parse time, so it contributes no
Jacobian term and can never want a slope. `getSegmentSlope()` is retained as
that pair's exact companion and is exercised only by
`Tests/ProcessLib/RichardsMechanics/StrainedFilmPotential.cpp`.]

**Parse-time precondition on `<prefactors>` (added 2026-08-31; UNCOMMITTED
working-tree state, `AGENTS.md` entry 33).** Value and tangent both evaluate
`ln(K_r/K_l)`, so every table entry must be strictly positive.
`parsePotentialExchangeParameters` (`CreateRichardsMechanicsProcess.cpp`)
enforces that at parse time and `OGS_FATAL`s on the first non-positive
entry, naming its index and its value. The non-positive cases do NOT share
one failure mode — probed 2026-08-31 on a standalone transcription of the
kernel (IEEE-754 double), and written out case by case in the source comment
at the guard: `K_l == 0` gives `+inf` out of the logarithm and a NaN only one
step later, at `0*inf`; `K_l < 0` is the only case the logarithm itself
catches; `K_r == 0` returns a clean, plausible value of `0` while poisoning
ONLY the companion slope — the silent case, and therefore the dangerous one;
and two knots of the same negative sign produce no NaN anywhere and would
carry a negative K straight through. The guard is not restricted to live
mode, so it widens parser rejection for frozen-K decks too; all seven
shipped tables are strictly positive, so no deck in the tree changes
behaviour. PRJ-authoring consequence: the `<prefactors>` bullet under "PRJ
interface" above.

Vinay's physical justification (Dixon 2023's exponential swelling-pressure-
vs-EMDD law; the convexity/Jensen argument against the straight chord) and
the node-preservation audit are recorded in commit 1bb414ac05 and in
`AGENTS.md` entry 24 (2026-08-26). Downstream MS33 consequences are tracked
there, not here.
