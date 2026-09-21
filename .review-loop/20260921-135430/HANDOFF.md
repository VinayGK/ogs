# Handoff — review–fix loop on `dsm_upstream_mr_2026-09-20`

Written 2026-09-21 on the laptop, to be picked up on the mac mini.
Run id `20260921-135430`. Full trace: `ledger.jsonl`; regenerated summary: `report.md`.

---

## Where the run stands

| | |
|---|---|
| Branch | `dsm_upstream_mr_2026-09-20` |
| Fixed point | `c41ada5f75` (the upstream master merge the branch sits on) |
| Scope | 70 files, +14073 / −135 |
| Tip at pause | `dea92ed01cadee4ee6550143ff443ea5a2caa8c2` |
| Rounds done | 2 of a 10-round budget — **8 rounds remain** |
| Last verdict | `continue` (round 2: new=14 fixed=4 open=0 verify=pass) |
| Verification | green both rounds: 46 tests, 44 pass, 2 pre-existing skips, 0 fail |
| Tree | clean |

Round 2 was itself a resume: the first round-2 agent died on a transport error
(ECONNRESET) after triage but before verifying. A fresh agent picked the ledger up,
recorded the four missing `resolve` events, verified, and closed the round. Nothing
was lost — that is what the append-only ledger is for.

## What the loop changed

Five `fixup!` commits on top of the four original branch commits:

```
dea92ed01c fixup! RichardsMechanics: add micro-macro potential exchange
c3bb85feb6 fixup! RichardsMechanics: add micro-macro potential exchange
7a9c21e206 fixup! web: benchmark page for the double-structure model
dfc61bc1c7 fixup! web: benchmark page for the double-structure model
bb10d046b0 fixup! RichardsMechanics: add micro-macro potential exchange
2680878e57 web: benchmark page for the double-structure model      <- original tip
86fd36c4ff Tests/Data: BEACON 1a01/1b/1c double-structure benchmarks
4f0541182c Tests: unit tests for the double-structure micro-macro exchange
569ba9dd15 RichardsMechanics: add micro-macro potential exchange
```

All seven applied fixes are **comment / prose / unit-annotation only** — no
expression, literal or parameter changed. The round-2 agent rehearsed the autosquash
in a throwaway worktree: `rebase --autosquash` onto `c41ada5f75` replays clean and
`git diff dea92ed01c HEAD` is empty, leaving the four original commits.

When you are ready to fold them in:

```bash
git -C <worktree> rebase --autosquash --interactive c41ada5f75
```

## Verification environment — read this before rebuilding

**The branch does not build out of the box on macOS 26.** `range-v3` (pinned at
c704) uses `_LIBCPP_TEMPLATE_VIS`, which the macOS 26 libc++ removed, so every
translation unit that includes it fails. The failing targets were the
`BaseLib` / `GeoLib` / `MathLib` / `NumLib` PCHs and `MaterialLib_Utils` — none of
them in the branch diff, which is confined to `ProcessLib/RichardsMechanics`,
`Tests`, `Tests/Data` and `web`. So this is a **pre-existing toolchain/dependency
incompatibility, not branch-introduced**; the older build dirs only still work
because their objects predate the SDK bump.

Worked around with build-dir-only compile definitions. **No source file was
touched, so this workaround does not travel with the branch** — you must
reapply it on the mac mini:

```bash
cmake -S <worktree> -B <builddir> -G Ninja \
  -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=ON \
  -DOGS_USE_MFRONT=ON -DOGS_BUILD_TESTING=ON -DOGS_BUILD_UTILS=ON \
  -DCPM_SOURCE_CACHE=<cpm cache> \
  -DCMAKE_CXX_FLAGS="-D_LIBCPP_TEMPLATE_VIS= -D_LIBCPP_ENUM_VIS= -D_LIBCPP_EXCEPTION_ABI= -D_LIBCPP_HIDDEN= -D_LIBCPP_FUNC_VIS= -D_LIBCPP_TYPE_VIS="
```

If the mac mini runs an older SDK that still defines these macros, drop
`CMAKE_CXX_FLAGS` and check whether it builds unaided — worth knowing either way.

Verification command (adjust the build dir):

```bash
cd <builddir> && ninja -j <cores> ogs testrunner && ./bin/testrunner --gtest_filter='RichardsMechanics*'
```

Expected: 46 tests from 4 suites, 44 passed, 2 skipped
(`RichardsMechanicsExactFilmPair.LiquidCarrierEnergyPressureConsistency`,
`.ExpulsionProbeDrainedRamp`), 0 failed. Both skips are pre-existing and carry
stated reasons. **Integration ctests are not part of the verify command** — the
loop never ran them.

## To resume the loop

```
/review-fix-loop on dsm_upstream_mr_2026-09-20 maxmimum 10 iterations
```

then tell it the run dir already exists and to continue at **round 3**:
`.review-loop/20260921-135430`. Do not re-`init` — that starts a second run and
orphans this trace.

## The ten parked questions — all yours, none of them the loop's to answer

These are the run's real open item. Every one is guardrail-blocked: the loop
confirmed the defect, could not repair it without making a call CLAUDE.md reserves
to you, and left the artifact honest rather than patched. **No further round can
clear these.** Full reasoning per question is in `report.md` §"Questions for you".

**q1 — the Pi/Pi′/Pi″ formulation call.** `applyMacroFloorCutoff`
(`RichardsMechanicsFEM-impl.h:693-697`) rescales `mu_lR`, `dmu_lR_dnl`,
`dmu_lR_drho_lR`, `dmu_lR_dnS`, `dmu_lR_drho_SR` by `g_cut` but never
`d2mu_lR_dnl2`, which is then consumed uncut at :880 as Pi″. The block comment at
:866-870 says Pi, Pi′, Pi″ are the *bare* van-der-Waals disjoining pressure and its
derivatives — yet Pi and Pi′ are built from the already-cut values. Code and comment
disagree in both directions, so there are two self-consistent repairs: make them all
bare, or make them all cut (adding `g_cut*d2 + 2*dg_dnl*dmu0_dnl + d2g_dnl2*mu0`,
with `d2g/dn_l² = -2/w²` inside the band). Which `mu` the mechanical partner is
conjugate to is physics, not chain-rule completion.
*Inert in every shipped deck* — all five registered beacon PRJs set
`macro_porosity_floor = 0.0`, and the cutoff returns early unless floor > 0.

**q2 — `specific_surface` m²/g vs m²/kg.** The same incident already recorded in
CLAUDE.md §1.3. `beacon_1a01_dsm_micromacro_inflow.prj:123` cites 523 m²/g
(Seiphoori, Ferrari & Laloui 2014, Tab. 1 p. 724) and self-discloses in the same
header that the code consumes it as m²/kg, entered unconverted, so h is a factor
1000 too large — and it ships in a registered ctest. Both repairs are yours:
523 → 523000 changes a cited material parameter (§1.1/§12.5), or re-annotating the
code to m²/g changes the consumed convention for every DSM deck in the tree.

**q3–q6 — four registered ctests with uncited §12.2 provenance.**
`beacon_1a01_..._smoke.prj`, `..._stressprobe.prj`, `beacon_1b_..._smoke.prj`,
`beacon_1c_..._smoke.prj`. Each declares material-parameter groups (elastic E/ν,
porosity split φ0/φ_tr, intrinsic permeability) as an explicitly uncited "working
value", which §12.2 calls non-compliant for a deck in `Tests.cmake`. The only
compliant repairs are supplying a citation from one of the eight §12.1 families
(§1.1 forbids the loop inventing one) or dropping the ctest registration (which
deletes the coverage the branch exists to add). The headers are honest — they
declare the gap rather than hide it.

**q7 — a tolerance gate that does not cover its own observation.**
`ProcessLib/RichardsMechanics/Tests.cmake:177`: the beacon_1c σ vtkdiff gate is
registered at abs 1e-11 while the adjacent comment records an observed 4.91e-11 on
another toolchain, root cause unresolved. §1.2 makes raw tolerance literals
user-approved; §9 routes a tolerance not derived from problem scale to ASK USER.

**q8 — two "BaselineHistory" tests that never regress against their baseline.**
`DSMMicroMacroSingleIntegrationPoint.cpp`, tests at :1602 and :~1767. They load
per-step reference rows (`n_l`, `mu_LR`, `mu_lR`, `rho_lR`, `rho_l_hat`,
`epsilon_sw`, `stress_xx`) from two CSVs, but exhaustive grep confirms the only
fields ever read are `row.pressure`, `row.saturation`, `row.epsilon_v_total`,
`row.delta_epsilon_v`, `row.phi`. Every assertion is a finiteness/positivity/range
check or the trivial `phi_m + phi_M == phi` identity. A regression that flipped the
sign or magnitude of the computed `n_l`, `mu_lR`, `rho_lR`, `rho_l_hat` or swelling
stress would pass both tests. The same file's `...ReferencePath` test (:665) shows
the correct `EXPECT_NEAR`-against-independent-reference pattern. The repair needs
expected values and tolerances — §1.1 hard-bans the loop writing those, §3 MUST-2
requires `TODO(user)`.

**q9 — mismatched chain rule in the analytic Jacobian (tangent-only).**
`dmu_lR_vdw_drho_lR` sums `micro_potential.dmu_lR_drho_lR` (w.r.t. bulk `rho_LR`)
and the Maxwell partner's `mc.dmu_lR_mech_drho_lR` (w.r.t. the confined micro
density `rho_mc`), then chains the sum through the **bulk** `drho_LR_dpL`, although
`drho_lR_exchange_input_dpL = rho_lR_state*beta_LR` is already computed and threaded
in. Is bulk pairing the intended convention for the mechanical partner too (as the
comment at :5314-5316 states for the film channel)?
*Left unchanged — residuals and the converged solution are byte-identical either
way; at most Newton converges more slowly when
`use_micro_liquid_density_for_micro_pressure` is on and `film_pressure_coupling`
is off.*

**q10 — parse-time fallback K uses the K-linear interpolant.** Every runtime site
uses `getValueLogLinear()`; the phi-less fallback resolves with the K-linear
`getValue()`. `PotentialExchangeParameters.h:21-30` says the K-linear pair is
deliberately retained for exactly this parse-time path. Confirm that is intended for
an off-knot `<dry_density>` in live mode too, where the fallback would differ from
the live value at the same ρ_d.
*Left unchanged and treated as documented design — all shipped decks sit on table
knots, where both interpolants agree exactly.*

Record answers as they come:

```bash
python3 ~/.claude/skills/review-fix-loop/scripts/ledger.py \
  --run <worktree>/.review-loop/20260921-135430 answer --qid q1 --text "..."
```

## Caveats on what "verified" means here

- Verification is **unit tests only**. Integration ctests were never run, so nothing
  in this run speaks to benchmark or reference-VTU status.
- The seven fixes are comment-only, so the green verification is a statement that
  nothing broke, not evidence that any physics was validated.
- The `_LIBCPP_*_VIS` workaround means the binary under test was built with
  slightly different visibility attributes than an unmodified toolchain would
  produce. It is the only way this branch builds on this host, but it is a
  deviation and is recorded as one.
- The loop has not converged. Round 2 still produced 14 new findings, so more
  remain to be found; 8 rounds of budget are left.
