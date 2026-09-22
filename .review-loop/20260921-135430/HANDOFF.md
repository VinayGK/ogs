# Handoff — review–fix loop on `dsm_upstream_mr_2026-09-20`

Written 2026-09-21 on the laptop, to be picked up on the mac mini.
Run id `20260921-135430`. Full trace: `ledger.jsonl`; regenerated summary:
`report.md`; compilation notes: `BUILD-NOTES.md`.

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

**Full detail is in `BUILD-NOTES.md`; read that before configuring anything.**
Summary:

**The branch does not build unaided on macOS 27 / SDK 27.** `range-v3` (pinned at
c704) references `_LIBCPP_TEMPLATE_VIS`, a visibility macro that libc++ no longer
defines, so every translation unit that includes it fails. The failing targets were
the `BaseLib` / `GeoLib` / `MathLib` / `NumLib` PCHs and `MaterialLib_Utils` — none
of them in the branch diff, which is confined to `ProcessLib/RichardsMechanics`,
`Tests`, `Tests/Data` and `web`. So this is a **pre-existing toolchain/dependency
incompatibility, not branch-introduced**; the older build dirs only still work
because their objects predate the SDK bump.

One build-dir-only flag fixes it. **No source file was touched, so it does not
travel with the branch** — reapply it on the mac mini *only if that host needs it*:

```bash
cmake -S <worktree> -B <builddir> -G Ninja \
  -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=ON \
  -DOGS_USE_MFRONT=ON -DOGS_BUILD_TESTING=ON -DOGS_BUILD_UTILS=ON \
  -DCPM_SOURCE_CACHE=<cpm cache> \
  -DCMAKE_CXX_FLAGS="-D_LIBCPP_TEMPLATE_VIS="
```

`OGS_BUILD_TESTING=ON` is required — the branch adds unit tests, and the older DSM
build dirs have it OFF. Try configuring **without** `CMAKE_CXX_FLAGS` first; if the
mini's SDK still defines the macro it builds unaided and the flag should not be
carried over. `BUILD-NOTES.md` §7 has the one-line check.

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
- The `-D_LIBCPP_TEMPLATE_VIS=` workaround means the binary under test was built
  with one visibility macro emptied relative to an unaffected toolchain. It is the
  minimum needed to compile on this host, but it is a deviation and is recorded as
  one in `BUILD-NOTES.md`.
- The loop has not converged. Round 2 still produced 14 new findings, so more
  remain to be found; 8 rounds of budget are left.

---

# Session 2 — mac mini, 2026-09-21/22

Rounds 3–10 run here; the 10-round budget is now spent. Run id and ledger are the
same; nothing was re-initialised.

| | |
|---|---|
| Host | `macmini.fritz.box`, macOS 27.0, Apple clang 21.0.0, SDK 27.0, 10 cores |
| Worktree | `~/git/ogs-worktrees/dsm_mr_s2_20260921` (fresh; local branch `mr_s2_2026-09-21`) |
| Build dir | `~/git/build/dsm_mr_s2_20260921`, `-D_LIBCPP_TEMPLATE_VIS=` (the SDK-27 check in BUILD-NOTES §7 fires on this host too) |
| Tip at close | `c268ac583d` — 25 commits over `c41ada5f75` |
| Rounds run | 3, 4, 5, 6, 7, 8, 9, 10 |
| Stop reason | **budget exhausted** — not converged |
| Verification | green in **every** round: 46 tests / 4 suites, 44 passed, 2 pre-existing skips, 0 failed |

## Rounds at a glance

| round | new | fixed | verify |
|---|---|---|---|
| 3 | 5 | 0 | pass |
| 4 | 8 | 6 | pass |
| 5 | 2 | 1 | pass |
| 6 | 6 | 6 | pass |
| 7 | 7 | 4 | pass |
| 8 | 9 | 5 | pass |
| 9 | 5 | 1 | pass |
| 10 | 1 | 1 | pass |

Session 2 added **14 fixup commits**. Every one is comment, prose, formatting,
`const`, or rename only — **no expression, literal, unit value, tolerance or
expected value was changed by any round.** Verified at the gate: across all of
session 2's commits, not one floating-point literal was removed or altered; the
only numbers added are §4.3 magnitude figures, and all of them sit inside
comments (checked: zero on a non-comment line).

## What session 2 covered that session 1 did not

- **The upstream `ogs-ai` review criteria** (`git@gitlab.opengeosys.org:ogs/inf/ogs-ai.git`
  @ `f5592e1`), folded in at Vinay's instruction on 2026-09-22 and applied in
  rounds 8–10. Treated as additive: it opened an axis, it did not retract any
  earlier finding, and its narrower scope rules (C++ only, `+` lines only) were
  applied to *new* findings only. Biggest yields: **183 added lines carrying
  non-ASCII** in comments (now transliterated; zero were in string literals),
  **Doxygen-invisible `//` documentation** on 9 new public functions, the
  **`OD` → `output_data`** rename, and 15 clang-format violations.
- The N/A sections were *verified*, not assumed: PETSc/MPI, OpenMP, exprtk,
  autocheck and staggered/monolithic are all genuinely absent from the diff.
  `IterationNumberBasedTimeStepping` **does** apply (beacon_1a01 inflow and
  stressprobe) and was checked line by line — no `1.0` multiplier.
- Rounds repeatedly audited their predecessors. Round 6 found three defects the
  loop had introduced itself (including a magnitude claim that did not
  reproduce); round 5 found that round 1's own fixup had leaked internal
  `CLAUDE.md` references into source destined for upstream; round 10 found round
  9's "const-correctness: zero findings" was incomplete and fixed 9 real cases.

## Open items you inherit

1. **q1–q24, all unanswered.** Session 1 left ten, session 2 added fourteen.
   Every one is guardrail-blocked — a formulation call, a cited parameter, an
   assertion value, a tolerance, a breaking config rename, or a precedence call
   between `CLAUDE.md` and `ogs-ai`. **No further round can clear them.** Full
   text per question is in `report.md`.
2. **The autosquash does NOT replay clean** — this supersedes session 1's note
   that it did (true then, before eleven more fixups landed). Rehearsed in a
   throwaway worktree: it stops on conflict **four times**, first on
   `b1f170c9b9` in `StrainedFilmPotential.cpp`, because round 4 and round 6 edit
   the same comment lines. Worse, resolving mechanically by taking the tip's
   whole-file content *does* reach an empty tree diff but **corrupts the commit
   structure** — three commits instead of four, with contents under the wrong
   messages. So an empty diff alone is not proof of a good squash. Per-hunk
   resolution is needed and is **yours**: which comment wording survives, and
   which commit each hunk belongs to, is an authorship call (§7, §9). The branch
   is therefore pushed with fixups **unsquashed**, which is a valid MR state.
3. **Two axes no round ever covered.** Integration/benchmark **ctests were never
   run** — verification was the unit binary only, so nothing here speaks to
   benchmark or reference-VTU status. And the five **beacon reference VTUs were
   never regenerated** or checked against a live run; they ship in the same MR as
   the decks that produce them.

## Harness note — read before resuming

The `review-fix-loop` skill and its `scripts/ledger.py` **do not exist on this
host**; they were never pushed off the laptop. Session 2 drove the loop manually
against this same append-only ledger using `ledger_s2.py`, committed next to it.
That reconstruction was validated adversarially: regenerating `report.md` from
the untouched round-1/2 ledger reproduced the laptop's file **byte for byte**.
Use it as `python3 ledger_s2.py --run <run dir> report|append|answer`. It has no
`init` subcommand by design.

One incident worth carrying: an interrupted agent left the loop worktree
**suspended mid-interactive-rebase**, having started the squash rehearsal in
place. It was recovered with `git rebase --quit` (not `--abort`, which would have
discarded two rounds of work) after tagging every sha as `rescue/*`. **The
rehearsal belongs in a throwaway worktree, always.**
