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

> **PARTLY SUPERSEDED (annotated 2026-09-23).** This section describes the state at
> the close of the loop (tip `c268ac583d`). Post-loop work on 2026-09-22 changed the
> tip and falsified five of its claims. They are left as written, per CLAUDE.md §6.3,
> and corrected here. **For the current state read "Session 2, continued" at the end.**
>
> 1. *"Tip at close `c268ac583d` — 25 commits"* — true at loop close; the tip is now
>    `e7a6cfe91d`, 30 commits over `c41ada5f75`.
> 2. *Round 9 "new = 5"* — the ledger holds **six** findings for round 9; its own
>    `round_end` undercounted and `report.md` copied it. Run total is 67.
> 3. *"Every one is comment, prose, formatting, const, or rename only"* — too narrow.
>    No numeric literal, unit, tolerance or expected value changed; that part holds.
>    But round 6's `b1f170c9b9` changed **four user-facing `OGS_FATAL` string literals**
>    (removing internal `DSM/*.md` references from runtime error text), and two
>    parameters were anonymised. Executable code, not comments. Nothing asserts on
>    that text.
> 4. *Open item 1, "every one is guardrail-blocked"* — overstated for about ten of the
>    24. And three are now closed (q14, q16(a), q22).
> 5. *Open item 2, the autosquash mechanism* — the four conflict stops are real and
>    reproduced independently, but the stated cause ("round 4 and round 6 edit the same
>    comment lines") is wrong: round 4's fixup never touches `StrainedFilmPotential.cpp`.
>    *Open item 3, "ctests were never run"* — no longer true; they were run post-loop
>    and pass.

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


---

# Session 2, continued — mac mini, 2026-09-22 (post-loop)

**This is the current state.** The loop's budget was spent at round 10; everything
below was done afterwards, each code change explicitly approved by Vinay.

| | |
|---|---|
| Tip | **`e7a6cfe91d`** — 30 commits over `c41ada5f75`, identical on `origin` and `github` |
| Carrier branch | `review-loop/20260921-dsm_upstream_mr`, holds this trace |
| Unit suite | 46 tests / 4 suites, 44 passed, 2 pre-existing skips, 0 failed |
| Beacon ctests | 9/9 (13 before q14 removed four duplicates), all vtkdiff gates green |
| RichardsMechanics ctests | 39/39 by label (43 before q14) |
| Build | clean, zero warnings |

## What changed after the loop

1. **Integration ctests run for the first time.** All RichardsMechanics tests pass,
   including the ~30 pre-existing registrations that exercise the shared assembler this
   branch rewrites (135 deleted lines in 7 pre-existing files) - the axis with a measured
   prior of breakage in `incident_dsm_noninertness_clamp_2026-09-20`. All four vtkdiff
   gates pass against the checked-in reference VTUs, so those VTUs reproduce from a live
   run at tip. **Trap for the next person:** a build dir made with `ninja ogs testrunner`
   alone lacks `vtkdiff` and the MFront behaviour libraries, so ctest reports failures
   that look like broken benchmarks and are missing binaries. Build those targets first.
   One unrelated failure remains, `ThermoRichardsMechanics/.../bentonite_column-LARGE-omp`:
   `large`-labelled (excluded by the documented `ctest -LE large` gate) and structurally
   outside this branch - ThermoRichardsMechanics has zero references to RichardsMechanics.
   **Not** labelled pre-existing: no base build was run to prove it.
2. **21 clang-tidy `WarningsAsErrors` hits fixed** (`aa78e267c6`). Both checks are in
   the repo's own `.clang-tidy` `WarningsAsErrors`, so they were `error:` and would have
   blocked CI. No round found them because no round ran the tool. 17 boolean guards plus
   4 unnamed `ioName` parameters. **De Morgan was deliberately NOT applied:** clang-tidy's
   own fixit rewrites `!(a && b)` as `!a || !b`, which changes NaN semantics, and five of
   the guards are NaN-reachable with no leading `isfinite` (a bisection bracket, a branch
   selector, a PRJ validation guarding an `OGS_FATAL`, the macro-floor cutoff, the 2x2 FD
   denominators). Each predicate was hoisted into a named `bool`, tokens unchanged.
   17/17 equivalence proven by mechanically un-hoisting and token-diffing against base.
3. **q22** (`28913af8a6`) - two `EXPECT_NEAR` replaced by `EXPECT_DOUBLE_EQ`; removes
   the tolerance literals 1e-14 and 1e-12, after re-verifying both exact-zero contracts.
4. **q16(a)** (`0a82f54e06`) - `std::call_once` WARN at the two silent Newton exits;
   no number or control flow touched. The policy half stays open.
5. **q14** (`bc2dac021d`) - four duplicate beacon ctest registrations removed.
6. **nu = 0.2 documented** in all five beacon deck headers (`e7a6cfe91d`), header text
   only. Source: BEACON D3.3, annex "D3.3 - BGR Results of Task 3.3", PDF p.66, Table 2.
   D3.3 relays it from Akesson et al. 2010, which **D3.3's own reference list (p.72)
   identifies as SKB TR-10-44** - not TR-10-11. A local text extraction of TR-10-44
   gives nu = 0.2 in Table 12-1 "Data used for different MX-80 buffer materials", printed
   p.81. It matches the live value, so it documents and changes nothing.
   **E could not be documented:** the decks hold 50e6 Pa (40e6 for 1b) while the two local
   Beacon sources give 180 MPa (D3.3 p.66) and 28.1 MPa (D3.3 p.214 at e0). They bracket
   the shipped value without supporting it. E stays explicitly uncited in every header.

Also produced, outside the repo (durable, mac mini only):
`~/ogs-models/scratch/2026-09-22_1130_q1-q24-decision-dossier_successful/` - a
per-question decision dossier, the refutation audit, and the correctly-paired refutations;
`~/ogs-models/scratch/2026-09-22_1047_ogs-ai-review-criteria_successful/` - a mirror of the
upstream `ogs-ai` criteria at `f5592e1`.

## State of q1-q24

| status | questions |
|---|---|
| **Closed** | q14, q16(a), q22 |
| **Ordinary calls - need only a preference** | q13, q18, q19, q20, q21, q24 |
| **Genuinely yours - physics, cited value, or tolerance** | q1, q2, q7, q8, q9, q10, q11, q12, q15, q17, q23 |
| **Blocked on absent literature** | q5, q6, and parts of q3/q4 |

Suggested order, from the dossier's critic: **q2** first (the only one that moves a
shipped number; it gates q17, q23 and the inflow deck's header). Then **one ruling:
are the five decks meant to exercise the augmentation at all?** K = 0 in all five, and
three set A = 1e-30, Sa = 1.0 - deliberate null cases. That single answer closes or
dissolves q17 and q23 and reframes q3-q6. Then **q11**, the only one the §6.7.5 gate
blocks a push on. q1, q12 and q15 are one ruling: which mu the mechanical partner is
conjugate to.

**The E/nu fork, if you take it up:** E = 28.11 MPa is only self-consistent coupled with
nu = 0.35 from the same LEI table - that changes two literals and rebaselines registered
ctests. The alternative, nu = 0.2 alone, is already done and costs nothing. They are
mutually exclusive.

## Citation refutation audit

Eight of ten citation candidates for q3-q6 were refuted; the refutations were then
themselves audited (§13.3(5)). **At least three kills were false**, and the cause was the
orchestrator's instruction to "default to REFUTED when uncertain", which led refuters to
conflate "adopting this would change a deck literal" with "the citation is invalid".
Firm: the E kill was false (it put the deck's nu = 0.2 into the source's own formula;
with the source's nu = 0.35, E_min = 18.0 MPa and the floor does not bind at e0), and both
nu kills were false. The per-file auditors and the synthesis **disagree** on the two
permeability candidates; both verdicts are recorded and neither is authoritative. q6 is
**blocked, not dead**.

## Errors made by the orchestrator in session 2 - all recorded in the ledger

- **A pairing bug** built audit files by parallel-completion order rather than
  submission order, gluing each candidate to another verdict's reasoning. It invalidated
  one whole audit pass, which was discarded and re-run.
- **Three false premises in agent briefs**, one failure mode (§13.3(3)): "Akesson et al.
  2010 is TR-10-11" (it is TR-10-44); "TR-10-11 is not on this machine" (a `find
  -maxdepth 6` missed it at depth 9); "TR-10-44 is not on this machine" (its extraction
  is named `SKB_SRSite_data_Akesson2010.txt`, with no "10-44" to match). The first reached
  all five headers before an adversarial audit caught it; that commit was replaced, never
  pushed.
- **A verifier told to change no file amended a commit** after the adversarial re-audit
  had passed a different sha. HEAD was re-verified against the primary source directly.
- **An agent installed Homebrew LLVM** to get clang-tidy without asking first. Vinay
  later confirmed that is fine.

## Still open

- **Autosquash** needs per-hunk resolution by hand - it stops four times, and whole-file
  resolution collapses four commits into three while still showing an empty tree diff.
- **The 6 `misc-const-correctness` findings** were not approved, so they stay unapplied.
- **The MBP has not been synced.** It went offline during session 2. Both remotes carry
  everything; one `git fetch` brings it across.
- **CLAUDE.md §12.1 family 8 names TR-10-11**, but D3.3's Akesson attribution is
  TR-10-44. Any other header routing a D3.3 Akesson claim to family 8 has the defect fixed
  here. Family 8's recorded path also omits `subtask41/`.
- **D3.3 lists V. Kumar among its contributors** - the annex the nu citation rests on is
  partly the repo owner's own work. Whether D3.3 counts as §12.1 family 1 is left as
  Vinay's call in the headers, not asserted.
- **To read on the MacBook:** Talandier D5.1.1 (2018) settles q5/q6 and the 1a01
  parameter set. The original TR-10-44 PDF would confirm the Table 12-1 row, which rests
  here on a single text extraction.
