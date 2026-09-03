# Working rules for AI coding agents on this repo

You implement, verify, and document. You do NOT make scientific or
numerical decisions. All physics, parameter values, and expected
results belong to the user (Vinay). When in doubt: ASK.

These rules apply to every agent and subagent operating in this
repo. They supersede default behavior.

---

## §0 Meta-protocol — notification and exemptions

### §0.1 You MUST notify the user when a guardrail fires

A guardrail "fires" when applying a rule from §1–§12 would change
what you would otherwise have done — typically: stopping to ASK
before writing a literal, leaving an expected value as TODO,
labeling a claim "predicted not verified," refusing to add a
trailer, or blocking a partial commit. (Default-state behaviors
like "no `Co-Authored-By` trailer" do not fire on each commit —
they only fire if you would otherwise have added the trailer.)

When a guardrail fires, you MUST tell the user explicitly,
naming the rule. Silent application is NOT acceptable.

Format:

> "Guardrail §<n.m> fires: <one-line reason>. <Action I would
>  take by default>. Approve, or grant one-off exemption?"

Examples:

> "Guardrail §1.1 fires: I need a literal for `A_Hamaker` here
>  and have no cited source. Default action: STOP and ask.
>  Approve a value with source, or grant one-off exemption?"

> "Guardrail §3 fires: you asked me to write the assertion
>  values for the new shear test. Default action: leave as
>  `TODO(user)` and ask. Approve a value with source, or grant
>  one-off exemption?"

> "Guardrail §5 fires: this commit message states 'reduces
>  iteration count by 30%' as fact, without a re-run. Default
>  action: relabel as 'expected, not yet verified.' Approve
>  the relabel, or grant one-off exemption?"

### §0.2 One-off exemptions

The user MAY grant a one-off exemption for a specific action.
Phrases that count as an exemption: "exempt", "skip the rule
here", "go ahead anyway", "override", "yes do it anyway",
"bypass §<n>", or equivalent. Ambiguous responses ("ok",
"fine", "sure") do NOT count — re-confirm before proceeding.

When an exemption is granted you MUST:

1. **Scope the exemption narrowly.** It applies only to the
   specific action discussed. Future similar actions in the
   same conversation need their own exemption — unless the
   user says "exempt all", "blanket exemption", or equivalent
   broader phrasing.

2. **Document the exemption in the artifact** where the rule
   would have fired. Format:

   ```
   GUARDRAIL EXEMPTION §<n.m> (YYYY-MM-DD): <one-line user
   justification, or "no justification provided">.
   ```

   Examples — inline code comment:

   ```cpp
   constexpr double dummy_E = 1.0e9;  // GUARDRAIL EXEMPTION
       // §1.1 (2026-05-27): smoke-test placeholder, not
       // physical; user approved.
   ```

   Example — commit message footer:

   ```
   Refactor mu_lR computation.

   GUARDRAIL EXEMPTION §5 (2026-05-27): describes expected
   5% runtime improvement as fact; user approved as the
   refactor is algorithmic with no numerical changes.
   ```

3. **Never silently exempt.** Without an explicit phrase from
   the user, the rule applies. If you skipped a rule and only
   noticed afterwards, surface this immediately ("I just
   realized §X fired and I didn't notify — here's what
   happened") rather than burying it.

4. **Exemptions do not carry across artifacts.** A §1.1
   exemption for one literal in `foo.cpp` does NOT exempt a
   different literal in `bar.cpp`, even in the same session.

5. **Exemptions do not retroactively bless past silent
   skips.** If you discover an old literal was written without
   citation or exemption, ASK now — do not annotate it as
   exempt without an explicit user confirmation.

---

## §1 Numerical literals — the hard rule

### §1.1 HARD BAN — params, expected values, BC magnitudes, signs

You MUST NOT write a literal for any of these without a citation:

- Material parameters (E, nu, A_Hamaker, K, lambda, phi, rho, alpha, ...)
- Expected values in test assertions
- Sign conventions, exponents in constitutive laws
- Boundary-condition magnitudes (sigma0, p0, T0)

The citation MUST be one of:

1. A quote from the user in this conversation.
2. A citable source with a locator:
   - paper: author + year + page/eq/table (e.g., "Villar 2007 Tab. 3")
   - textbook: title + section
   - online resource: URL + date accessed
   - reference database: name + entry (e.g., "NIST WebBook, water,
     2024-03-12")
   - vendor datasheet: vendor + part + revision/date
3. A prior commit hash traceable to user-approved work.
4. An analytical derivation in the same file.

If none exist → STOP and ASK. When proposing a value, label it as a
candidate with source, and wait for explicit approval:

> "I propose A_Hamaker = 2.2e-20 J from Israelachvili & Adams 1978.
>  Approve?"

For sources from item 2, an inline comment on the line where the value
appears MUST record the source. Examples:

```cpp
constexpr double A_Hamaker = 2.2e-20;  // J. Israelachvili & Adams 1978,
    // J. Chem. Soc. Faraday Trans. I, vol. 74, p. 975, Table 2
    // (mica-water-mica SFA).

constexpr double rho_water_293K = 998.21;  // kg/m^3. NIST WebBook,
    // water, T=293.15 K, p=101.325 kPa, retrieved 2026-05-27.
```

When multiple literals come from a single source (e.g., a full PRJ
block transcribed from a paper, or a parameter table from one
publication), a single citation comment at the top of the block is
sufficient — per-line citation is not required. The block boundary
must be unambiguous (a contiguous edit, a clearly delimited region).

> **Incident — A_Hamaker is not a knob:** A = 2.2e-20 J is a material
> constant from Israelachvili & Adams 1978 (mica-water-mica SFA), NOT
> a fitting parameter. Only K and lambda are calibrated to Villar.
> Tuning A silently would discard the literature anchor.
> See [[ogs_rm_dsm_potential_physics]].

> **Incident — micro vs bulk density:** The rho_lR in the vdW
> denominator is the *confined micro-liquid density* (~1100 kg/m^3,
> from micro EOS), NOT bulk free-water (~1000 kg/m^3). Density
> literals MUST cite which density they refer to (micro vs macro).
> See [[ogs_rm_dsm_potential_physics]].

### §1.2 Allowed — scoped numeric defaults

These may be written without explicit approval but should be flagged
in the response:

- Array sizes, loop bounds, iteration caps
- File-format constants, version numbers
- Tolerances *derived in the same file* from problem scale
  (e.g., `tol = 1e-10 * scale`)

Raw tolerance literals require user approval — see §3 tolerance
incident.

### §1.3 Declared unit vs consumed unit — the mechanical check

Approved by Vinay 2026-09-02.

Every dimensioned PRJ parameter MUST carry its unit in the
`ogs-prj-audit` registry (`audit/properties/bands.yaml`), and the audit
gate MUST fail when that declared unit disagrees with the unit the
consuming code annotates. A parameter whose unit is declared in only
one of the two places, or in neither, is non-compliant.

A dimensional-consistency comment is NOT this check. §4.2's unit
annotation verifies that exponents balance; it is blind to scale, and
balances identically under a wrong unit. §1.3 compares the two
declarations against each other. Both are required.

> **Incident — specific_surface m²/g vs m²/kg (2026-09-02):**
> `specific_surface = 523` is m²/g (Seiphoori, Ferrari & Laloui 2014,
> Géotechnique 64(9):721-734, Tab. 1 p. 724; corroborated Tang &
> Cui 2005 CGJ Tab. 3 / Saiyouri et al. 1998 at
> 522 m²/g). `PotentialExchange.h:223` annotated it `[m^2/kg]` and
> consumed it as such. The DSM van der Waals core carries Sa³ with no
> compensating λ, so the mismatch suppressed it by exactly 1e9 while
> ξ = h/λ stayed invariant and hid it. `bands.yaml:42` had already
> declared `m^2/g` — correctly — with no rule comparing it to the
> code, and with null bands so it never fired. The infrastructure
> existed and was inert. OGS's own house convention for the tag name
> is per gram
> (`Documentation/ProjectFile/prj/chemical_system/surface/site/t_specific_surface_area.md`).

---

## §2 Calibration vs verification — NEVER cross the streams

A test MUST NOT both fit parameter P AND assert on a quantity derived
from P. That is circular and hides bugs for weeks.

If a test calibrates K, its assertion must use a quantity *independent*
of K — a conservation law, an analytical limit, a symmetry, or a value
measured upstream of calibration.

When a test crosses this line: FLAG before running. Do NOT "fix" it by
adjusting the expected value (see §3, no-tuning rule).

> **Incident — active_nS denominator (2026-05-26):** A wrong
> `(1-phi_M)` instead of `(1-n_l)` survived for weeks because Model I
> calibration absorbed the error (K values shifted to fit Villar
> targets). The bug only surfaced when a different verification config
> showed an 18% drift at the dry endpoint.
> See [[ogs_rm_dsm_potential_physics]].

> **Incident — EOS bypass via small `a`:** Model I PRJs set
> `micro_liquid_density_a = 1e-16`, making the micro EOS effectively
> constant. Under that setting K calibration is invariant to many
> micro-physics changes (incl. active_nS fix) because rho_lR is
> decoupled from p_L_m. A test that passes under a=1e-16 may FAIL
> under a=50 (real EOS active). Tests MUST run under the same EOS
> regime as the physics they verify; PRJs MUST flag
> `micro_liquid_density_a` when tuned to bypass physics.
> See [[ogs_rm_dsm_potential_physics]].

---

## §3 Test generation rules

When asked to generate or suggest tests:

### MUST

1. Generate test STRUCTURE only (name, input config, what is probed).
2. Leave expected values as `// TODO(user): supply expected value` or
   propose a candidate per §1.1 and wait for approval.
3. State the physics anchor from this closed list:
   (a) analytical limit, (b) conservation law, (c) frame indifference,
   (d) symmetry, (e) published benchmark with citation, (f) regression
   baseline previously approved by user (cite commit). No other
   justification is acceptable. If none fits → ASK.
4. Flag overlap with existing benchmarks. The user's test matrix is
   authoritative — new tests are SUPPLEMENTS, never replacements.
5. If verification of a step is impossible, the test MUST use
   `GTEST_SKIP()` with a TODO comment. NEVER a trivially-passing
   assertion.

### MUST NOT

- Iterate "change expected value until test passes." A failing test
  goes back to the user with a root-cause hypothesis.
- Generate >5 tests in one batch without user confirmation.
- Claim a test is "complete" or "sufficient." Only the user does.

### Output format for proposed tests

```
Name:            <test name>
Physics anchor:  <one of the closed list above>
Input config:    <strain/stress/state/loading>
Expected:        TODO(user) — proposed: <value> from <source>
Justification:   <why this is discriminating>
Catches:         <failure mode>
Overlap:         <existing benchmark name | none>
```

> **Incident — sigma0 free-boundary (RM/DSM):** Default `sigma0=0` at
> high suction caused unphysical drift. Correct rule:
> `sigma0_eff = alpha * chi(S_w) * p_L_initial` — never defaulted.
> See [[feedback_ogs_rm_sigma0_free_boundary]].

> **Incident — pressure tolerance floor:** abstol = 5e-8 Pa at 100 MPa
> suction is below machine epsilon for the pressure scale; Newton
> stagnated, every time-step rejected. Use ≥ 1 Pa at the 100 MPa
> scale. Tolerances must be derived from problem scale.
> See [[feedback_ogs_rm_pressure_tolerance]].

---

## §4 Code-edit discipline

### §4.1 Additive accumulators

For accumulator variables (potentials, fluxes, residuals), you MUST
use `+=`. Replacement assignment (`=`) on an accumulator is FORBIDDEN
— it silently discards prior contributions.

Awareness: the "is this an accumulator?" heuristic can misfire on
single-use intermediates that look like accumulators (e.g., a local
`sum` initialized to zero, used once, then assigned out). When
unsure, check whether the variable is read elsewhere in the same
scope before final write. If still unclear → ASK rather than silently
picking `=` or `+=`.

> **Incident — overwriting vdW base:**
> ```cpp
> out.mu_lR += mu_aug;   // correct
> out.mu_lR  = mu_aug;   // FORBIDDEN; discards Hamaker base.
> ```
> Calibration would absorb the missing term silently.
> See [[ogs_rm_dsm_potential_physics]].

### §4.2 Unit discipline on physical quantities

After ANY change to a potential, flux, or stress expression, you MUST
verify units in a line comment on the changed line:

```cpp
out.mu_lR += A * Sa3 * nS3 * rhoSR3 / (6.0 * pi * nl3 * rhoLR);  // J/kg
```

If you cannot annotate the units, the change is not finished.

> **Incident — vdW term in Pa, not J/kg (commit 0d579e8aeb,
> 2026-05-22):** The vdW base term was dimensionally wrong (Pa
> instead of J/kg) for an extended period. The exchange equation
> `rho_hat_l = alpha_M_eff * (mu_LR - mu_lR)` REQUIRES J/kg on both
> sides; mixing units produced a silent dimensional error.
> See [[ogs_rm_dsm_potential_physics]].

### §4.3 Powers of unit-bearing quantities — state the magnitude

Approved by Vinay 2026-09-02.

Any expression that raises a unit-bearing quantity to a power ≥ 2 MUST
carry, in addition to §4.2's dimensional annotation, the **physical
magnitude** its inputs imply — evaluated at a named state, with the
quantity a reader can sanity-check:

```cpp
// Units: A · Sa^3 · nS^3 · rho_SR^3 / (n_l^3 · rho_lR) = J/kg  ✓  (§4.2)
// Magnitude at dd1600 endpoint (n_l = 0.4245, Sa = 5.23e5 m^2/kg):
//   h = n_l/(nS*rho_SR*Sa) = 5.07e-10 m = 0.51 nm ~ 2 water layers  (§4.3)
```

Rationale: a dimensional check cannot detect a scale error, and a power
≥ 2 amplifies one. The magnitude line makes it human-visible. A film
thickness of 507 nm written next to the word "interlayer" is caught in
one second; `[m^2/kg]^3` balancing to J/kg is not.

Where the derived magnitude has an accepted physical range, say so and
cite it. If you cannot state the magnitude, the change is not finished.

---

## §5 Documenting consequences — predicted ≠ verified

You MUST NOT write "consequences," "expected impact," or "effects" of
a code change in commits, docs, papers, or AGENTS.md as if they were
established.

- Predicted impact MUST be labeled as predicted (e.g., "expected to
  shift K by ~5%, not yet verified").
- Upgrade to "verified" ONLY after re-running the affected benchmarks
  and confirming the prediction.
- If a benchmark contradicts the prediction, update the doc AND the
  prediction reasoning.

### §5.1 Carve-out — numerical inline-comment claims

Any inline code comment making a numerical claim about runtime
behavior (iteration counts, time-step sizes, convergence rates,
memory usage) MUST label whether it is expected or measured:

```cpp
// expected: 3-5 Newton iterations per step (not benchmarked)
// measured: 3-5 Newton iterations per step (run X, commit Y)
```

This carve-out applies to *numerical* claims only. Qualitative inline
comments ("converges quickly") don't trigger it but also shouldn't be
written.

Conversational discussion is exempt from this section.

> **Incident — paper note correction:** When the active_nS fix
> landed, paper_DSM.tex initially claimed "K must be re-calibrated."
> Re-running Model I showed K values preserved to <0.05% under
> `micro_liquid_density_a=1e-16`. The note was corrected to "table
> values preserved under the fix." The original wording was a
> plausible but unverified consequence claim.
> See [[ogs_rm_dsm_potential_physics]].

---

## §6 Operational discipline

1. **Read full implementation context** before any change — call
   graph, test suite, and PRJ usage. No shortsighted edits.

2. **Never delete VTU mesh or PRJ files.** They are tracked assets.
   Use `git checkout --` immediately if accidentally removed.

3. **Never delete .md files.** If content is superseded, annotate as
   historical/deprecated within the file.
   - **They must be committed, not merely on disk.** CLAUDE.md, every
     AGENTS.md, and all guardrail/worklog `.md` files are tracked and
     committed on a durable branch (master), because an untracked or
     ignored copy is silently absent from every fresh `git worktree` and
     from `master`. Therefore: (a) never list them in `.gitignore` or
     `.git/info/exclude`; if a pattern (`CLAUDE.md`, `AGENTS.md`,
     `AGENTS_*.md`) ever excludes them, remove it and re-add the files;
     (b) if any is untracked, `git add` and commit it on master so it
     propagates to all worktrees, not just the throwaway one you are in;
     (c) verify with `git ls-files CLAUDE.md` — empty output means the
     rule is being violated right now. These files may only be *cleaned
     up* (trimmed, reflowed, annotated historical) — never removed from
     disk or index. The two occasions this file went missing (2026-05-29,
     2026-06-06) are recorded in [[feedback_research_guardrails]].

4. **Never remove agent instructions from AGENTS.md.** Mark completed
   steps as DONE with a YYYY-MM-DD timestamp and a one-line outcome.
   Historical record is intentional.

5. **Partial commits OK; no broken states.** Iterative WIP commits on
   a feature branch are allowed, but each commit MUST leave the tree
   in a *reversible, functional state*:
   - The build must compile.
   - Tests that were passing before must still pass.
   - Half-applied refactors, uncompilable code, or trees with
     newly-failing tests MUST NOT be committed.

   Failing unit tests MUST be root-caused; "pre-existing" is not an
   acceptable explanation for a new failure introduced by the commit.
   Pre-existing failures that pre-date your work and remain
   pre-existing — note them, do not claim them as yours.

6. **Resolved TODOs must be marked DONE with date.** Never leave a
   resolved item open; never mark an item DONE unless verified.

7. **Provenance & traceability gate — before EVERY commit and push.**
   This gate fires on every commit/push; failing any check is a
   FLAG, STOP (announce per §0.1) — never commit/push through it.
   Before staging, verify and report:

   1. **Header–value consistency.** For every DSM PRJ (or
      parameter-bearing artifact) in the change, each literal cited in
      its §12 provenance block / header (fitted K, target sigma,
      densities, E, ν, …) MUST equal the live `<…>` value in the same
      file. A header naming a different number than the value it
      documents is a STOP. (Incident 2026-06-10: III/IV headers cited
      K=43182 while the value and the frozen run used 103879.)

   2. **Cross-artifact consistency.** A committed value MUST agree
      across the PRJ, its calibration record (`_calib_result*.json`),
      and any report / manifest / metrics that quote it. A value that
      differs between the PRJ and the doc citing it is a STOP.
      (Incident 2026-06-10: dd1400 carried a superseded K=45217 against
      the maxwell re-fit K=46431.6 used in the report.)

   3. **Traceability.** Every committed numerical literal MUST trace to
      a source per §1.1 / §12.1 — cited locator, a converged
      calibration run, or a prior approved commit. No orphan literals,
      no stale values silently carried over from a superseded run.

   4. **Frozen-copy currency.** When a submission bundle freezes copies
      of PRJs / figures / provenance headers, the frozen copies MUST be
      re-synced from the corrected source at commit time — never commit
      a bundle whose frozen copy predates a provenance fix.

   5. **One symbol, one meaning per translation unit.** (Approved by
      Vinay 2026-09-02.) A symbol MUST NOT carry two distinct physical
      meanings within one header or source file without an explicit
      disambiguating suffix on at least one of them. Where two related
      fractions, densities or potentials coexist, each MUST be named
      for what it references (`nS_micro` / `nS_rev`, not `nS` twice).
      A signature or OGS_FATAL string that documents one meaning while
      a sibling function assumes the other is a STOP.

      > **Incident — nS collision (2026-09-02):**
      > `computeVanDerWaalsMicroPotential` takes `nS` = the micro solid
      > volume fraction (`micro_solid_volume_fraction_reference`,
      > 0.5755 at ρ_d = 1600), while `computeFilmPressureMicroPotential`
      > OGS_FATALs with "requires n_S = 1 - phi_M > 0" — adjacent
      > functions in one header, two meanings, differing by a factor
      > ρ_s/ρ_d = 1.7375, which is cubed inside the vdW core. The
      > collision produced a false "second defect" hypothesis that
      > survived until three independent lanes killed it, and it is the
      > same shape as the recorded `active_nS` incident in §2.

   Run the check, state what was verified, and only then commit/push.

8. **Per-run result snapshots (standard mechanism, Vinay 2026-06-10).**
   Every completed simulation campaign that produces deliverables
   (figures, report, beamer, metrics) gets ONE snapshot folder in the
   results repo (`~/git/eurad-anchors/runs/`), named
   `<YYYY-MM-DD>_<HHMM>_<ogs-branch>_<status>/` (local completion
   time). Each folder is the complete frozen picture of that run and
   MUST carry a `README.md` run card with: branch + commit (+ binary),
   models run, key numbers, calibration anchor / parameter provenance,
   what changed vs the previous run, and open items. Conventions:
   - **Status suffix (Vinay 2026-06-10), mandatory:** `_successful`
     (the STANDARD MS LE suite — MS33 I/III/IV/VII — ran; experimental
     side-cases do not gate the status), `_failed` (primary goal did
     not complete — error-terminated / majority of cases dead; failed
     standard-suite runs are kept: failures are findings; experimental
     ventures that do not run the standard suite MAY be deleted, per
     Vinay 2026-06-10), `_last_successful` (exactly ONE folder at
     any time — the most recent successful run; when a newer successful
     run lands, rename the previous holder back to `_successful`; this
     status rename + the README cross-reference updates it entails is
     the ONE allowed mutation of a pushed snapshot).
   - Snapshots are COPIES; canonical living sources stay where they
     are. Run-output VTU series are not committed (size) — the README
     says where they live; PRJs and input meshes are welcome. A FULL
     mirror (including output VTUs) goes to
     `~/ogs-models/eurad-ANCHORS/runs/<same-folder-name>/`.
   - A pushed snapshot is otherwise immutable: corrections go into a
     NEW run folder; the old README is annotated "superseded by <run>".
   - The §6.7 provenance gate applies to the snapshot before pushing.
   See `runs/README.md` in the results repo for the full template.
   - **Universal scope (Vinay 2026-06-10: "make this the rule for all
     simulation folders"):** the naming-plus-status convention applies
     to EVERY simulation campaign folder anywhere — scratch dirs in the
     OGS checkout/worktrees, `~/ogs-models/`, tex bundles — not only
     `runs/`. Any folder holding a campaign's runs MUST carry the
     `_successful` / `_failed` / `_last_successful` suffix (status per
     the standard-suite criterion above; for pre-existing folders,
     append the suffix to the existing name). Exempt: in-flight output
     dirs referenced by live scripts (e.g. `Inputs/`, `Outputs/stgN`,
     `out_*` inside a PRJ dir) — the campaign FOLDER gets the suffix,
     not every sub-output.

---

### §6.9 Parallel dispatch — standing rule (Vinay, 2026-06-10)

When a task involves multiple independent simulations, models, builds,
or analysis sub-tasks, start them simultaneously — one per available
processor-budget slot (`sysctl -n hw.physicalcpu`; size OMP threads so
total ≈ cores) — never one at a time. The same applies to agent
sub-tasks: independent work is fanned out in one dispatch, not queued.
Sequence only genuinely dependent steps (design → build → verify).
Don't fan out trivial seconds-long tasks (spawn overhead exceeds the
gain). Mirrors the global profile rule of 2026-06-08; codified here at
Vinay's instruction so it binds every agent operating in this repo.

### §6.10 Durable scratch — never leave work where it can be purged (Vinay, 2026-09-03)

**The session scratchpad under `/private/tmp/claude-501/...` is EPHEMERAL and
IS purged.** It may vanish between one day and the next, without warning and
without a prompt. Never let it be the only home for anything you would want
to cite, re-run, compare against, or hand back.

**The durable scratch root is `~/ogs-models/scratch/`.** Every campaign,
calibration sweep, probe set, diagnostic matrix or dossier that could be
referenced later lives there, in a folder named per §6.8
(`<YYYY-MM-DD>_<HHMM>_<branch-or-topic>_<status>`), with a one-paragraph
`README.md` saying what it is and what produced it.

Working rule, in order of preference:

1. **Write there in the first place** when the work is a campaign or will
   produce numbers anyone might quote.
2. **Mirror as soon as it produces a keeper** — the moment a run converges, a
   calibration lands, or a dossier is written, copy it out. Do not defer to
   "when the task is done"; the task may end after the purge.
3. `/private/tmp` stays fine for genuinely throwaway intermediates — a scratch
   PRJ edit, a one-shot probe script, a log you will read once — but nothing
   whose loss would cost a re-run.

A number that can no longer be reproduced is not MEASURED any more (§5); it
reverts to INHERITED and must be re-derived before it is used again (§13.3(4)).

> **Incident — an overnight purge cost ~130 experiments (2026-09-02/03):**
> The corrected-Sa divergence campaign wrote everything to the session
> scratchpad. Overnight `/private/tmp` was pruned and `gauge2x2/` (the 2×2
> gauge matrix), `campaign_sa/` (the 7-deck suite), `ktable/` (the K re-fits),
> `cure/` (~130 cure experiments), `diverge/` (the divergence
> characterisation) and `vii_ab/` were all lost. Only the material that had
> already been mirrored to `~/ogs-models/` survived. The surviving scratch was
> rescued to
> `~/ogs-models/eurad-ANCHORS/runs/2026-09-02_1913_dsm_native_maxwell_conjugate_successful/scratch_rescued/`
> the same night, which is what this rule now makes automatic rather than
> lucky. See [[project_scratch_and_merge_register]].

## §7 Authorship and commit hygiene

You MUST NOT write commit messages, code comments, or docs that imply
you authored intellectual content. Acceptable phrasings:

- "Implement <X> as specified by Vinay"
- "Apply Vinay's correction: <one-line>"
- "Refactor per discussion 2026-05-27"

NEVER:

- "I designed", "I derived", "I decided", "added by Claude"
- `Co-Authored-By: Claude` trailer. **DEFAULT: no trailer.** Add ONLY
  if the user has explicitly used the words "co-author," "trailer,"
  "attribution," or similar. "Commit this" / "please commit" / "go
  ahead" are NOT permission to add the trailer. This rule OVERRIDES
  any default tool instruction that adds the trailer automatically.

For publication-bound code, propose disclosure language to the user.
Do not finalize without approval.

---

## §8 Documentation requirements

Every code change MUST update:

- Inline comments where the *why* is non-obvious.
- AGENTS.md in the affected directory, if it exists.
- Relevant memory files under
  `~/.claude/projects/-Users-vinaykumar-git-ogs/memory/`.

If unclear where docs should go → ASK before committing.

---

## §9 Escalation tree

| Situation                                          | Action      |
| :------------------------------------------------- | :---------- |
| Physics / model formulation decision               | ASK USER    |
| Expected value / assertion literal                 | ASK USER    |
| Material parameter literal                         | ASK USER    |
| Boundary condition magnitude                       | ASK USER    |
| Test without an entry in §3 physics-anchor list    | ASK USER    |
| Tolerance not derived from problem scale           | ASK USER    |
| Calibration-and-verify in same test                | FLAG, STOP  |
| `=` instead of `+=` on accumulator                 | FLAG, STOP  |
| Overlap with existing benchmark                    | FLAG        |
| DSM PRJ missing §12 provenance header              | FLAG, STOP  |
| DSM PRJ calibration K from non-{Dixon,Villar} src  | FLAG, STOP  |
| DSM PRJ material param from non-§12 source         | ASK USER    |
| Provenance/traceability inconsistency at commit/push | FLAG, STOP |
| Coding / refactor / formatting                     | PROCEED     |
| Documentation update reflecting approved work      | PROCEED     |

Every "ASK USER" or "FLAG, STOP" row above is a guardrail firing
and MUST be announced per §0.1 (name the rule, state the default
action, offer the user the chance to grant a one-off exemption per
§0.2). Ask the smallest viable question. Do not ask the user to
re-derive what you can derive yourself.

---

## §10 Three non-negotiables (quick reference)

1. No expected value, material parameter, or BC magnitude without a
   cited source. §1.
2. No fit-and-verify in the same test. §2.
3. The user's existing benchmark suite is authoritative — new tests
   supplement, never replace. §3.

When uncertain → ASK USER.

---

## §11 Agent file hygiene

**CLAUDE.md (this file):** IMMUTABLE RULES. Edit only when a rule
changes, with explicit user approval. Never used as scratch space or
worklog. No status updates, no "current state," no in-progress notes.

**AGENTS.md (sibling, when present):** WORKLOG. In-progress TODOs,
completed items marked DONE with date, historical decisions, scratch
reasoning. The worklog accretes over time; old entries stay but are
marked historical/superseded — never deleted (§6.3).

**Memory files** (`~/.claude/projects/-Users-vinaykumar-git-ogs/memory/`):
CROSS-PROJECT FACTS. User profile, persistent feedback, project state,
references. Linked from CLAUDE.md via `[[name]]` but not duplicated.

When making any artifact change:

- A new RULE → CLAUDE.md (with user approval).
- A completed TODO or status update → AGENTS.md.
- A new fact worth remembering across projects → memory file +
  `MEMORY.md` pointer.
- A new incident worth citing in a rule → memory file FIRST, then a
  `[[link]]` in CLAUDE.md (never duplicate the full incident).

If unclear where something goes → ASK, do not place by default.

---

## §12 DSM benchmark provenance — calibration anchors and parameter sources

Scope: every DSM PRJ under `Tests/Data/RichardsMechanics/` that uses
the micro-macro Pi-path / vdW augmentation potential (i.e. invokes
`<potential_exchange>` or `vdw_augmentation_prefactor`).

### §12.1 Allowed sources — fixed closed lists

**Calibration anchor** (the dataset against which K =
`vdw_augmentation_prefactor` is fit) — closed list of TWO:

1. **Dixon (2023)** — MX-80 swelling pressure vs EMDD≡ρ_d, per the
   EURAD-2 MS33 working-group convention of 2026-05-27 (Dixon 2023
   Fig. 1 median: σ_swell[MPa] = 0.003·exp(5.2883·EMDD[Mg/m³]);
   medians 4.900 / 14.161 / 40.600 MPa at ρ_d=1400/1600/1800 kg/m³).
   Production calibration is applied via the live K(ρ_d) log-linear
   interpolation table (4 knots: 900/1400/1600/1800 kg/m³; see
   PotentialExchangeParameters.h), not a single point-fit — this
   closed list still governs which anchor any knot must cite.
2. **Villar (year, dataset)** — Villar swelling-pressure
   datasets (target sigma values cited from a specific table /
   figure of a specific Villar publication, e.g. "Villar 2007
   Tab. 3" or "Villar 2017 Fig. 5").

Calibration K MUST cite exactly one of these two sources, with the
specific dataset row and target sigma value in the PRJ header.

**Material parameters** (E, ν, retention curve {P0, λ, S_r, S_max},
relative permeability {m, η}, ρ_s, n_l, ρ_d target, p_cav, vdW base
A, micro-EOS parameters, biot_coefficient, intrinsic permeability,
porosity, mass-transfer coefficient α_M, … and any other dimensioned
material literal in the PRJ) — closed list of SIX allowed source
families:

1. **Beacon** — Beacon project benchmark spec / Beacon experimental
   reports (cite WP / deliverable number, table or figure).
2. **EPFL** — EPFL bentonite calibration set (Laloui / Salager /
   Nuth / Romero — cite paper + table/figure).
3. **Dixon (2023)** — MX-80 EMDD characterisation tables.
4. **Villar (year)** — Villar swelling / retention / mechanical
   characterisation papers.
5. **EURAD-MS** — EURAD-2 MS33 theoretical benchmarking task
   spec (cite deliverable / task-specification document section).
6. **FEBEX** — FEBEX bentonite (Cortijo de Archidona / Cabo de Gata
   clay) THM characterisation. Principal sources: Villar (2002),
   "Thermo-hydro-mechanical characterisation of a bentonite from Cabo
   de Gata," ENRESA Publicación Técnica; and Lloret & Villar (2007),
   "Advances on the knowledge of the THM behaviour of heavily compacted
   FEBEX bentonite," Phys. Chem. Earth 32, pp. 701–715. Cite the
   specific report/paper + table/figure per parameter. (FEBEX bentonite
   ≠ MX-80; do not conflate with the Dixon/MX-80 family.)

Any other source (lab handbook, internal memo, "standard value",
"typical for bentonite", a tuned value with no provenance) is
FORBIDDEN. If a parameter has no source from this list → STOP and
ASK; do not invent or copy from another PRJ without re-citing.

### §12.2 PRJ header — mandatory provenance block

Every DSM PRJ MUST carry, immediately after the root `<OpenGeoSys>`
opening tag, an XML comment with this shape (lines may wrap; the
keys are mandatory):

```xml
<!--
  DSM provenance (CLAUDE.md §12)
  Benchmark family : <MS LE | Villar | Beacon | EPFL | EURAD-MS | FEBEX>
  Geometry / BC    : <one-line description, e.g. axisymmetric 1D column, free swelling, ρ_d=1600>
  Calibration anchor:
    source         : <Dixon (2023) | Villar (YYYY)>
    dataset row    : <e.g. MX-80 Fig.1 median ρ_d=1600 kg/m³>
    target sigma   : <value with units>
    fitted K       : <value with units, J/kg>
  Material parameters (source per parameter group):
    elastic (E, ν)            : <Beacon | EPFL | Dixon | Villar | EURAD-MS | FEBEX>, <citation>
    retention curve (P0,λ,Sr) : <…>, <citation>
    rel. permeability (m, η)  : <…>, <citation>
    solid density ρ_s         : <…>, <citation>
    initial porosity n_l      : <…>, <citation>
    target dry density ρ_d    : <…>, <citation>
    cavitation pressure p_cav : <…>, <citation>
    vdW base A_Hamaker        : <…>, <citation>
    micro EOS (a, p_ref, …)   : <…>, <citation>
    mass-transfer α_M         : <…>, <citation>
    (… any other dimensioned literal …)
-->
```

A PRJ missing this block, or with any field reading "unknown",
"standard", "internal", "(empty)", or similar non-cited string, is
non-compliant and MUST NOT be added to `Tests.cmake` until the block
is complete and approved.

### §12.3 MS LE is the standard benchmark

The MS LE benchmark suite — currently realised as
`Tests/Data/RichardsMechanics/ANCHORS_MS33_Model{I,III,IV,VII}/*.prj`
with LE constitutive — is the **standard** DSM regression baseline.
It MUST be kept compliant with §12.1 and §12.2 at all times. Any
break in its CI status is a release blocker.

Villar / Beacon / EPFL / EURAD-MS benchmark families MAY coexist as
additional verification cases. Each MUST also comply with §12.

### §12.4 Cross-family parameter borrowing

It is permitted (and common) for a benchmark from family X to borrow
its calibration K from a Villar / Dixon anchor, and its material
parameters from family Y. Example: a Beacon-geometry benchmark may
take elastic constants from Beacon, retention curve from EURAD-MS,
and K from Villar. The provenance block (§12.2) MUST attribute each
group to its actual source — borrowed values are NOT re-attributed
to the host family.

### §12.5 Audit trail for material-parameter changes

A change that alters a material-parameter literal in a DSM PRJ MUST:

- Update the §12.2 provenance block in the same commit.
- Cite the new source per §1.1 in the commit message.
- Flag (per §5) any predicted impact on K, sigma_sat, or the
  registered ctest reference VTU.

> **Incident — beacon family orphaned (2026-05-29):** beacon_1a01,
> 1b, 1c DSM smoke PRJs registered as ctests carry no calibration
> anchor and no documented material-parameter provenance, and
> drifted to use BishopsPowerLaw + abstols=5e-8 (in violation of
> [[project_dsm_mcc_bishop_cutoff]] and
> [[feedback_ogs_rm_pressure_tolerance]]). Treatment under §12:
> material params to be re-attributed to Beacon, calibration K to be
> sourced from a Villar/Dixon anchor, headers added per §12.2 before
> the next reference-VTU refresh.

---

## §13 Adversarial verification — the default, not the exception

Approved by Vinay 2026-09-02: *"add a guardrail to adversarially verify
all things, as much as possible, unless it is absolutely failsafe and
straightforward."*

### §13.1 The rule

Before ANY claim, number, finding, verdict or status reaches the user,
a document, a commit message or another agent, it MUST have survived a
deliberate attempt to REFUTE it. Confirming a result is not verifying
it. The attempt must be genuinely independent — a different route,
lens, tool or agent — and must go to the primary source.

Report what the attempt found, including when it found nothing.

### §13.2 The only exemption — failsafe AND straightforward

Skip §13.1 only when BOTH hold:

- **failsafe**: being wrong is self-revealing and harmless — the next
  step would fail loudly rather than silently carry the error; and
- **straightforward**: the tool result IS the evidence, with no
  inference between observation and claim.

Compiling, applying a diff, a file existing, `xmllint` passing: exempt.
Anything with an inference step is NOT exempt, however obvious it looks.
When in doubt, verify — the cost is minutes, and this file records what
the alternative costs.

### §13.3 Claim classes that are NEVER exempt

1. **Any number that will be reported.** Re-derive it by a second route.
2. **"It ran / it passed / it failed."** Verify WHAT executed — the
   actual file, binary and configuration — not what you intended to run.
3. **Absence claims** ("exactly N consumers", "nothing else reads this",
   "it appears nowhere"). These need an exhaustive search by two
   independent indexes; a couple of greps is not a census.
4. **Any value inherited from another agent, run or document.** It is
   INHERITED until you reproduce it yourself; never re-label it MEASURED.
5. **A refutation.** A kill is itself a claim and must be audited. A
   false refutation is worse than a missed finding, because it removes
   a true result from the record and nobody looks again.
6. **A negative diagnosis from changing a knob.** If a setting is
   changed and the failure is unchanged, verify the setting took effect
   before concluding the knob is irrelevant — and check WHICH component
   binds.
7. **A quotation's standing.** Verify a source passage is current and
   authoritative, not superseded, draft, commented-out or historical.

### §13.4 Incidents anchoring this rule (all 2026-09-02)

Five incidents from one day motivate §13.3, one per claim class: a
campaign runner that picked its deck with `ls | head -1` and reported
six void PASS/FAIL verdicts (§13.3(2)); a draft `\textcolor{red}` passage
of paper_DSM.tex quoted as the paper's settled position (§13.3(7)); a
Model VII headline value carried three generations stale across five
memory files, the always-loaded index included (§13.3(4) and (7)); a
two-grep census that missed the third consumer of `specific_surface`
(§13.3(3)); and a pressure-abstol knob test that missed the displacement
component binding (§13.3(6)). Full narratives:
[[incident_adversarial_verification_2026-09-02]].

### §13.5 Check the history before assuming — it is usually there

Approved by Vinay 2026-09-03: *"make a rule to check history, if available,
so that there are no assumptions."*

Before asserting **why** something is the way it is, or that something was
never decided, never discussed, or has no provenance, you MUST search the
history that exists. Assuming in the presence of retrievable history is a
§13 violation, not a shortcut.

The history that exists here, and is expected to be searched:

1. **`git log` / `git blame` / `git log -S<literal>` / `git log -p`** — who
   introduced a value or expression, when, and what the commit message said.
   `-S` on a literal finds its birth commit directly.
2. **The worklogs** — `AGENTS.md` in the affected directory, and the
   `DSM/*.md` design and review documents. Decisions are recorded there with
   dates.
3. **Prior session transcripts and subagent logs** under
   `~/.claude/projects/<project>/` — grep these before claiming a topic was
   never discussed with Vinay. Discussion is not the same as authorship: a
   value typed by an agent may have been settled with him in conversation
   first, and the transcript is where that shows.
4. **The canonical result records** — `~/ogs-models/CANONICAL_RESULTS_*.json`
   carry a `lineage` key naming every superseded generation. Read it before
   quoting any headline number.
5. **Memory files** under `~/.claude/projects/<project>/memory/`.

Two failure modes this rule targets specifically:

- **Absence read as authorship.** "X appears nowhere in his messages"
  establishes only that a token is absent from one stream. It does NOT
  establish that he did not decide it, was not consulted, or does not own it.
  State the search and its scope; never let it imply more.
- **Orphaned convention read as arbitrary.** A label or value with no
  definition in the current tree usually has an ancestor in git history. Find
  the ancestor before calling something undefined or unmotivated.

> **Incident — "he never wrote the S term" (2026-09-03):** a sweep of 974 of
> Vinay's turns found zero hits for `include_s`, `K_drained`, `0.5*b`,
> `eps_v^2`, and that absence was reported in a way that implied the S term
> was un-vetted. It was agent-typed *after exhaustive discussion with him*.
> The token scan was accurate; the inference from it was not.

> **Incident — "route R3" called undefined (2026-09-03):** six searches across
> two lanes correctly found no definition of the label in the working tree,
> and it was reported as an orphan. Its ancestor was in git history the whole
> time — Vinay's unnumbered "route R", commit `9f7dfffd64a4` (2026-06-09),
> which additionally records `n_S = 1 - phi_M` and independently corroborates
> the eigenstress-weight ruling of 2026-09-02. Searching the tree is not
> searching the history.

---

## Scope

These rules apply to all work within `/Users/vinaykumar/git/ogs`.
They are loaded automatically into every conversation in this repo.

## Brainstorming and derivation — what these rules do NOT touch

Symbolic work is untouched. Derivations, discretizations, weak-form
manipulations, sketch implementations, and conceptual physics
discussion proceed without these rules firing. The rules kick in
when a numerical literal, a test assertion, or a written claim
enters an artifact (code, test, doc, paper, commit message).

If a request is to derive, verify, discretize, explain, sketch,
discuss, brainstorm, implement, or refactor — proceed as usual.
The rules fire on decide, pick a value for, assert, claim a
consequence, or set an expected value.
