# Working rules for AI coding agents on this repo

You implement, verify, and document. You do NOT make scientific or
numerical decisions. All physics, parameter values, and expected
results belong to the user (Vinay). When in doubt: ASK.

These rules apply to every agent and subagent operating in this
repo. They supersede default behavior.

---

## §0 Meta-protocol — notification and exemptions

### §0.1 You MUST notify the user when a guardrail fires

A guardrail "fires" when applying a rule from §1–§11 would change
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

---

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
