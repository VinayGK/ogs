> PROVENANCE OF THIS DOCUMENT. Produced 2026-08-12 by an 11-agent mining run
> (5 tracers, one per deck -> 5 adversarial verifiers -> 1 synthesis), on
> branch `dsm_native_maxwell_conjugate` @ `40551b6aad`. Companion to
> `DSM_BENCHMARK_DEVIATION_AUDIT_2026-06-06.md`, which it extends for the
> beacon family. Read-only: no file was edited to produce it.
>
> CAVEAT ON SECTION 0. The run overlapped a `git stash` of the beacon working
> set (the restore + inert-floor edit was parked while a clean-tree ctest
> baseline ran). Section 0's "absent from branch tip" finding is an artifact of
> that timing, not a repository defect; the decks were committed to
> `dsm_native_maxwell_conjugate` the same day. Everything else in the document
> stands, and all line numbers refer to the version that was committed.

# BEACON DSM PRJ — Consolidated Provenance Ledger

**Scope:** 5 decks — `beacon_1a01_dsm_micromacro_smoke.prj`, `beacon_1a01_dsm_micromacro_stressprobe.prj`, `beacon_1a01_dsm_micromacro_inflow.prj`, `beacon_1b_dsm_micromacro_smoke.prj`, `beacon_1c_dsm_micromacro_smoke.prj`.
**Rule applied:** a literal carries either a citation whose verification verdict was CONFIRMED, or the literal string `TODO(Vinay): UNSOURCED`. Every `probable`, `REFUTED`, and `UNVERIFIABLE` attribution has been collapsed to `TODO(Vinay): UNSOURCED`.

---

## 0. BLOCKING PROCESS FINDING — read before transcribing anything

All five files are **absent from branch tip `40551b6aad`** and from the working tree of `/Users/vinaykumar/git/ogs-worktrees/dsm_native_maxwell_conjugate_wt`. `git cat-file -e HEAD:<path>` fails; `git ls-files Tests/Data/RichardsMechanics | grep -i beacon` is empty. `Tests/Data` is a plain tree (no submodule) and the paths are not gitignored.

They resolve only in:
- **`stash@{0}` = commit `131a2c9599`** ("On dsm_native_maxwell_conjugate: beacon_wip_2026-08-12"). This is the only copy carrying the 2026-08-12 inert-floor pair. `git branch -a --contains 131a2c9599` returns nothing.
- Committed copies on **deprecated branches only** (`deprecated/dsm_native_hierarchical`, `deprecated/dsm_native`, `deprecated/dsm_native_tuller_macro_film`, `deprecated/dsm_native_tuller_review`, `dsm-nb-transition`, `archive/claude_dazzling-rubin-9d242d`; commit `72f4f3a192`). **No committed version contains the floor keys**, and the canonical binary requires `macro_porosity_floor`/`micro_water_content_floor` as a pair (parser change `71366ac0d3`) — a git-restored copy will not parse.
- On-disk audit pins `/Users/vinaykumar/git/ogs-worktrees/_dsm_audit_2026-06-06/{hierarchical,maxwell_conjugate}/…` (md5 `25c216537ed9c75e24e92edf8e1d7ec1` for 1a01 smoke) — a **different** version from the stash (`2d5e02c10024c3c628c05af55eb4bd91`).
- A **divergent third copy** at `/Users/vinaykumar/git/GitHub/ogs` (branch `dsm-nb-transition`): `biot_coefficient=0.6` not 1.0, no floors, and it carries `<mode>full_potential</mode>` which the stash blob lacks.

**Consequences.** (a) `ProcessLib/RichardsMechanics/Tests.cmake` still registers all five (lines 76–121, 85–93, 96 + 125, 147–173, 176–203) — registered ctests point at inputs absent from the branch, and `beacon_1a01.gml` / `beacon_1a01_domain.vtu` / `beacon_1c_domain.vtu` are missing too. (b) Any §12.2 header written from this ledger documents **uncommitted stashed WIP** and must name the artifact (stash `131a2c9599`) it describes, or the deck must be committed first. (c) A §6.7(2) cross-artifact inconsistency already exists between the stash copy and the `dsm-nb-transition` copy.

All line numbers below are the **stash `131a2c9599`** version. The 2026-06-06 audit's line numbers are ~7 lower (the floor comment block shifted them).

---

## 1. Calibration anchor — plainly stated

**None of the five decks has a calibration anchor.**

- The tag named in the brief, `vdw_augmentation_prefactor`, **does not exist on this branch** — grep over the whole source tree returns zero hits. The current tags are `potential_augmentation_prefactor` (`CreateRichardsMechanicsProcess.cpp:555-557`) and the table variant `potential_augmentation_prefactor_vs_dry_density` (`:502`).
- **Neither tag appears in any of the five decks.** K therefore defaults to `0.0`; `PotentialExchange.h:259` gates the augmentation on `K > 0` and skips the branch; `PotentialExchange.h:118` states "Setting K = 0 (default) reduces exactly to the pure vdW form."
- **No deck asserts a σ_sat / swelling-pressure target.** The `swelling_pressures` literals (1.3e7 / 1.1e7 / 0.8e7 Pa) are *uncited inputs* to the macro `SaturationDependentSwelling` law — they are **not** calibration targets and must not be presented as a Dixon/Villar anchor row.
- In the three **smoke** decks the vdW base is *also* switched off (§5), so the entire Π-path potential is inactive: they exercise micro–macro exchange *plumbing* only. In **stressprobe** (A=6e-20, Sa=100) and **inflow** (A=5.1e-21, Sa=523) the vdW base is live, but K is still absent.

**Mandatory §12.2 "Calibration anchor" block for all five decks — transcribe verbatim:**

```
Calibration anchor:
  source          : none — no calibration performed
  dataset row     : n/a
  target sigma    : n/a
  fitted K        : none — no <potential_augmentation_prefactor> tag;
                    K defaults to 0.0 and the augmentation branch is
                    skipped (PotentialExchange.h:259). This deck is a
                    micro-macro exchange PLUMBING SMOKE TEST, not a
                    calibrated benchmark. §12.1's Dixon/Villar rule is
                    vacuous here, not violated.
```

---

## 2. The hamaker_constant = 1e-30 / specific_surface = 1.0 finding

**Verdict: deliberate vdW switch-off placeholder, mechanism CONFIRMED from source; author intent UNDOCUMENTED.**

Applies to the **three smoke decks only**: `beacon_1a01_…smoke.prj:21-22`, `beacon_1b_…smoke.prj:17-18`, `beacon_1c_…smoke.prj:17-18` (1c per-medium overrides Sa = 0.8 block `:30`, 1.5 pellet `:38`). It does **not** apply to stressprobe (6e-20 / 100) or inflow (5.1e-21 / 523).

Evidence:

1. **Both parameters enter only through one multiplicative prefactor.** `ProcessLib/RichardsMechanics/ConstitutiveRelations/PotentialExchange.h:218-220`: `prefactor = A · Sa³ / (6π)`; `:232-234`: `mu_lR_vdW = prefactor · nS³ · ρ_SR³ / (n_l³ · ρ_lR)` [J/kg]. Nothing else in the model reads A or Sa. Driving the prefactor to zero removes the entire microscale vdW potential and, via `p_L_m = −ρ_lR·mu_lR`, the micro pressure with it.
2. **Zero is forbidden by the code, so a tiny positive value is the only way to disable it.** `CreateRichardsMechanicsProcess.cpp:304-316` (`get_positive_required_or_default`) `OGS_FATAL`s unless the value is `> 0`, applied at `:362-364` (hamaker) and `:365-367` (specific_surface); `PotentialExchange.h:158-171` repeats the guard at evaluation time. An author wanting A = 0 **cannot write it**. `A=1e-30, Sa=1.0` is the minimal positivity-satisfying stand-in.
3. **Magnitude — numerically dead, not merely small.** Against the physical MS33 pairing (A=2.2e-20 J, Sa=523): `A·Sa³` = 1e-30 vs 3.15e-12, ratio **3.2e-19** taking Sa as entered; **~3e-28** if Sa is first converted to SI (5.23e5 m²/kg). Evaluated on the smoke decks' own state (nS=0.6, ρ_SR=2650, n_l=0.1, ρ_lR≈1e3) `mu_lR_vdW ≈ 2.1e-22 J/kg`, i.e. `p_L_m ≈ −2e-19 Pa`, against a macro potential `mu_LR ≈ p/ρ ≈ 1e3 J/kg` and a 1e6 Pa suction scale — **~25 orders down**. (Arithmetic from file values; **predicted, not verified by a run** — §5.)
4. **The value pair is the DSM unit-test harness state.** The identical quintuple (A=1e-30, Sa=1.0, `micro_solid_density_reference`=2650, `micro_solid_volume_fraction_reference`=0.6, `initial_micro_water_content`=0.1) appears verbatim at `Tests/ProcessLib/RichardsMechanics/DSMMicroMacroSingleIntegrationPoint.cpp:663-667`, `:753-757`, `:1793-1797` — tests of the exchange/tangent machinery, not of vdW physics. The smoke decks' `potential_exchange` block is that harness state transplanted into a PRJ.
5. **Sa = 1.0 m²/kg is physically absurd for bentonite** (MX-80 is O(1e5–1e6) m²/kg), confirming it is a unit-valued dummy paired with A, not a measured area.
6. **The standing audit reads it the same way:** `DSM_BENCHMARK_DEVIATION_AUDIT_2026-06-06.md:197,199` — "hamaker_constant=1e-30 (L17, not the literature 2.2e-20 J Israelachvili & Adams 1978 — smoke placeholder)".

**Limit on the evidence, which must not be papered over:** no commit message and no in-file comment states the intent. `git log -S'1e-30'` on these paths surfaces only the 2026-04-10 rename commits and the stash. The standing audit still treats both literals as *material* literals requiring either a §12.1 family or a §0.2 exemption.

**Header treatment:**
```
vdW base A_Hamaker  : NUMERICAL PLACEHOLDER, term switched OFF.
    A=1e-30 J with Sa=1.0 m^2/kg drives the vdW prefactor
    A*Sa^3/(6*pi) to ~1e-19 of its physical value; exact zero is
    rejected by the parser (PotentialExchange.h:158-171), so a
    minimal positive pair is the only way to disable the branch.
    NOT a material parameter; NO §12.1 citation applies.
    The literature constant it stands in for is A = 2.2e-20 J
    (Israelachvili & Adams 1978) — DO NOT write that attribution
    here. Requires a §0.2 GUARDRAIL EXEMPTION note from Vinay
    recording the switch-off intent; until granted, the fallback
    classification is TODO(Vinay): UNSOURCED.
```

---

## 3. Per-PRJ ledgers

### 3.1 `beacon_1a01_dsm_micromacro_smoke.prj`

**MATERIAL**

| tag | value | unit | attribution |
|---|---|---|---|
| `micro_porosity/mass_exchange_coefficient` | 1e-13 | undocumented in OGS (see §6, unit item) | TODO(Vinay): UNSOURCED |
| `potential_exchange/micro_solid_density_reference` | 2650 | kg/m³ | TODO(Vinay): UNSOURCED |
| `potential_exchange/micro_solid_volume_fraction_reference` | 0.6 | dimensionless | TODO(Vinay): UNSOURCED |
| `potential_exchange/initial_micro_water_content` | 0.1 | dimensionless (n_l0) | TODO(Vinay): UNSOURCED |
| `AqueousLiquid/viscosity` | 1e-3 | Pa·s | TODO(Vinay): UNSOURCED |
| `AqueousLiquid/density` | 1e3 | kg/m³ | TODO(Vinay): UNSOURCED |
| `Solid/density` | 2780 | kg/m³ | TODO(Vinay): UNSOURCED |
| `swelling_stress_rate/swelling_pressures` | 1.3e7 ×3 | Pa | TODO(Vinay): UNSOURCED |
| `swelling_stress_rate/exponents` | 1 1 1 | dimensionless | TODO(Vinay): UNSOURCED |
| `swelling_stress_rate/lower_saturation_limit` | 0 | dimensionless | TODO(Vinay): UNSOURCED |
| `swelling_stress_rate/upper_saturation_limit` | 1 | dimensionless | TODO(Vinay): UNSOURCED |
| `biot_coefficient` | 1.0 | dimensionless | TODO(Vinay): UNSOURCED |
| `permeability` | 1e-18 | m² | TODO(Vinay): UNSOURCED |
| `reference_temperature` | 293.15 | K | TODO(Vinay): UNSOURCED |
| `relative_permeability` (Constant) | 1 | dimensionless | TODO(Vinay): UNSOURCED |
| `saturation/residual_liquid_saturation` | 0.0 | dimensionless | TODO(Vinay): UNSOURCED |
| `saturation/residual_gas_saturation` | 0 | dimensionless | TODO(Vinay): UNSOURCED |
| `saturation/exponent` (van Genuchten m) | 0.5 | dimensionless | TODO(Vinay): UNSOURCED |
| `saturation/p_b` | 1e6 | Pa | TODO(Vinay): UNSOURCED |
| `saturation_micro/residual_liquid_saturation` | 0 | dimensionless | TODO(Vinay): UNSOURCED |
| `saturation_micro/residual_gas_saturation` | 0 | dimensionless | TODO(Vinay): UNSOURCED |
| `saturation_micro/exponent` | 0.5 | dimensionless | TODO(Vinay): UNSOURCED |
| `saturation_micro/p_b` | 1e7 | Pa | TODO(Vinay): UNSOURCED |
| `bishops_effective_stress/exponent` (BishopsPowerLaw) | 1 | dimensionless | TODO(Vinay): UNSOURCED |
| `parameter[E]` | 50e6 | Pa | TODO(Vinay): UNSOURCED |
| `parameter[nu]` | 0.2 | dimensionless | TODO(Vinay): UNSOURCED |
| `parameter[phi0]` | 0.4 | dimensionless | TODO(Vinay): UNSOURCED |
| `parameter[phi_tr0]` | 0.3 | dimensionless | TODO(Vinay): UNSOURCED |

**vdW pair (reclassified — see §2):** `hamaker_constant` 1e-30 J and `specific_surface` 1.0 m²/kg → NUMERICAL placeholder, needs a §0.2 exemption note, **not** a citation; if the placeholder reading is not ratified, both revert to `TODO(Vinay): UNSOURCED`.

**BC_IC**

| tag | value | unit | attribution |
|---|---|---|---|
| `specific_body_force` | 0 0 | m/s² | TODO(Vinay): UNSOURCED |
| `parameter[sigma0]` (4 components) | 0, 0, 0, 0 | Pa | TODO(Vinay): UNSOURCED |
| `parameter[displacement0]` | 0 0 | m | TODO(Vinay): UNSOURCED |
| `parameter[pressure_ic]` | −1e6 | Pa (1 MPa suction) | TODO(Vinay): UNSOURCED |
| `parameter[dirichlet0]` | 0 | m | TODO(Vinay): UNSOURCED |
| `parameter[top_pressure]` | 2e3 | Pa | TODO(Vinay): UNSOURCED |

Notes: §3 σ0 free-boundary rule does **not** fire — all four displacement boundaries carry a normal Dirichlet constraint (rollers, `:143-146`), so the cell is confined; the audit produced no bucket-H finding for this file. Live deviations to carry into the header: pressure `abstol = 5e-8 Pa` at a 1 MPa suction scale with no `<reltols>` (`:106`, audit bucket E); `BishopsPowerLaw` where the DSM standard is `BishopsSaturationCutoff(cutoff=1)` (`:95`, audit bucket F). Registered ctest with a **frozen reference VTU** (`beacon_1a01_reference_t_1000.000000.vtu`, tolerances to 1e-16) — giving the placeholders physical values will require a reference refresh (§12.5).

---

### 3.2 `beacon_1a01_dsm_micromacro_stressprobe.prj`

**MATERIAL**

| tag | value | unit | attribution |
|---|---|---|---|
| `mass_exchange_coefficient` | 1e-13 | undocumented | TODO(Vinay): UNSOURCED |
| `hamaker_constant` | 6e-20 | J | TODO(Vinay): UNSOURCED |
| `specific_surface` | 100 | m²/kg | TODO(Vinay): UNSOURCED |
| `micro_solid_density_reference` | 2650 | kg/m³ | TODO(Vinay): UNSOURCED |
| `micro_solid_volume_fraction_reference` | 0.6 | dimensionless | TODO(Vinay): UNSOURCED |
| `initial_micro_water_content` | 0.1 | dimensionless | TODO(Vinay): UNSOURCED |
| `micro_liquid_density_reference` | 1300 | kg/m³ | TODO(Vinay): UNSOURCED |
| `micro_liquid_density_a` | 1.3 | dimensionless | TODO(Vinay): UNSOURCED |
| `micro_liquid_density_b` | 1.0 | dimensionless | TODO(Vinay): UNSOURCED |
| `AqueousLiquid/viscosity` | 1e-3 | Pa·s | TODO(Vinay): UNSOURCED |
| `AqueousLiquid/density` | 1e3 | kg/m³ | TODO(Vinay): UNSOURCED |
| `Solid/density` | 2780 | kg/m³ | TODO(Vinay): UNSOURCED |
| `swelling_pressures` | 1.3e7 ×3 | Pa | TODO(Vinay): UNSOURCED |
| `swelling exponents` | 1 1 1 | dimensionless | TODO(Vinay): UNSOURCED |
| `lower/upper_saturation_limit` | 0 / 1 | dimensionless | TODO(Vinay): UNSOURCED |
| `biot_coefficient` | 1.0 | dimensionless | TODO(Vinay): UNSOURCED |
| `permeability` | 1e-18 | m² | TODO(Vinay): UNSOURCED |
| `reference_temperature` | 293.15 | K | TODO(Vinay): UNSOURCED |
| `relative_permeability` | 1 | dimensionless | TODO(Vinay): UNSOURCED |
| `saturation` S_r / S_gr / m / p_b | 0.0 / 0 / 0.5 / 1e6 | –, –, –, Pa | TODO(Vinay): UNSOURCED |
| `saturation_micro` S_r / S_gr / m / p_b | 0 / 0 / 0.5 / 1e7 | –, –, –, Pa | TODO(Vinay): UNSOURCED |
| `bishops_effective_stress/exponent` | 1 | dimensionless | TODO(Vinay): UNSOURCED |
| `parameter[E]` | 50e6 | Pa | TODO(Vinay): UNSOURCED |
| `parameter[nu]` | 0.2 | dimensionless | TODO(Vinay): UNSOURCED |
| `parameter[phi0]` | 0.4 | dimensionless | TODO(Vinay): UNSOURCED |
| `parameter[phi_tr0]` | 0.3 | dimensionless | TODO(Vinay): UNSOURCED |

**BC_IC**

| tag | value | unit | attribution |
|---|---|---|---|
| `specific_body_force` | 0 0 | m/s² | TODO(Vinay): UNSOURCED |
| `parameter[sigma0]` ×4 | 0, 0, 0, 0 | Pa | TODO(Vinay): UNSOURCED |
| `parameter[displacement0]` | 0 0 | m | TODO(Vinay): UNSOURCED |
| `parameter[pressure_ic]` | −1e6 | Pa | TODO(Vinay): UNSOURCED |
| `parameter[dirichlet0]` | 0 | m | TODO(Vinay): UNSOURCED |
| `parameter[top_pressure]` | 2e3 | Pa | TODO(Vinay): UNSOURCED |

Notes: **A = 6e-20 J lies outside the range the source file itself documents** — `PotentialExchange.h:97-101`: "A = 2.2e-20 J (Israelachvili & Adams 1978…)", "Range: 1–5e-20 J (smectite, DLVO literature)", "DO NOT calibrate A". The vdW term is **live** here (prefactor ≈ 3.2e-15 vs 5.3e-32 in the smoke decks) — the switch-off reading does **not** transfer. **§2 EOS warning:** `micro_liquid_density_a = 1.3`, not the `1e-16` bypass used by the MS33 form-(a) decks — the micro EOS is ACTIVE, so results and calibrations obtained under the bypass do not transfer either way; state this in the header. The only textual "source" for the 1300/1.3/1.0 triple is **circular**: `beacon_1c_…smoke.prj:22-24` says it was taken from the 1a01 case, which has no header. Registered ctest, **no DIFF_DATA** (run-only), so §2 fit-and-verify does not fire.

---

### 3.3 `beacon_1a01_dsm_micromacro_inflow.prj`

**MATERIAL**

| tag | value | unit | attribution |
|---|---|---|---|
| `mass_exchange_coefficient` | 1e-13 | undocumented | TODO(Vinay): UNSOURCED |
| `hamaker_constant` | 5.1e-21 | J | TODO(Vinay): UNSOURCED |
| **`specific_surface`** | **523** | **as entered, unitless numeral; source states m²/g; code consumes m²/kg — see §6** | **EPFL — Seiphoori, Ferrari & Laloui (2014), *Géotechnique* 64(9):721–734, Table 1, p.724, "Specific surface area, s: m²/g" = 523, MX-80 (Wyoming) granular bentonite, footnote "Data from Plötze & Weber (2007)". doi:10.1680/geot.14.P.017. Verified verbatim.** |
| `micro_solid_density_reference` | 2650 | kg/m³ | TODO(Vinay): UNSOURCED |
| `micro_solid_volume_fraction_reference` | 0.6 | dimensionless | TODO(Vinay): UNSOURCED |
| `initial_micro_water_content` | 0.01 | dimensionless | TODO(Vinay): UNSOURCED |
| `AqueousLiquid/viscosity` | 1e-3 | Pa·s | TODO(Vinay): UNSOURCED |
| `AqueousLiquid/density` | 1e3 | kg/m³ | TODO(Vinay): UNSOURCED |
| `Solid/density` | 2780 | kg/m³ | TODO(Vinay): UNSOURCED |
| `swelling_pressures` | 1.3e7 ×3 | Pa | TODO(Vinay): UNSOURCED |
| `swelling exponents` | 1 1 1 | dimensionless | TODO(Vinay): UNSOURCED |
| `lower/upper_saturation_limit` | 0 / 1 | dimensionless | TODO(Vinay): UNSOURCED |
| `biot_coefficient` | 1.0 | dimensionless | TODO(Vinay): UNSOURCED |
| `permeability` | 1e-18 | m² | TODO(Vinay): UNSOURCED |
| `reference_temperature` | 293.15 | K | TODO(Vinay): UNSOURCED |
| `relative_permeability` | 1 | dimensionless | TODO(Vinay): UNSOURCED |
| `saturation` S_r / S_gr / m / p_b | 0.0 / 0 / 0.5 / 1e6 | –, –, –, Pa | TODO(Vinay): UNSOURCED |
| `saturation_micro` S_r / S_gr / m / p_b | 0 / 0 / 0.5 / 1e7 | –, –, –, Pa | TODO(Vinay): UNSOURCED |
| `bishops_effective_stress/exponent` | 1 | dimensionless | TODO(Vinay): UNSOURCED |
| `parameter[E]` | 50e6 | Pa | TODO(Vinay): UNSOURCED |
| `parameter[nu]` | 0.2 | dimensionless | TODO(Vinay): UNSOURCED |
| `parameter[phi0]` | 0.4 | dimensionless | TODO(Vinay): UNSOURCED |
| `parameter[phi_tr0]` | 0.3 | dimensionless | TODO(Vinay): UNSOURCED |

No micro-EOS literals: `local_nonlinear_solve_mode = scalar_microstate_storage_mode` → `LocalNonlinearSolveMode::ScalarReferenceStorage` (`CreateRichardsMechanicsProcess.cpp:57-60`), which makes the micro-liquid-density EOS keys optional (`:373-406`). Header line: `micro EOS : not used (mode = scalar_microstate_storage_mode)`.

**BC_IC**

| tag | value | unit | attribution |
|---|---|---|---|
| `specific_body_force` | 0 0 | m/s² | TODO(Vinay): UNSOURCED |
| `parameter[sigma0]` ×4 | 0, 0, 0, 0 | Pa | TODO(Vinay): UNSOURCED |
| `parameter[displacement0]` | 0 0 | m | TODO(Vinay): UNSOURCED |
| `parameter[pressure_ic]` | +2e3 | Pa | TODO(Vinay): UNSOURCED |
| `parameter[dirichlet0]` | 0 | m | TODO(Vinay): UNSOURCED |
| `parameter[top_pressure]` | 2e3 | Pa | TODO(Vinay): UNSOURCED |

Notes:
- **`hamaker_constant = 5.1e-21` is a stale, superseded value in a registered ctest.** Commit `0d579e8aeb` (2026-05-22) states verbatim: "All 27 PRJ files: hamaker_constant 5.1e-21 -> 2.2e-20 J (Israelachvili & Adams 1978, mica-water-mica SFA, standard smectite proxy); A is a material constant — NOT calibrated to swelling pressure". **The beacon decks were missed by that sweep.** The only artifact ever describing 5.1e-21 is the header of the since-deleted sibling `beacon_1a01_dsm_micromacro_augmented_inflow.prj` at commit `8daad6f894`: "A_lit = 5.1e-21 J (literature value, NOT an effective multiplier)" — a bare "literature value" claim with no author/year/table, i.e. non-citing under §1.1. **Writing "Israelachvili & Adams 1978" next to 5.1e-21 would be a fabricated attribution** — that source produced 2.2e-20.
- The deck therefore pairs a *sourced* Sa (523, EPFL) with an *unsourced, superseded* A. Flag this inconsistency in the header.
- `initial_micro_water_content = 0.01` is declared a tuned knob by the file's own comment (`:18-22`): "make the attractive vdW microscale term larger by lowering the initial micro water content" — a modelling choice, no measured w₀ cited.
- `pressure_ic = +2e3` is likewise declared a numerical choice (`:18-20`, "keep the macro branch on the saturated helper side (mu_LR = 0) by using positive liquid pressure"); it equals `top_pressure`, so the hydraulic drive is zero at t=0.
- Registered **twice** in Tests.cmake (`:96` run, `:125` reference), the latter diffing six fields against `beacon_1a01_dsm_micromacro_inflow_reference_t_100000.000000.vtu`.

---

### 3.4 `beacon_1b_dsm_micromacro_smoke.prj`

**MATERIAL**

| tag | value | unit | attribution |
|---|---|---|---|
| `mass_exchange_coefficient` | 1e-13 | undocumented | TODO(Vinay): UNSOURCED |
| `micro_solid_density_reference` | 2650 | kg/m³ | TODO(Vinay): UNSOURCED |
| `micro_solid_volume_fraction_reference` | 0.6 | dimensionless | TODO(Vinay): UNSOURCED |
| `initial_micro_water_content` | 0.1 | dimensionless | TODO(Vinay): UNSOURCED |
| `AqueousLiquid/viscosity` | 1e-3 | Pa·s | TODO(Vinay): UNSOURCED |
| `AqueousLiquid/density` | 1e3 | kg/m³ | TODO(Vinay): UNSOURCED |
| `Solid/density` | 2780 | kg/m³ | TODO(Vinay): UNSOURCED |
| `swelling_pressures` | 1.1e7 ×3 | Pa | TODO(Vinay): UNSOURCED |
| `swelling exponents` | 1 1 1 | dimensionless | TODO(Vinay): UNSOURCED |
| `lower/upper_saturation_limit` | 0 / 1 | dimensionless | TODO(Vinay): UNSOURCED |
| `biot_coefficient` | 1.0 | dimensionless | TODO(Vinay): UNSOURCED |
| `permeability` | 5e-16 | m² | TODO(Vinay): UNSOURCED |
| `reference_temperature` | 293.15 | K | TODO(Vinay): UNSOURCED |
| `relative_permeability` | 1 | dimensionless | TODO(Vinay): UNSOURCED |
| `saturation` S_r / S_gr / m / p_b | 0.0 / 0 / 0.5 / 5e5 | –, –, –, Pa | TODO(Vinay): UNSOURCED |
| `saturation_micro` S_r / S_gr / m / p_b | 0 / 0 / 0.5 / 5e6 | –, –, –, Pa | TODO(Vinay): UNSOURCED |
| `bishops_effective_stress/exponent` | 1 | dimensionless | TODO(Vinay): UNSOURCED |
| `parameter[E]` | 40e6 | Pa | TODO(Vinay): UNSOURCED |
| `parameter[nu]` | 0.2 | dimensionless | TODO(Vinay): UNSOURCED |
| `parameter[phi0]` | 0.5 | dimensionless | TODO(Vinay): UNSOURCED |
| `parameter[phi_tr0]` | 0.4 | dimensionless | TODO(Vinay): UNSOURCED |

**vdW pair:** `hamaker_constant` 1e-30 J / `specific_surface` 1.0 m²/kg → NUMERICAL placeholder per §2 (§0.2 exemption required; otherwise TODO(Vinay): UNSOURCED).

**BC_IC**

| tag | value | unit | attribution |
|---|---|---|---|
| `specific_body_force` | 0 0 | m/s² | TODO(Vinay): UNSOURCED |
| `parameter[sigma0]` ×4 | 0, 0, 0, 0 | Pa | TODO(Vinay): UNSOURCED |
| `parameter[displacement0]` | 0 0 | m | TODO(Vinay): UNSOURCED |
| `parameter[pressure_ic]` | −1e6 | Pa | TODO(Vinay): UNSOURCED |
| `parameter[dirichlet0]` | 0 | m | TODO(Vinay): UNSOURCED |
| `parameter[top_pressure]` | 0 | Pa | TODO(Vinay): UNSOURCED |
| `parameter[bottom_pressure]` | 1e4 | Pa | TODO(Vinay): UNSOURCED |

Notes: `E = 40e6` is unique to this deck in the whole `Tests/Data/RichardsMechanics` tree. **The [F] Bishop deviation was introduced, not inherited:** at the root import the ancestor `beacon_1b_notebook_smoke.prj` carried `BishopsSaturationCutoff cutoff_value=1`; commit `8712c6d7ee` (2026-04-15) replaced it with `BishopsPowerLaw exponent=1`, predating the 2026-05-29 standing decision. This deck appears in audit bucket **H LOW** (σ0 = 0 under −1 MPa suction) — unlike stressprobe, which the audit deliberately excluded.

---

### 3.5 `beacon_1c_dsm_micromacro_smoke.prj` — two media (0 = block, 1 = pellet)

**MATERIAL — global `potential_exchange` block**

| tag | value | unit | attribution |
|---|---|---|---|
| `mass_exchange_coefficient` | 1e-13 | undocumented | TODO(Vinay): UNSOURCED |
| `micro_solid_density_reference` (global, both media inherit) | 2650 | kg/m³ | TODO(Vinay): UNSOURCED |
| `micro_solid_volume_fraction_reference` (global) | 0.6 | dimensionless | TODO(Vinay): UNSOURCED |
| `initial_micro_water_content` (global) | 0.1 | dimensionless | TODO(Vinay): UNSOURCED |
| `micro_liquid_density_reference` (global + both media) | 1300 | kg/m³ | TODO(Vinay): UNSOURCED |
| `micro_liquid_density_a` (global + both media) | 1.3 | dimensionless | TODO(Vinay): UNSOURCED |
| `micro_liquid_density_b` (global + both media) | 1.0 | dimensionless | TODO(Vinay): UNSOURCED |
| `medium[0]/micro_solid_volume_fraction_reference` | 0.65 | dimensionless | TODO(Vinay): UNSOURCED |
| `medium[0]/initial_micro_water_content` | 0.08 | dimensionless | TODO(Vinay): UNSOURCED |
| `medium[1]/micro_solid_volume_fraction_reference` | 0.45 | dimensionless | TODO(Vinay): UNSOURCED |
| `medium[1]/initial_micro_water_content` | 0.14 | dimensionless | TODO(Vinay): UNSOURCED |

**vdW pair:** `hamaker_constant` 1e-30 J, `specific_surface` 1.0 (global) / 0.8 (block) / 1.5 (pellet) m²/kg → NUMERICAL placeholder per §2. Per-medium overrides do **not** restore physicality (Sa³ scaling leaves μ_lR at 2.1e-22 and 8.5e-23 J/kg respectively).

**MATERIAL — medium[0] block / medium[1] pellet**

| tag | block | pellet | unit | attribution (both) |
|---|---|---|---|---|
| `AqueousLiquid/viscosity` | 1e-3 | 1e-3 | Pa·s | TODO(Vinay): UNSOURCED |
| `AqueousLiquid/density` | 1e3 | 1e3 | kg/m³ | TODO(Vinay): UNSOURCED |
| `Solid/density` | 2780 | 2780 | kg/m³ | TODO(Vinay): UNSOURCED |
| `swelling_pressures` | 1.3e7 ×3 | 0.8e7 ×3 | Pa | TODO(Vinay): UNSOURCED |
| `swelling exponents` | 1 1 1 | 1 1 1 | dimensionless | TODO(Vinay): UNSOURCED |
| `lower/upper_saturation_limit` | 0 / 1 | 0 / 1 | dimensionless | TODO(Vinay): UNSOURCED |
| `biot_coefficient` | 1.0 | 1.0 | dimensionless | TODO(Vinay): UNSOURCED |
| `permeability` | 1e-18 | 1e-15 | m² | TODO(Vinay): UNSOURCED |
| `reference_temperature` | 293.15 | 293.15 | K | TODO(Vinay): UNSOURCED |
| `relative_permeability` | 1 | 1 | dimensionless | TODO(Vinay): UNSOURCED |
| `saturation` S_r / S_gr / m / p_b | 0.0 / 0 / 0.5 / 1e6 | 0.0 / 0 / 0.5 / 5e5 | –,–,–,Pa | TODO(Vinay): UNSOURCED |
| `saturation_micro` S_r / S_gr / m / p_b | 0 / 0 / 0.5 / 1e7 | 0 / 0 / 0.5 / 5e6 | –,–,–,Pa | TODO(Vinay): UNSOURCED |
| `bishops_effective_stress/exponent` | 1 | 1 | dimensionless | TODO(Vinay): UNSOURCED |

**MATERIAL — parameters**

| tag | value | unit | attribution |
|---|---|---|---|
| `parameter[E]` | 50e6 | Pa | TODO(Vinay): UNSOURCED |
| `parameter[nu]` | 0.2 | dimensionless | TODO(Vinay): UNSOURCED |
| `parameter[phi0_block]` | 0.35 | dimensionless | TODO(Vinay): UNSOURCED |
| `parameter[phi_tr0_block]` | 0.25 | dimensionless | TODO(Vinay): UNSOURCED |
| `parameter[phi0_pellet]` | 0.675 | dimensionless | TODO(Vinay): UNSOURCED |
| `parameter[phi_tr0_pellet]` | 0.55 | dimensionless | TODO(Vinay): UNSOURCED |

**BC_IC**

| tag | value | unit | attribution |
|---|---|---|---|
| `specific_body_force` | 0 0 | m/s² | TODO(Vinay): UNSOURCED |
| `parameter[sigma0]` ×4 | 0, 0, 0, 0 | Pa | TODO(Vinay): UNSOURCED |
| `parameter[displacement0]` | 0 0 | m | TODO(Vinay): UNSOURCED |
| `parameter[pressure_ic]` | −1e6 | Pa | TODO(Vinay): UNSOURCED |
| `parameter[dirichlet0]` | 0 | m | TODO(Vinay): UNSOURCED |
| `parameter[top_pressure]` | 1e4 | Pa | TODO(Vinay): UNSOURCED |
| `parameter[bottom_pressure]` | 0 | Pa | TODO(Vinay): UNSOURCED |

Notes:
- **The only in-file "citation" in the whole family is inadmissible.** Lines 22-24 attribute the micro-EOS metadata (1300 / 1.3 / 1.0) to "the 1a01 dsm_micromacro-role alignment case". `beacon_1a01` has **no §12 header**. §12.1 forbids exactly this ("do not invent or copy from another PRJ without re-citing"). Circular; both ends UNSOURCED.
- **Micro EOS is ACTIVE** (`a = 1.3`), not the `1e-16` bypass of the MS33 form-(a) decks — a physics difference, not only a provenance gap.
- `phi0_pellet = 0.675` is **close to but not equal to** the MS33 pellet value 0.676259 (= 1 − 900/2780). Different number → UNSOURCED. Do not round into a match.
- **Per-material parametrization question for Vinay:** block and pellet currently share one grain density (2780), one E, one ν, one biot, one Bishop law, one liquid phase. Whether that is intended must be asked, not assumed. Nagra NAB 14-053 (Gaus et al. 2014) tabulates 2,700 kg/m³ for bentonite **pellets** and **blocks** separately.
- Registered as **two** ctests (`Tests.cmake:176`, `:186`) with a 9-field reference-VTU comparison (`:195-203`).

---

## 4. CROSS-PRJ CONSISTENCY CHECK

Values appearing in more than one deck. Attribution must be identical across all occurrences.

| value / tag | decks | attribution | consistent? |
|---|---|---|---|
| α_M `mass_exchange_coefficient` = 1e-13 | all 5 | TODO(Vinay): UNSOURCED | ✅ |
| `micro_solid_density_reference` = 2650 kg/m³ | all 5 | TODO(Vinay): UNSOURCED | ✅ |
| Solid `density` = 2780 kg/m³ | all 5 (6 instances; 1c ×2) | TODO(Vinay): UNSOURCED | ⚠️ **DISAGREEMENT — see below** |
| `biot_coefficient` = 1.0 | all 5 | TODO(Vinay): UNSOURCED | ✅ |
| `bishops_effective_stress/exponent` = 1 (BishopsPowerLaw) | all 5 | TODO(Vinay): UNSOURCED | ✅ |
| liquid `viscosity` = 1e-3 Pa·s | all 5 | TODO(Vinay): UNSOURCED | ✅ |
| liquid `density` = 1e3 kg/m³ | all 5 | TODO(Vinay): UNSOURCED | ✅ |
| `reference_temperature` = 293.15 K | all 5 | TODO(Vinay): UNSOURCED | ✅ value; ⚠️ **classification** (MATERIAL in 4 traces, BC_IC in 1c) |
| `relative_permeability` = 1 | all 5 | TODO(Vinay): UNSOURCED | ✅ |
| VG `exponent` = 0.5 (macro + micro) | all 5 | TODO(Vinay): UNSOURCED | ✅ |
| `nu` = 0.2 | all 5 | TODO(Vinay): UNSOURCED | ✅ |
| swelling `exponents` = 1 1 1 | all 5 | TODO(Vinay): UNSOURCED | ✅ |
| swelling limits 0 / 1 | all 5 | TODO(Vinay): UNSOURCED | ⚠️ **classification** (MATERIAL in 4 traces, NUMERICAL in 1b) |
| `E` = 50e6 Pa | 1a01 smoke, stressprobe, inflow, 1c | TODO(Vinay): UNSOURCED | ✅ (1b uses 40e6, also UNSOURCED) |
| `permeability` 1e-18 m² | 1a01 ×3, 1c block | TODO(Vinay): UNSOURCED | ✅ (1b 5e-16, 1c pellet 1e-15, both UNSOURCED) |
| `swelling_pressures` 1.3e7 Pa | 1a01 ×3, 1c block | TODO(Vinay): UNSOURCED | ✅ (1b 1.1e7, 1c pellet 0.8e7, both UNSOURCED) |
| retention `p_b` 1e6 / 1e7 Pa | 1a01 ×3, 1c block | TODO(Vinay): UNSOURCED | ✅ (1b + 1c pellet use 5e5 / 5e6, also UNSOURCED) |
| micro EOS 1300 / 1.3 / 1.0 | stressprobe, 1c | TODO(Vinay): UNSOURCED | ✅ (circular provenance in both) |
| `initial_micro_water_content` 0.1 | 1a01 smoke, stressprobe, 1b, 1c global | TODO(Vinay): UNSOURCED | ✅ (inflow 0.01; 1c media 0.08/0.14) |
| `pressure_ic` −1e6 Pa | 1a01 smoke, stressprobe, 1b, 1c | TODO(Vinay): UNSOURCED | ✅ (inflow +2e3) |
| `hamaker_constant` | 1e-30 (3 smoke) / 6e-20 (stressprobe) / 5.1e-21 (inflow) | placeholder / UNSOURCED / UNSOURCED | ⚠️ **three values, one family** |
| `specific_surface` | 1.0 (3 smoke, +0.8/1.5 in 1c) / 100 (stressprobe) / **523 (inflow)** | placeholder / UNSOURCED / **EPFL-Seiphoori 2014** | ⚠️ **same tag, one citation, four TODOs** |

### ⚠️ Flag 1 — solid density 2780 kg/m³: verifiers disagree

Four independent verifications **REFUTED** the EURAD-MS attribution (1a01 smoke, stressprobe, 1b, 1c ×2); one **CONFIRMED** it (inflow). **Ledger resolution: `TODO(Vinay): UNSOURCED` in all six instances.** Grounds:

1. **The sibling citation has no locator.** `ANCHORS_MS33_ModelI/ms33_modelI_dd1600.prj:50` reads "solid density rho_s=2780 : EURAD-2 MS33 spec" — no section, table, or figure. Unverifiable under §1.1 item 2 / §12.2.
2. **The named table does not contain the quantity.** `/Users/vinaykumar/tex/eurad2_MS34/MSXX/c_theoretical_benchmarking/theoretical_benchmarking.tex` — Table `tab:TH-G_model_parameters` (caption L108, body L110-130) holds only hydraulic entries (P0=27 MPa, α=0.45, k0=5.6e-21, φ0=0.42, rel-perm A/n, tortuosity, D). Solid density appears **only in prose at L102**: "The solid density is considered to be 2.78~g/cm³." The string `2780` occurs **nowhere** in the spec.
3. **Chronology is backwards.** The beacon decks carry 2780 from the root import (2026-04-10); `ms33_modelI_dd1600.prj` acquires 2780 at `35a1009d3c` (2026-05-20) and its header line first exists at `33eb1cac8e` (2026-05-31). A header six weeks younger cannot be the provenance.
4. **A competing transmission path exists and is documented.** Upstream OGS `Tests/Data/RichardsMechanics/DoubleStructureBenchmark/double_porosity_swelling_RM.prj:87` carries `<value>2780</value>`, and its whole Solid block (`SaturationDependentSwelling`, `13e6 13e6 13e6` ≡ 1.3e7, `exponents 1 1 1`, limits 0/1) plus `nu=0.2` is reproduced verbatim in the beacon decks. Merged to master in `550b3606c8` on **2025-08-11**, eight months before the beacon decks. **An upstream deck with no §12 header is a forbidden source under §12.1.**
5. **The standing audit already classified it.** `DSM_BENCHMARK_DEVIATION_AUDIT_2026-06-06_findings.json`, bucket "C §12.1 material-param source", severity HIGH, listing "solid density=2780" among smoke-test placeholders; report line 200 for 1c explicitly demands separate "pellet vs block" attribution.
6. **Internal incoherence.** Every beacon deck sets `micro_solid_density_reference = 2650` alongside solid density 2780 — two grain densities in one file, 4.7% apart. The genuine MS33 decks are coherent (2780.0 in both slots).
7. **The literature for this configuration disagrees.** Nagra NAB 14-053 (Gaus et al. 2014) gives 2,700 kg/m³ for bentonite blocks *and* pellets; GRS-202 gives 2700 (Serrata) / 2670 (FoCa).

**Available to Vinay, not asserted here:** the value 2.78 g/cm³ = 2780 kg/m³ *is* attested with a real locator at `theoretical_benchmarking.tex:102` for the **MS33/MSXX exercises**. Adopting it for the BEACON decks is a §12.4 cross-family borrowing decision that is Vinay's to make; if made, the sibling header at `ms33_modelI_dd1600.prj:50` should be corrected to carry the L102 locator (and note that L20's implication that it sits in Table TH-G is wrong).

### ⚠️ Flag 2 — `hamaker_constant`: three different values across one benchmark family

1e-30 J (smoke ×3, switch-off placeholder) / 6e-20 J (stressprobe, **outside** the code's own documented DLVO range 1–5e-20 at `PotentialExchange.h:98-101`) / 5.1e-21 J (inflow, **superseded** by the `0d579e8aeb` sweep to 2.2e-20 which missed these decks). No two decks in the family share a Hamaker constant, none equals the branch's approved literature value, and none is sourced.

### ⚠️ Flag 3 — `specific_surface`: one citation, four TODOs, and a 1000× unit question

523 (inflow) is the **only** sourced material literal in the family. 100 (stressprobe) and 1.0 / 0.8 / 1.5 (smoke decks) are different numbers and **must not inherit the Seiphoori citation**. See §6 item 20 for the unit defect that qualifies even the surviving citation.

### ⚠️ Flag 4 — classification inconsistencies to normalise before transcription

`reference_temperature` (MATERIAL vs BC_IC) and swelling `lower/upper_saturation_limit` (MATERIAL vs NUMERICAL) were classified differently by different traces. Pick one convention and apply it to all five headers.

---

## 5. WHAT IS ACTUALLY UNSOURCEABLE

Parameters for which **no §12.1-family source exists anywhere** — not in the PRJ tree, not in git history, not in the local thm-lit / Mendeley libraries. Vinay must supply a source, or grant a §0.2 exemption declaring the decks non-physical smoke harnesses. Grouped by why they are unsourceable.

**A. No source exists anywhere in the tree — including in the §12-compliant MS33 decks**

1. **Young's modulus E** (50e6 Pa ×4 decks; 40e6 Pa in 1b). No E anywhere in `Tests/Data/RichardsMechanics` carries a §12.1 source. The MS33 working value 52e6 is branded in its own header (`ms33_modelI_dd1600.prj:39-46`): "*** WARNING: UNCITED WORKING VALUE — NOT A §12.1 MATERIAL CONSTANT ***… ABSENT from the EURAD-2 MS33/MSXX spec… and from the thm-lit library." There is nothing to borrow.
2. **Mass-transfer α_M = 1e-13.** The MS33 headers explicitly decline to source it (`ms33_modelI_dd1600.prj:57`, wording varies by file): "phenomenological; not identified by sigma_sw (fast-exchange limit)". Also declared unidentifiable in `project_alpha_M_identifiability.md`. Reuse that wording; do not invent a family.
3. **Micro EOS ρ_l0 = 1300 kg/m³, a = 1.3, b = 1.0** (stressprobe, 1c). MS33 headers still read "micro EOS (a, b, …) : TODO(Vinay): primary source". The only in-tree "citation" is circular (1c → 1a01, which has no header). Additionally inconsistent with the memory record that the confined micro-liquid density is ~1100 kg/m³.
4. **Liquid viscosity 1e-3 Pa·s and bulk liquid density 1e3 kg/m³** (all 5). Rounded, **not** the IAPWS-95 values at 293.15 K (1.002e-3 Pa·s, 998.21 kg/m³). A NIST/IAPWS locator is obtainable but **exists in no artifact today**; no headered sibling carries a fluid-property group at all.
5. **Reference temperature 293.15 K** (all 5). Present in MS33 bodies but in no MS33 header. Inert in these isothermal decks; still uncited.

**B. BEACON-specific — and no BEACON specification exists to cite**

`thm_search` for the BEACON 1a01 / MX-80 swelling-pressure benchmark and `mendeley_search` for a BEACON deliverable both returned **no BEACON document** (only unrelated DECOVALEX / GRS / Nagra hits). The §12.1 "Beacon" family therefore **cannot be written into any header** without Vinay supplying the WP/deliverable + table locator himself.

6. **Swelling law** — `swelling_pressures` 1.3e7 / 1.1e7 / 0.8e7 Pa, `exponents 1 1 1` (§1.1 bans uncited constitutive-law exponents explicitly), `lower/upper_saturation_limit` 0 / 1.
7. **Intrinsic permeability** 1e-18 / 5e-16 / 1e-15 m². Three orders from the MS33 anchor (k_ref = 5.6e-21 m² at φ_ref = 0.42) and a different law type — no borrowing possible.
8. **Macro and micro retention sets** — S_r, S_gr, VG m = 0.5, p_b = 1e6/1e7 and 5e5/5e6 Pa. MS33's cited set is P0 = 27 MPa with α = 0.45, a different parametrisation *and* different numbers.
9. **Relative permeability = 1** (law disabled). MS33's `S_e³` citation ("EURAD-2 MS33 spec, A=1, n=3") does not transfer to a constant law.
10. **Initial porosity / transport porosity splits** — 0.4/0.3; 0.5/0.4; 0.35/0.25 (block); 0.675/0.55 (pellet). None is derived in-file from ρ_d/ρ_s; none of the decks declares a target dry density at all.
11. **Micro-structure reference state** — ρ_SR^m = 2650 kg/m³ (all 5; occurs **nowhere else** in `Tests/Data/RichardsMechanics`), nS_ref = 0.6 / 0.65 / 0.45, n_l0 = 0.1 / 0.01 / 0.08 / 0.14. *Trap:* 2650 kg/m³ does appear in GRS-202 as the grain density of **Boom Clay** — wrong material, and GRS is not a §12.1 family. Do not use it.
12. **Solid density 2780 kg/m³** (block *and* pellet) — see §4 Flag 1. Requires Vinay's ruling: adopt the `theoretical_benchmarking.tex:102` locator as an explicit §12.4 cross-family borrow, or leave UNSOURCED. Separately: **do block and pellet share one grain density?** — an ask, not an assumption.
13. **Poisson ratio 0.2** (all 5). The only cited ν in the tree is 0.3 ("EURAD-MS / CIMNE-UPC, CIMNE-UPC.tex Tab. tab:CIMNE_BBM L283"). Different value → citation does not transfer.
14. **Hamaker constant, live values** — 6e-20 J (stressprobe; outside the code's own 1–5e-20 DLVO range) and 5.1e-21 J (inflow; superseded by `0d579e8aeb`, and its only gloss is a bare "literature value" with no author/year/table).
15. **Specific surface 100 m²/kg** (stressprobe) — 3–4 orders below any MX-80 SSA; a scaled probe value, not a measurement.

**C. Modelling assumptions, where "cite a §12.1 family" may be the wrong instrument — Vinay's call**

16. **`biot_coefficient` = 1.0** (all 5). Commit `5a0792dbb0` (2026-06-01) verified verbatim and its diff touches these files; it records an *incompressible-mineral assumption*, not a material measurement, and supplies no §12.1 family. Caveat: the commit is **not an ancestor of HEAD** (reachable only via `deprecated/dsm_native_hierarchical`), which weakens it as a durable §1.1 item-3 citation. A §0.2 note recording the assumption is the honest treatment.
17. **Bishop exponent 1 / `BishopsPowerLaw`** (all 5). Commit `8712c6d7ee` (2026-04-15) is a **one-line commit with no body** that merely restates the value and names no source; it is also a non-ancestor of HEAD. It must **not** be written into a §12.2 header as a citation. It is simultaneously a live audit deviation [F, MEDIUM]: DSM decks are expected to use `BishopsSaturationCutoff(cutoff=1)` per the 2026-05-29 standing decision, and 1b's ancestor originally complied.
18. **Zero-magnitude BC/IC** — `specific_body_force` 0 0, `sigma0` 0×4, `displacement0` 0 0, `dirichlet0` 0. Definitional zeros/structural switches; per the binary rule they are logged `TODO(Vinay): UNSOURCED`, but the realistic remedy is a one-line declaration ("gravity neglected; undeformed reference configuration; fully constrained cell"), not a literature hunt. Note §3's σ0 free-boundary rule does **not** fire in 1a01 smoke/stressprobe/inflow (each boundary constrained in its normal component — rollers, fully confined); it *is* flagged [H, LOW] for 1b and 1c.
19. **Non-zero BC/IC magnitudes with no source** — `pressure_ic` −1e6 Pa (×4) / +2e3 Pa (inflow); `top_pressure` 2e3 / 0 / 1e4 Pa; `bottom_pressure` 1e4 / 0 Pa. Audit bucket **D, HIGH**: "Cite/confirm, or annotate as non-physical smoke placeholders under a §0.2 exemption note." MS33's spec initial suction is 100 MPa (−1e8 Pa) — 100× different, nothing to inherit.

**D. Unit questions that must be settled before any header states a unit**

20. **`specific_surface` unit convention — 1000× discrepancy, affects the MS33 gating decks too.** `PotentialExchange.h:257` documents `h = n_l / (nS · ρ_SR · Sa)` [m], which requires Sa in **m²/kg**. Seiphoori Tab. 1 gives 523 **m²/g** = 5.23e5 m²/kg. The PRJs feed 523 raw. This is already known in-repo: the deleted `augmented_inflow` header at commit `8daad6f894` lines 11-15 — "specific_surface = 523 m2/g (not converted to SI m2/kg). The code-convention film thickness h_code = n_l/(nS*rho_SR*523) is 1000x larger than physical h_SI" — while `ProcessLib/RichardsMechanics/DSM/LITERATURE_1W_FLOOR.md:16` mislabels the tag "523.0 m2/kg" and `ms33_modelI_dd1600.prj:59` calls it "m2/g". **The surviving Seiphoori citation matches the numeral, not the quantity as consumed.** Vinay must rule on the convention before it is frozen into headers.
21. **α_M unit is undocumented by OGS.** `Documentation/…/t_mass_exchange_coefficient.md` is a bare `\copydoc`. The three traces derived three different units from the residual balance (kg m⁻³; kg·s·m⁻⁵; kg² m⁻³ s⁻¹ J⁻¹). Do not assert a unit in the header until this is settled.

**E. Memory-base defects to correct independently (they will otherwise seed a wrong citation)**

22. Three memory records disagree on 523 m²/g. `feedback_piexact_E52_approved_warning.md:31-33` (EPFL / Seiphoori 2014 / MX-80) is the one the primary source supports and is now verbatim-verified. `ogs_rm_dsm_potential_physics.md:131-132` calls it "523 m2/g of FEBEX montmorillonite" — **wrong bentonite**, a cross-material conflation §12.1(6) explicitly warns against. `project_dsm_run_status_board.md:413` records it "PROPOSED (Madsen/Karnland), §12.1 BLOCKED" — superseded by the 2026-06-15 pin.

---

## 6. Header items that ARE covered (no §12.1 source required)

- **`macro_porosity_floor = 0.0` / `micro_water_content_floor = 0.0`** (all five, stash version only). NUMERICAL, inert. Covered by the in-file justification comment and Vinay's approval of 2026-08-12 (§1.1 item 1), plus the parser requirement since commit `71366ac0d3`. Verified inert in source: `PotentialExchange.h:199-205` (`clamped = (n_l_floor > 0.0) && …`) — floor 0 is byte-identical to the unfloored path. Both keys must be present as a **pair**. **No committed version of any of these decks carries them.**
- All solver tolerances, iteration caps, time-stepping tables, output schedules, integration order, component counts/orders, BC component indices, porosity clamps 0/1, `mass_lumping`, `fd_jacobian_for_exchange`, linear-solver settings: NUMERICAL, out of §12.1 scope.

---

## 7. Guardrails firing (§0.1 announcement for the header author)

- **§12.2 / §12.3 — FLAG, STOP.** All five decks have **no provenance block at all** (`<OpenGeoSysProject>` is followed immediately by `<mesh>`), and all five are **already registered in `Tests.cmake`**, which §12.2 forbids until the block is complete and approved. Recorded as a BLOCKER in `DSM_BENCHMARK_DEVIATION_AUDIT_2026-06-06.md:68, 69, 70, 142, 143` and the findings JSON (cluster `beacon_1a01_family`, `beacon_1b_1c_smoke`, `applies_to` includes `dsm_native_maxwell_conjugate`). Also on record as the 2026-05-29 "beacon family orphaned" incident in CLAUDE.md §12.5.
- **§12.1 — ASK USER.** Every material literal in all five decks except `specific_surface = 523` (inflow) is UNSOURCED. Default action taken: report UNSOURCED. Specifically declined: ν=0.2 ← the CIMNE-UPC locator that covers ν=0.3; viscosity/liquid density ← IAPWS (values differ); 2650 ← the Boom Clay grain density in GRS-202; 6e-20 / 5.1e-21 ← Israelachvili & Adams 1978; 2780 ← the bare "EURAD-2 MS33 spec" string; 0.675 ← the MS33 pellet 0.676259; E=40e6 ← a fuzzy Nagra nab14-053 table hit.
- **§1.1 — ASK USER.** All BC/IC magnitudes (suction ICs, hydraulic Dirichlet values, σ0, body force) are uncited; audit bucket D, HIGH.
- **§5 — predicted ≠ verified.** The vdW-suppression magnitudes in §2 are arithmetic from file values, not run-verified. Any later decision to give the placeholders physical values will require a reference-VTU refresh for the 1a01 smoke, 1a01 inflow, and 1c reference ctests (§12.5).
- **Live numerical deviations to carry into the headers, not paper over:** pressure `abstol = 5e-8 Pa` at a 1 MPa suction scale with no `<reltols>` in all five (audit bucket E; `feedback_ogs_rm_pressure_tolerance`); `BishopsPowerLaw` in place of `BishopsSaturationCutoff(cutoff=1)` in all five (audit bucket F).

**Realistic resolution, as the standing audit itself proposes:** these are non-physical plumbing harnesses derived from the upstream OGS `double_porosity_swelling` benchmark with perturbed drivers. Either (a) Vinay supplies §12.1 locators, (b) the decks are **de-registered** from `Tests.cmake`, or (c) each carries an explicit `GUARDRAIL EXEMPTION §0.2 (YYYY-MM-DD)` block declaring the parameter set non-physical smoke placeholders. That call is Vinay's. Nothing was edited in producing this ledger.