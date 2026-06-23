# Literature-anchored 1W hydration floor for smectite → `micro_water_content_floor`

Status: IN PROGRESS (2026-06-11). PROPOSED values only — §1.1 approval pending (Vinay).
Agent worklog; appended as findings land.

## 0. Source PRJ parameters (read 2026-06-11)

NOTE: the prompt's path `/Users/vinaykumar/git/ogs/validation_2026-06-09/ebs_task13/...`
does not exist; the actual files live at
`/Users/vinaykumar/ogs-models/EBS/runs/validation_2026-06-09_failed/ebs_task13/`.

From `prj-common/process-common-schema20260610.xml` (`<potential_exchange>` block):

| Tag | Value | Symbol |
|---|---|---|
| `specific_surface` | 523.0 m²/kg | Sa |
| `micro_solid_density_reference` | 2780.0 kg/m³ | rho_SR |
| `micro_solid_volume_fraction_reference` | 0.6296296296296297 | nS |
| `hamaker_constant` | 2.2e-20 J | A |
| `potential_augmentation_prefactor` | 58000 (J/kg) | K |
| `potential_augmentation_exponent` | 7.5e-7 m | lambda |
| `micro_liquid_density_reference` | 100.0 kg/m³ (a=1e-16 → EOS bypassed, rho_lR const) | rho_lR (PRJ) |
| `initial_micro_water_content` | 1.3162e-3 | n_l(t=0) |

Caveat to carry into the cross-check: the PRJ's micro rho_lR is 100 kg/m³ (EOS
bypassed); the task-specified cross-check uses rho_lR ~ 1100 kg/m³ (confined
micro-liquid density, per memory ogs_rm_dsm_potential_physics). Both are
reported below.

## 1. Literature values

### 1.1 d001 ranges per hydration state (XRD, Ferrage-family nomenclature)

Day-Stirrat & Bryndzia (2020), "Hydration behavior by X-ray diffraction
profile fitting of smectite-bearing minerals", E3S Web of Conferences 205,
04009 (ICEGT 2020), doi:10.1051/e3sconf/202020504009 — PDF retrieved
2026-06-11 from e3s-conferences.org. Verbatim (Sec. 3, p. 3, re Fig. 1;
modeling per Ferrage):

> "The basal spacings of smectite (Fig. 1) are as follows: dehydrated
> [0W, d001 9.7-10.2 Å], mono-hydrated [1W, d(001) 11.6-12.9 Å],
> bi-hydrated [2W, d(001) 14.9-15.8 Å], and tri-hydrated [3W+, d(001)
> 18-20 Å] layers."

Same paper, gravimetric water contents (Sec. 5, p. 4):

> "in 1W, 2W, and 3W water layer states water is coordinated as 4, 6,
> and 8 molecules (molar mass 18 g/mol), respectively, and that the
> molar mass of smectite and illite-smectite is 757 g/mol. Yielding
> water contents for 1W, 2W, and 3W states as 9.5%, 14%, and 19%
> (molar mass), respectively."

i.e. 1W gravimetric water content ≈ 0.095 g/g (optional cross-anchor).

These ranges trace to Ferrage, Lanson, Sakharov & Drits (2005),
"Investigation of smectite hydration properties by modeling experimental
X-ray diffraction patterns: Part I. Montmorillonite hydration properties",
American Mineralogist 90(8-9), 1358-1374, doi:10.2138/am.2005.1776.
Norrish (1954), Discuss. Faraday Soc. 18, 120-134, is the classic anchor;
full text not retrievable online today — secondary MD literature anchored
on it reports stable Na-montmorillonite states at 9.7 / 12.0 / 15.5 /
18.3 Å (search result paraphrase, 2026-06-11).

### 1.2 Additional 0W anchor

Sun et al. (2019), "Effect of Layer Charge Density on Hydration Properties of
Montmorillonite: Molecular Dynamics Simulation and Experimental Study",
PMC6720539, https://pmc.ncbi.nlm.nih.gov/articles/PMC6720539/ (accessed
2026-06-11), Sec. 2.2.1: "when the structural unit layer is anhydrous,
c = 0.960 nm" — i.e. d001_dry = 9.60 Å.

Norrish (1954) Discuss. Faraday Soc. 18, 120–134 (doi:10.1039/df9541800120):
full text paywalled (RSC); not directly quotable today. Secondary literature
anchored on it (MD validation chains in the 2026-06-11 search results) gives
stable Na-montmorillonite basal spacings 9.7 / 12.0 / 15.5 / 18.3 Å.
Ferrage et al. (2005) Am. Min. 90, 1358–1374 (doi:10.2138/am.2005.1776) is the
primary modern XRD source behind the Sec. 1.1 ranges; abstract not retrievable
through publisher today (HAL + De Gruyter blocked fetches), cited via
Day-Stirrat & Bryndzia (2020) who use its layer-type ranges.

### 1.3 PROPOSED representative values (within all quoted ranges)

- d001_1W = 12.4 Å  (Ferrage-family range 11.6–12.9 Å; canonical Na-mont 1W)
- d001_dry = 9.7 Å  (range 9.7–10.2 Å; Sun et al. 9.60 Å; Norrish-chain 9.7 Å)

PROPOSED, not approved — §1.1 sign-off by Vinay pending.

## 2. Derivation h_1W → n_l_floor

Geometric assumption (stated, not derived): the interlayer expansion from the
dehydrated to the 1W state is carried by ONE water monolayer shared between
the two facing clay surfaces, so the per-face film thickness is half the
basal-spacing increment:

    h_1W = (d001_1W − d001_dry)/2 = (12.4 − 9.7)/2 Å = 1.35 Å = 1.35e-10 m

Code relation (PotentialExchangeParameters.h:137-139, worktree
dsm_maxwell_jac_parallel_wt — same on the schema20260610 branch):

    h = n_l / (nS · rho_SR · Sa)   [mean water film thickness, m]

Inverting with the Task-13 PRJ's own parameters (Sec. 0):

    n_l_floor = h_1W · nS · rho_SR · Sa
              = 1.35e-10 · 0.6296296296296297 · 2780 · 523
              = 1.35e-10 · 915443.7 [1/m]
              = 1.2358e-4

**PROPOSED: micro_water_content_floor = 1.236e-4** (dimensionless n_l,min).
Sanity: initial_micro_water_content = 1.3162e-3 → floor is 9.4 % of the
initial micro water content (initial film h = 1.44 nm vs floor film 0.135 nm);
the floor only engages deep in dehydration, as intended.

Internal-consistency note (honest caveat): the PRJ's Sa = 523 m²/kg is a
model-scale parameter, far below montmorillonite's physical total (interlayer)
surface (~7.5e5 m²/kg). The conversion is nevertheless exact IN-MODEL because
the code defines h through the same Sa: the floor reproduces h_1W = 1.35 Å in
the disjoining law by construction. The gravimetric cross-anchor closes only
with the physical surface: w_1W ≈ h_1W·Sa_phys·rho_w ≈ 1.35e-10·7.5e5·1100
≈ 0.111 g/g vs the literature 0.095 g/g (Sec. 1.1) — agreement within ~15 %,
supporting the h_1W geometry; do NOT use the PRJ Sa for gravimetric checks.

## 3. Pi-cap cross-check

vdW core in pressure form: Pi_vdW = A·Sa³·nS³·rho_SR³/(6π·n_l³·rho_lR)·rho_lR
= A/(6π·h³) (rho_lR cancels when converting the J/kg potential to Pa):

    h³ = (1.35e-10)³ = 2.460e-30 m³
    Pi_vdW(floor) = 2.2e-20 / (6π · 2.460e-30) = 4.74e8 Pa

Augmentation: mu_aug = K·exp(−h_1W/λ) = 58000·exp(−1.35e-10/7.5e-7)
≈ 58000 J/kg (exponent ≈ 1.8e-4, exp ≈ 1).
  - with rho_lR = 1100 kg/m³ (confined micro-liquid, per memory):
    Pi_aug = 6.38e7 Pa → Pi_cap ≈ 5.38e8 Pa ≈ 0.54 GPa
  - with the PRJ's bypassed-EOS rho_lR = 100 kg/m³:
    Pi_aug = 5.8e6 Pa → Pi_cap ≈ 4.80e8 Pa

Verdict: Pi_cap ≈ 5e8 Pa. The task's target window was ~1e9–1e10 Pa; the
result sits a factor ~2 below the 1e9 lower edge but squarely on the
hydration-force scale (Israelachvili, Intermolecular and Surface Forces,
3rd ed. 2011, Ch. 15: hydration-force amplitudes ~1e8–1e9 Pa at sub-nm
films) and FOUR HUNDRED times below the unphysical ~2e11 Pa the unfloored
law diverges to. Documented honestly: slightly below the stated window,
not outside the physics. The vdW core scales as h^-3, so the cap is
sensitive to the d001 choice: taking d001_1W−d001_dry = 2.3 Å (12.3−10.0)
gives h = 1.15 Å and Pi_vdW = 7.7e8 Pa — same order.

## 4. Files written (2026-06-11)

All under /Users/vinaykumar/ogs-models/EBS/runs/validation_2026-06-09_failed/ebs_task13/:

- prj-common/process-common-litfloor.xml — copy of
  process-common-schema20260610.xml + <micro_water_content_floor>1.236e-4</...>
  with citation comment, placed after <potential_augmentation_exponent>
  (parse is name-based and the parameter is optional-with-default, verified in
  CreateRichardsMechanicsProcess.cpp:601-609).
- stg1/12a_t-13_MCC_0.5MPa_litfloor.prj — copy of 12a_..._schema20260610.prj
  pointing at the new common file.
- stg1/12b_t-13_MCC_6MPa_litfloor.prj — same for 12b.

NOT run — PRJs are run-ready only. Value remains PROPOSED pending §1.1
approval.
