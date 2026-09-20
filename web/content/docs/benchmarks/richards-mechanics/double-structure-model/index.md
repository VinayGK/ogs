+++
project = ["RichardsMechanics/beacon_1a01_dsm_micromacro_smoke.prj", "RichardsMechanics/beacon_1a01_dsm_micromacro_stressprobe.prj", "RichardsMechanics/beacon_1a01_dsm_micromacro_inflow.prj", "RichardsMechanics/beacon_1b_dsm_micromacro_smoke.prj", "RichardsMechanics/beacon_1c_dsm_micromacro_smoke.prj"]
author = "Vinay Kumar, BGR (Federal Institute for Geosciences and Natural Resources)"
date = "2026-09-20"
title = "Double-structure model for compacted bentonite"
weight = 154
image = ""
+++

The double-structure model (DSM) is an optional extension of the
`RICHARDS_MECHANICS` process for materials, such as compacted bentonite, whose
pore space is organised at two scales. The total porosity is split into a
macro (inter-aggregate) pore space, which carries the usual Richards flow and
the mechanical coupling, and a micro (intra-aggregate) pore space represented
by a single scalar micro water content per integration point. Water is
exchanged between the two scales not by a pressure difference but by a
difference of chemical (mass) potential, and the micro-scale potential
carries a van der Waals / disjoining-pressure contribution that represents
the force between clay platelets at small water-film thickness.

{{< data-link >}}

## Switching it on

The feature is enabled through the optional `<potential_exchange>` subtree of
the `RICHARDS_MECHANICS` process definition (parsed in
`CreateRichardsMechanicsProcess.cpp`,
`ProcessLib/RichardsMechanics/ConstitutiveRelations/PotentialExchange.h` and
`PotentialExchangeParameters.h`). With `enabled = true`, each integration
point solves a local micro state (`n_l`, the micro water content) whose
disjoining potential

$$
\mu_{lR}^{\mathrm{vdW}} = \frac{A\, S_a^3}{6\pi}\,
\frac{n_S^3\,\rho_{SR}^3}{n_l^3\,\rho_{lR}}
$$

(Hamaker constant $A$, specific surface $S_a$, micro solid volume fraction
$n_S$ and micro solid/liquid densities $\rho_{SR}$, $\rho_{lR}$) drives a mass
exchange with the macro liquid phase, $\hat\rho_l = \alpha_M\,(\mu_{LR} -
\mu_{lR})$, with $\alpha_M$ the (phenomenological) mass-transfer coefficient.
An optional exponential augmentation term, scaled by a prefactor `K`
(`potential_augmentation_prefactor`, J/kg, default `0.0`) and a decay length
`lambda` (`potential_augmentation_exponent`, m), can be added to the bare van
der Waals term; setting `K` to its default of `0.0` reduces exactly to the
plain van der Waals form (`PotentialExchange.h:124`). None of the five decks
below configures `potential_augmentation_prefactor`, so all five run in this
default, unaugmented regime. A per-material `<medium id="...">` block below
`<potential_exchange>` lets an individual medium override its own micro
reference state — used for the block/pellet media in
`beacon_1c_dsm_micromacro_smoke.prj`.

An excerpt of the block as used in `beacon_1a01_dsm_micromacro_inflow.prj`,
the one deck of the five with a literature-cited `specific_surface` (see
below):

```xml
<potential_exchange>
    <enabled>true</enabled>
    <pressure_tolerance>0.0</pressure_tolerance>
    <hamaker_constant>5.1e-21</hamaker_constant>
    <specific_surface>523</specific_surface>
    <micro_solid_density_reference>2650</micro_solid_density_reference>
    <micro_solid_volume_fraction_reference>0.6</micro_solid_volume_fraction_reference>
    <initial_micro_water_content>0.01</initial_micro_water_content>
    <local_nonlinear_solve_mode>scalar_microstate_storage_mode</local_nonlinear_solve_mode>
    <fd_jacobian_for_exchange>false</fd_jacobian_for_exchange>
    <micro_potential_convention>negative_attractive</micro_potential_convention>
    <macro_porosity_floor>0.0</macro_porosity_floor>
    <micro_water_content_floor>0.0</micro_water_content_floor>
</potential_exchange>
```

`micro_water_content_floor` and `macro_porosity_floor` are mandatory once the
block is present (a floor of `0.0` reproduces the unfloored behaviour
exactly); they cap the disjoining pressure and the macro-to-micro porosity
split respectively.

## Benchmark cases

Five decks exercise the feature: micro-macro exchange plumbing tests on
axisymmetric cell geometries in the BEACON lineage (cases 1a01, 1b, 1c). Each
deck's own header notes that no BEACON specification document with a citable
work-package or table/figure locator has been identified for its parameter
set, so the family attribution is geometry lineage only, not a specification
citation.

| Case | Geometry | Held / free | What it demonstrates |
| :--- | :--- | :--- | :--- |
| `beacon_1a01_dsm_micromacro_smoke` | Axisymmetric cell, one material | All four boundaries roller-confined (fully confined); suction IC −1 MPa, top pressure 2 kPa; $t_\mathrm{end}=10^3$ s | Micro–macro exchange plumbing with the Pi-path potential numerically switched off (placeholder Hamaker constant/specific surface) — exercises the exchange machinery, not bentonite physics |
| `beacon_1a01_dsm_micromacro_stressprobe` | Same geometry | Same full confinement; $t_\mathrm{end}=10^5$ s | Same plumbing with the van der Waals term active (no augmentation, $K=0$); run-only ctest, no VTU regression comparison |
| `beacon_1a01_dsm_micromacro_inflow` | Same geometry | Fully confined; suction IC equal to the top pressure (zero initial hydraulic drive); $t_\mathrm{end}=10^5$ s | Inflow / resaturation probe with the van der Waals term live |
| `beacon_1b_dsm_micromacro_smoke` | Axisymmetric cell, one material, softer and more permeable than 1a01 | Suction IC −1 MPa; through-flow drive (top 0 Pa, bottom 10 kPa); $t_\mathrm{end}=10^3$ s | Exchange-plumbing smoke test on a softer, more permeable medium; van der Waals term off |
| `beacon_1c_dsm_micromacro_smoke` | Axisymmetric cell, two materials (block and pellet) across four medium blocks | Suction IC −1 MPa; top 10 kPa, bottom 0 Pa; $t_\mathrm{end}=10^3$ s | Two-material exchange-plumbing smoke test; van der Waals term off |

## Where the parameters come from

Only one material parameter across the five decks carries a literature
citation in its header. Reproduced here exactly as the header states it:

- **Specific surface** `specific_surface` = 523 — attributed by the
  `beacon_1a01_dsm_micromacro_inflow.prj` header to "EPFL. Seiphoori, Ferrari
  and Laloui (2014), Geotechnique 64(9):721-734, Table 1, p.724, 'Specific
  surface area, s (m2/g)' = 523, MX-80 Wyoming granular bentonite (Table 1
  footnote: data from Ploetze and Weber 2007), doi:10.1680/geot.14.P.017."
  The same header flags an unresolved unit inconsistency for this value (see
  the TODO note below).

Every other dimensioned quantity in the five decks — elastic modulus,
porosity split, intrinsic permeability, retention-curve parameters
(macro and micro), the `SaturationDependentSwelling` law, solid and
micro-reference densities, fluid properties, `biot_coefficient` and the
mass-transfer coefficient $\alpha_M$ — carries no literature citation in any
of the five headers; each header records its own values as uncited working
values (in the headers' own words: "no literature source identified",
"not sourced to a citable reference," "uncited working value," or "not tied
to a cited source," depending on the deck). The Hamaker constant
`hamaker_constant` is likewise uncited in all five decks and is a
placeholder used to switch the van der Waals term numerically on
(`beacon_1a01_dsm_micromacro_stressprobe.prj`, 6e-20 J;
`beacon_1a01_dsm_micromacro_inflow.prj`, 5.1e-21 J) or off
(`beacon_1a01_dsm_micromacro_smoke.prj`, `beacon_1b_dsm_micromacro_smoke.prj`,
`beacon_1c_dsm_micromacro_smoke.prj`, 1e-30 J), not a measured quantity;
correspondingly, `specific_surface` is also a placeholder (1.0, 100, or a
1.0/0.8/1.5 per-medium split) in every deck except
`beacon_1a01_dsm_micromacro_inflow.prj`.

<!-- TODO(author): the unit of the mass-transfer coefficient alpha_M
     (mass_exchange_coefficient) is not documented in the OGS parameter
     reference, and none of the five deck headers assert one. Confirm
     with Vinay whether a unit should be stated here before this page
     asserts one. -->

<!-- TODO(author): the specific_surface value 523 is cited by the
     beacon_1a01_dsm_micromacro_inflow.prj header as m^2/g (Seiphoori et al.
     2014, Table 1), but PotentialExchange.h derives the film thickness as
     h = n_l/(nS*rho_SR*Sa), which is dimensionally m^2/kg. The deck's own
     header calls this an "unresolved unit inconsistency" and states that
     the numeral 523 is entered unconverted, so the citation matches the
     number but not the physical quantity as consumed, making the derived
     film thickness a factor of 1000 too large. Ask Vinay which unit
     convention this page should state for specific_surface. -->
