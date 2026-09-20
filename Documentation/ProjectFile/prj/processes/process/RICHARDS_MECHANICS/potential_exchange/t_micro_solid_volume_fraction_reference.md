Reference micro-scale solid volume fraction \f$n_S\f$.

Type: `double`, read as `std::optional<double>` through the getter shared
with `hamaker_constant`, `specific_surface`,
`micro_solid_density_reference`, `micro_liquid_density_reference`,
`micro_liquid_density_a` and `micro_liquid_density_b`
(`ProcessLib/RichardsMechanics/CreateRichardsMechanicsProcess.cpp`).

Requiredness: mandatory (`> 0`, else
`OGS_FATAL("RichardsMechanics: {} micro_solid_volume_fraction_reference
must be > 0, got {:g}.")`) when `enabled` = `true`; optional (if given,
still `> 0`, else `OGS_FATAL("... must be > 0 if provided, got {:g}.")`),
defaulting to `0.0`, when `enabled` = `false`. A `<medium>` override
inherits the enclosing block's value when omitted.

Consumed as \f$n_S\f$ in
`ProcessLib/RichardsMechanics/RichardsMechanicsFEM-impl.h`:
- The value returned by `computeActiveMicroSolidVolumeFraction()` /
  `computePreviousMicroSolidVolumeFraction()` (clamped to at least
  \f$10^{-16}\f$) whenever
  \ref ogs_file_param__prj__processes__process__RICHARDS_MECHANICS__potential_exchange__micro_solid_volume_fraction_mode "&lt;micro_solid_volume_fraction_mode&gt;"
  = `"reference"`.
- Through that active \f$n_S\f$, the denominator term of the gravimetric
  water content \f$\omega_l = n_l\,\rho_{lR}/(n_S\,\rho_{SR})\f$ of the
  confined micro-liquid-density EOS
  (`computeReducedMicroLiquidDensity()`), whenever `"reference"` mode
  supplies the active \f$n_S\f$.
- The fallback for the micro-saturation reference water content
  \f$n_{l,\mathrm{ref}}\f$ (see
  \ref ogs_file_param__prj__processes__process__RICHARDS_MECHANICS__potential_exchange__initial_micro_water_content "&lt;initial_micro_water_content&gt;")
  used to form \f$S_{L,m} = n_l / n_{l,\mathrm{ref}}\f$ in
  `computeCompatibilityMicroHydraulicOutput()`, when
  `initial_micro_water_content` is not given.
- Required `> 0` by the full micro-disjoining-pressure swelling-stress
  increment (`computeReferenceMicroPorositySwellingStressIncrement()`)
  whenever `micro_solid_volume_fraction_mode` = `"reference"`.
