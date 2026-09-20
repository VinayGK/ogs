Reference micro-scale (grain-level) solid density \f$\rho_{SR}\f$.

Type: `double`, read as `std::optional<double>` through the getter shared
with `hamaker_constant`, `specific_surface`,
`micro_solid_volume_fraction_reference`, `micro_liquid_density_reference`,
`micro_liquid_density_a` and `micro_liquid_density_b`
(`ProcessLib/RichardsMechanics/CreateRichardsMechanicsProcess.cpp`).

Requiredness: when `enabled` = `true`, the selected value (this tag if
given, else an inherited `<medium>`-override default, else `0.0`) must be
`> 0`, else `OGS_FATAL("RichardsMechanics: {} micro_solid_density_reference
must be > 0, got {:g}.")`. When `enabled` = `false`, the tag is optional
and, if given, must still be `> 0`
(`OGS_FATAL("... must be > 0 if provided, got {:g}.")`); the selected value
defaults to `0.0` when omitted. A `<medium>` override inherits the
enclosing block's value when omitted.

Consumed as \f$\rho_{SR}\f$ in
`ProcessLib/RichardsMechanics/RichardsMechanicsFEM-impl.h` and
`ProcessLib/RichardsMechanics/PotentialExchangeParameters.h`:
- In the denominator of the gravimetric water content
  \f$\omega_l = n_l\,\rho_{lR}/(n_S\,\rho_{SR})\f$ of the confined
  micro-liquid-density EOS (`computeReducedMicroLiquidDensity()`).
- As the `rho_SR` argument of the van-der-Waals micro-potential evaluator
  `computeVanDerWaalsMicroPotential()`.
- In the live-\f$K(\rho_d)\f$ evaluation path, where the current dry
  density is computed as \f$\rho_d = \rho_{SR}\,(1-\phi)\f$
  (`effectiveAugmentationPrefactor()` and its \f$\phi\f$-derivative
  companion).
- The full micro-disjoining-pressure swelling-stress increment
  (`computeReferenceMicroPorositySwellingStressIncrement()`) OGS_FATALs
  up front, independently of the `enabled`-branch check above, unless
  `hamaker_constant > 0`, `specific_surface > 0` AND
  `micro_solid_density_reference > 0`; the same three-way check is
  repeated for the `film_energy_route = "exact"` branch of that function.
