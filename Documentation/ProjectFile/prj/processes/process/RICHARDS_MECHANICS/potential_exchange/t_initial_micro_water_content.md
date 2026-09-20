Optional micro-scale water content \f$n_l\f$ reference/initial value.

Type: `std::optional<double>` (kept as an optional; unlike the other
micro-state parameters it is never resolved to a `0.0` fallback at parse
time). Read through the same "must be `> 0` if provided" getter as the
other micro-state parameters, so a given value must be `> 0`, else
`OGS_FATAL("RichardsMechanics: {} initial_micro_water_content must be > 0
if provided, got {:g}.")`. A `<medium>` override defaults to the
enclosing block's value (which may itself be unset) when omitted.

Used in two places in
`ProcessLib/RichardsMechanics/RichardsMechanicsFEM-impl.h`:
- `setInitialConditionsConcrete()`: when `<potential_exchange>` is
  enabled for the medium, the integration point's initial micro water
  content is set to this value if given (subsequently clamped into
  \f$[10^{-12},\,\phi]\f$ with \f$\phi\f$ the initial porosity), in place
  of the value that would otherwise be derived from the current
  total/transport-porosity split, or from the medium's
  `saturation_micro` property when that property is defined.
- `computeCompatibilityMicroHydraulicOutput()`: the reference water
  content \f$n_{l,\mathrm{ref}}\f$ used to form the micro-scale
  saturation \f$S_{L,m} = n_l / n_{l,\mathrm{ref}}\f$ is this value if
  given; when omitted it falls back to
  \ref ogs_file_param__prj__processes__process__RICHARDS_MECHANICS__potential_exchange__micro_solid_volume_fraction_reference "&lt;micro_solid_volume_fraction_reference&gt;"
  (i.e. the reference solid volume fraction \f$n_S\f$, not a water-content
  quantity, is used as the fallback denominator when this tag is
  omitted).
