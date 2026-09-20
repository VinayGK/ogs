Exponent prefactor \f$a_\rho\f$ in the confined micro-scale liquid-density
equation of state, see
\ref ogs_file_param__prj__processes__process__RICHARDS_MECHANICS__potential_exchange__micro_liquid_density_reference "&lt;micro_liquid_density_reference&gt;":
\f[
  \rho_{lR} = \rho_{LR} + \rho_{l0}\,\exp\!\left(-a_\rho\,\omega_l^{\,b_\rho}\right)
\f]
Consumed in `computeReducedMicroLiquidDensity()`
(`ProcessLib/RichardsMechanics/RichardsMechanicsFEM-impl.h`), where the
value is clamped to at least \f$10^{-16}\f$ before use in the exponent and
in the analytic Jacobian of the EOS residual.

Type: `double`, read as `std::optional<double>` through the same shared
getter as `micro_liquid_density_reference`.

Same requiredness rule as `micro_liquid_density_reference`: mandatory
(`> 0`, else `OGS_FATAL("RichardsMechanics: {} micro_liquid_density_a must
be > 0, got {:g}.")`) when `local_nonlinear_solve_mode` =
`"scalar_micro_macro_mass_storage_mode"` and `enabled` = `true`; otherwise
optional (if given, still `> 0`, else `OGS_FATAL("... must be > 0 if
provided, got {:g}.")`), defaulting to `0.0` when omitted.

A `<medium>` override that omits this tag inherits the enclosing
`<potential_exchange>` block's value.
