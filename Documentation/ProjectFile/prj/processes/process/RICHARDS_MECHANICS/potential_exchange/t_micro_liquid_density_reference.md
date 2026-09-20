Reference micro-scale liquid density \f$\rho_{l0}\f$ entering the confined
micro-scale liquid-density equation of state solved in
`computeReducedMicroLiquidDensity()`
(`ProcessLib/RichardsMechanics/RichardsMechanicsFEM-impl.h`):
\f[
  \rho_{lR} = \rho_{LR} + \rho_{l0}\,\exp\!\left(-a_\rho\,\omega_l^{\,b_\rho}\right),
  \qquad
  \omega_l = \frac{n_l\,\rho_{lR}}{n_S\,\rho_{SR}}
\f]
where \f$\rho_{LR}\f$ is the bulk liquid density, \f$\omega_l\f$ is the
gravimetric water content (the code's own name for this ratio), and
\f$a_\rho\f$, \f$b_\rho\f$ are
\ref ogs_file_param__prj__processes__process__RICHARDS_MECHANICS__potential_exchange__micro_liquid_density_a "&lt;micro_liquid_density_a&gt;"
/
\ref ogs_file_param__prj__processes__process__RICHARDS_MECHANICS__potential_exchange__micro_liquid_density_b "&lt;micro_liquid_density_b&gt;".
The equation is solved for \f$\rho_{lR}\f$ by a Newton iteration on the
residual \f$\rho_{lR}-\mathrm{rhs}(\rho_{lR})\f$ (at most 30 iterations,
convergence tested against \f$10^{-14}\max(1,|\rho_{lR}|)\f$); on
non-convergence the last iterate is used and a one-time warning is logged.

Type: `double`, read as `std::optional<double>` through a getter shared
with `hamaker_constant`, `specific_surface`,
`micro_solid_density_reference`, `micro_liquid_density_a` and
`micro_liquid_density_b`.

Requiredness depends on `local_nonlinear_solve_mode`:
- If `local_nonlinear_solve_mode` = `"scalar_micro_macro_mass_storage_mode"`
  and `enabled` = `true`, the selected value (this tag if given, else an
  inherited `<medium>`-override default, else `0.0`) must be `> 0`, else
  `OGS_FATAL("RichardsMechanics: {} micro_liquid_density_reference must be
  > 0, got {:g}.")`.
- Otherwise the tag is optional: if given it must still be `> 0`
  (`OGS_FATAL("... must be > 0 if provided, got {:g}.")`), and the selected
  value defaults to `0.0` when omitted.

A `<medium>` override that omits this tag inherits the enclosing
`<potential_exchange>` block's value.

With \f$\rho_{l0}=0\f$ (its default outside the mass-storage mode) the EOS
above degenerates to \f$\rho_{lR}\approx\rho_{LR}\f$, i.e. the confined
liquid density collapses onto the bulk value.
