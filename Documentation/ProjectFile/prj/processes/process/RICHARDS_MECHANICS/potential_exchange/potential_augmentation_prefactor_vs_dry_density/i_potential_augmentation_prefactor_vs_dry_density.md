Optional piecewise-linear table \f$K = K(\rho_d)\f$ giving the augmentation
prefactor
\ref ogs_file_param__prj__processes__process__RICHARDS_MECHANICS__potential_exchange__potential_augmentation_prefactor "&lt;potential_augmentation_prefactor&gt;"
as a function of dry density, in place of a single scalar value. A config
subtree, optional; when absent it is inherited from the enclosing block's
default (`nullptr`/absent at the top level).

Two required child lists:
\ref ogs_file_param__prj__processes__process__RICHARDS_MECHANICS__potential_exchange__potential_augmentation_prefactor_vs_dry_density__dry_densities "&lt;dry_densities&gt;"
and
\ref ogs_file_param__prj__processes__process__RICHARDS_MECHANICS__potential_exchange__potential_augmentation_prefactor_vs_dry_density__prefactors "&lt;prefactors&gt;".
Giving this table together with the scalar
`<potential_augmentation_prefactor>` fails project-file parsing with an
error, as the two are mutually exclusive.

Resolving \f$K(\rho_d)\f$ needs a dry density. The source comment records
that this table is evaluated at the material's
`<dry_density>` (initial/target \f$\rho_d\f$, kg/m\f$^3\f$) unless
\ref ogs_file_param__prj__processes__process__RICHARDS_MECHANICS__potential_exchange__potential_augmentation_prefactor_live_dry_density "&lt;potential_augmentation_prefactor_live_dry_density&gt;"
is `true`; giving the table with neither a `<dry_density>` nor that live
mode fails project-file parsing with an error.

Two resolution paths, per the source comment:

- **Frozen (default)**: at parse time, \f$K\f$ is resolved once via the
  table's K-linear interpolant (`getValue()`, clamped/endpoint-held outside
  the knot range) evaluated at the given `<dry_density>`, and stored back
  into the scalar `potential_augmentation_prefactor`. \f$K\f$ is then
  constant in time, so no Jacobian term for it is introduced.
- **Live** (`potential_augmentation_prefactor_live_dry_density == true`):
  the table is not frozen; it is re-evaluated at every evaluation site that
  has the current total porosity in scope, via the table's log-linear
  interpolant (`getValueLogLinear()`), together with its analytic
  \f$dK/d\phi\f$ tangent (`getSegmentSlopeLogLinear()`).

The shared table set at the top level is inherited by each `<medium>`
override as its default, while each medium may supply its own
`<dry_density>`.
