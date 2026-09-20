Selects whether the augmentation prefactor \f$K\f$ is re-evaluated live
against the evolving dry density, instead of being frozen at parse time.
`bool`, optional, defaults to `false` (or the value inherited from the
enclosing block).

`false` (the default) is the parse-time freeze: bit-for-bit the existing
behavior, where
\ref ogs_file_param__prj__processes__process__RICHARDS_MECHANICS__potential_exchange__potential_augmentation_prefactor "&lt;potential_augmentation_prefactor&gt;"
is resolved once, when the project file is read.

`true` requires a
\ref ogs_file_param__prj__processes__process__RICHARDS_MECHANICS__potential_exchange__potential_augmentation_prefactor_vs_dry_density "&lt;potential_augmentation_prefactor_vs_dry_density&gt;"
table to be present; otherwise project-file parsing fails with an error. When
set, the table is not frozen: at every evaluation site that has the current
total porosity `phi` in scope, `effectiveAugmentationPrefactor()` re-evaluates
\f$K\f$ at the evolving dry density `rho_d = rho_SR * (1 - phi)`
(`rho_SR = micro_solid_density_reference`) via the table's log-linear
interpolant, with the source-derived analytic tangent `dK/dphi` wired into the
corresponding Jacobian blocks. Evaluation sites without a porosity in scope
fall back to the parse-time scalar
`<potential_augmentation_prefactor>`, which under this mode serves only as
that fallback (resolved, if a `<dry_density>` is given, from the table at that
initial/target density).
