Scalar augmentation amplitude \f$K\f$ [J/kg] for the optional lumped
exponential term added to the van der Waals micro potential. `double`
(read as `std::optional<double>` and then resolved).

One of two mutually exclusive ways to set \f$K\f$ — the source comment
labels them "(1) scalar `<potential_augmentation_prefactor>`" and "(2) table
\ref ogs_file_param__prj__processes__process__RICHARDS_MECHANICS__potential_exchange__potential_augmentation_prefactor_vs_dry_density "&lt;potential_augmentation_prefactor_vs_dry_density&gt;"".
Giving both in the same block fails project-file parsing with an error.

Resolution, in order:

- If a `<potential_augmentation_prefactor_vs_dry_density>` table is present,
  this tag must be absent; \f$K\f$ is instead resolved from the table (see
  that tag's page). It is an error to give the table without either a
  `<dry_density>` or
  \ref ogs_file_param__prj__processes__process__RICHARDS_MECHANICS__potential_exchange__potential_augmentation_prefactor_live_dry_density "&lt;potential_augmentation_prefactor_live_dry_density&gt;"
  `== true`.
- Otherwise, this tag's value is used if given, else it defaults to the value
  inherited from the enclosing block (`0.0` at the top level).

The resolved value must be `>= 0`, or project-file parsing fails with an
error.

It enters the micro potential additively: `mu_lR_aug = K * exp(-h / lambda)`
with `h = n_l / (nS * rho_SR * Sa)` the mean water film thickness [m] and
`lambda` the
\ref ogs_file_param__prj__processes__process__RICHARDS_MECHANICS__potential_exchange__potential_augmentation_exponent "&lt;potential_augmentation_exponent&gt;".
Total: `mu_lR = sign * (mu_lR_vdW + mu_lR_aug)` — the augmentation never
replaces the van der Waals term; `K = 0` (the default) reduces exactly to the
pure van der Waals form.
