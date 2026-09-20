Dry-density knots \f$\rho_d\f$ [kg/m\f$^3\f$] of the
\ref ogs_file_param__prj__processes__process__RICHARDS_MECHANICS__potential_exchange__potential_augmentation_prefactor_vs_dry_density "&lt;potential_augmentation_prefactor_vs_dry_density&gt;"
table. `std::vector<double>`, required (no default) whenever the enclosing
table is present.

Must have at least 2 entries and the same length as
\ref ogs_file_param__prj__processes__process__RICHARDS_MECHANICS__potential_exchange__potential_augmentation_prefactor_vs_dry_density__prefactors "&lt;prefactors&gt;",
paired one-to-one; otherwise project-file parsing fails with an error naming
both list lengths.

Outside \f$[\rho_{d,\min}, \rho_{d,\max}]\f$ the table's value is held flat
at the corresponding endpoint (clamped, not extrapolated) — the same
convention for both the K-linear (`getValue`) and log-linear
(`getValueLogLinear`) evaluation paths.
