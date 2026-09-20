Augmentation-prefactor knots \f$K\f$ [J/kg] of the
\ref ogs_file_param__prj__processes__process__RICHARDS_MECHANICS__potential_exchange__potential_augmentation_prefactor_vs_dry_density "&lt;potential_augmentation_prefactor_vs_dry_density&gt;"
table, paired one-to-one with
\ref ogs_file_param__prj__processes__process__RICHARDS_MECHANICS__potential_exchange__potential_augmentation_prefactor_vs_dry_density__dry_densities "&lt;dry_densities&gt;".
`std::vector<double>`, required (no default) whenever the enclosing table is
present, and must have the same length (at least 2 entries) as
`<dry_densities>`.

Every entry must be `> 0`, or project-file parsing fails with an error naming
the offending entry and its value. The source comment explains why this is
enforced unconditionally, for every such table, and not only when
\ref ogs_file_param__prj__processes__process__RICHARDS_MECHANICS__potential_exchange__potential_augmentation_prefactor_live_dry_density "&lt;potential_augmentation_prefactor_live_dry_density&gt;"
is set: the live evaluation path's log-linear interpolant evaluates
`ln(K_r/K_l)` between adjacent knots, which is undefined or sign-corrupting
for a non-positive knot (`K_l == 0` gives `0*inf`; `K_l < 0` gives `NaN`
throughout; `K_r == 0` gives a finite-looking value with a `NaN` slope, the
comment calls this "the silent case, and the dangerous one"). The older
K-linear evaluation path carries no such precondition on its own, but the
check is applied to every table regardless of which path is selected.
