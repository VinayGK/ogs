The material's dry density \f$\rho_d\f$ (initial/target rho_d, kg/m^3, per
the source comment).

Type: `std::optional<double>`, read directly with
`config.getConfigParameterOptional<double>("dry_density")`
(`ProcessLib/RichardsMechanics/CreateRichardsMechanicsProcess.cpp`); no
`> 0` or other numeric validation is applied to this tag itself. If the
tag is absent and a `<medium>` override's `defaults` block is present,
the value is inherited from `defaults->dry_density`.

Consumed only at parse time, to resolve the scalar augmentation prefactor
K from a `<potential_augmentation_prefactor_vs_dry_density>` table (when
one is given and the live dry-density mode is not selected):
`potential_augmentation_prefactor_vs_dry_density->getValue(*dry_density)`.
Parsing `OGS_FATAL`s if such a table is given, this tag is absent, and
the live dry-density mode is also not selected. Outside that table
resolution it is not read anywhere else: it does not appear in
`ProcessLib/RichardsMechanics/RichardsMechanicsFEM-impl.h`.
`ProcessLib/RichardsMechanics/PotentialExchangeParameters.h` stores it
only as `std::optional<double> dry_density`, carried in the struct solely
so that a per-`<medium id>` override can supply its own `<dry_density>`
against a shared, inherited `<potential_augmentation_prefactor_vs_dry_density>`
table.
