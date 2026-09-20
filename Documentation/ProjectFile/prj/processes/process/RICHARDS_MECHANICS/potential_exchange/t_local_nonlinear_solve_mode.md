Selects which local (integration-point) nonlinear solve/Jacobian path closes
the micro/macro potential exchange. `std::string`, converted to the enum
`LocalNonlinearSolveMode`.

Allowed values, exactly as checked by the parser:

- `local_nonlinear_solve_mode="scalar_exchange"` (default) -> enum
  `ScalarExchange`
- `local_nonlinear_solve_mode="scalar_microstate_storage_mode"` -> enum
  `ScalarReferenceStorage`
- `local_nonlinear_solve_mode="scalar_micro_macro_mass_storage_mode"` -> enum
  `ScalarReferenceMassStorage`

Any other string is rejected with `OGS_FATAL` ("unsupported
potential_exchange local_nonlinear_solve_mode ..."). Optional; defaults to
`"scalar_exchange"` at the top level, and to the enclosing top-level value in
a
\ref ogs_file_param__prj__processes__process__RICHARDS_MECHANICS__potential_exchange__medium
"&lt;medium&gt;" override that omits it.

Selecting `"scalar_micro_macro_mass_storage_mode"` sets the create-time flag
`uses_micro_liquid_density_eos`, which makes the micro liquid-density EOS
parameters (`micro_liquid_density_reference`, `micro_liquid_density_a`,
`micro_liquid_density_b`) required-with-inherited-default in this same block;
they are parsed but left unused by the other two modes. Throughout
`RichardsMechanicsFEM-impl.h` the resolved mode value branches which
micro-state storage variables and Jacobian terms a given integration point
carries.
