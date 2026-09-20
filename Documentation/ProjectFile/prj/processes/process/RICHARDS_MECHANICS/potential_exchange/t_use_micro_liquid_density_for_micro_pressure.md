Selects which liquid density is used as the *micro-pressure* density: the
density that converts the specific micro potential \f$\mu_{lR}\f$ into a
pressure \f$\Pi\f$ (disjoining pressure, used both in the swelling
eigenstress and in the micro-hydraulic compatibility output). `bool`.

- `true` (default) — use the confined micro-liquid density `rho_lR` from the
  micro EOS.
- `false` — use the bulk liquid density `rho_LR` instead. Not fatal, but
  logged with a `WARN` ("use_micro_liquid_density_for_micro_pressure=false
  selected; micro pressure will use bulk rho_LR instead of confined
  rho_lR.").

Optional. Defaults to `true` at the top level, and to the enclosing
top-level value in a
\ref ogs_file_param__prj__processes__process__RICHARDS_MECHANICS__potential_exchange__medium
"&lt;medium&gt;" override that omits it.

The flag gates the same choice consistently everywhere it is read in
`RichardsMechanicsFEM-impl.h`: the swelling-eigenstress evaluation of
\f$\Pi\f$ is deliberately kept mirrored to the hydraulic micro-pressure
evaluation ("MIRROR the hydraulic p_L_m density choice exactly", per the
source comment at the swelling-stress call site).
