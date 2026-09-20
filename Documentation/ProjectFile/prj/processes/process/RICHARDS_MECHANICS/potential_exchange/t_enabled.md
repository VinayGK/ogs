Master on/off switch for the double-structure (micro/macro) potential-exchange
formulation. `bool`.

Optional. At the top level of
\ref ogs_file_param__prj__processes__process__RICHARDS_MECHANICS__potential_exchange
"&lt;potential_exchange&gt;" it defaults to `false` when omitted; in a
\ref ogs_file_param__prj__processes__process__RICHARDS_MECHANICS__potential_exchange__medium
"&lt;medium&gt;" override it defaults to the enclosing top-level block's
resolved value.

`isPotentialExchangeEnabled` treats a medium as DSM-active only when this
flag is `true` in the `PotentialExchangeParameters` resolved for that medium
(the top-level block, or its own `<medium>` override when one exists for
that material id). If `enabled` is `true` anywhere — at the top level or in
any `<medium>` override — the process additionally requires a
\ref ogs_file_param__prj__processes__process__RICHARDS_MECHANICS__micro_porosity
"&lt;micro_porosity&gt;" block; `OGS_FATAL` otherwise
("potential_exchange.enabled=true requires a `<micro_porosity>` process
block.").
