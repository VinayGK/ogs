Optional subtree configuring the double-structure (micro/macro)
potential-exchange formulation of RICHARDS_MECHANICS. Read with
`getConfigSubtreeOptional`, so omitting the tag leaves the whole block
unconfigured.

Presence of the block alone does not switch the formulation on; that is done by
\ref ogs_file_param__prj__processes__process__RICHARDS_MECHANICS__potential_exchange__enabled
"&lt;enabled&gt;". If `enabled` resolves to `true` in this block or in any
\ref ogs_file_param__prj__processes__process__RICHARDS_MECHANICS__potential_exchange__medium
"&lt;medium&gt;" override, the process additionally requires a
\ref ogs_file_param__prj__processes__process__RICHARDS_MECHANICS__micro_porosity
"&lt;micro_porosity&gt;" block; otherwise the run stops with
`OGS_FATAL("RichardsMechanics: potential_exchange.enabled=true requires a
<micro_porosity> process block.")`.

Every child tag except the two floors is optional and carries a default. The
two exceptions,
\ref ogs_file_param__prj__processes__process__RICHARDS_MECHANICS__potential_exchange__macro_porosity_floor
"&lt;macro_porosity_floor&gt;" and
\ref ogs_file_param__prj__processes__process__RICHARDS_MECHANICS__potential_exchange__micro_water_content_floor
"&lt;micro_water_content_floor&gt;", are read with the no-default overload here
and must therefore be declared explicitly in this top-level block (an explicit
`0.0`, meaning no floor, is permitted).

Zero or more
\ref ogs_file_param__prj__processes__process__RICHARDS_MECHANICS__potential_exchange__medium
"&lt;medium&gt;" subtrees may follow, each overriding this block for one
material id. They are parsed by the same function as this block, with the
already-resolved top-level parameters supplied as defaults, so any tag a
`<medium>` omits is inherited from here.
