A per-material override of the
\ref ogs_file_param__prj__processes__process__RICHARDS_MECHANICS__potential_exchange
"&lt;potential_exchange&gt;" block, for one material (matched against
`MaterialIDs` via its
\ref ogs_file_attr__prj__processes__process__RICHARDS_MECHANICS__potential_exchange__medium__id
"id" attribute). Zero or more `<medium>` subtrees may appear inside
`<potential_exchange>`.

The same child tags as for
\ref ogs_file_param__prj__processes__process__RICHARDS_MECHANICS__potential_exchange
"&lt;potential_exchange&gt;" can be used, here — both are parsed by the same
function (`parsePotentialExchangeParameters`). Any child tag omitted in a
given `<medium>` is inherited from the enclosing, already-resolved,
top-level `<potential_exchange>` block.

\attention
\ref ogs_file_param__prj__processes__process__RICHARDS_MECHANICS__potential_exchange__macro_porosity_floor
"&lt;macro_porosity_floor&gt;" and
\ref ogs_file_param__prj__processes__process__RICHARDS_MECHANICS__potential_exchange__micro_water_content_floor
"&lt;micro_water_content_floor&gt;" are the two exceptions: they MUST be
declared in the top-level `<potential_exchange>` block (the parser no longer
defaults them there), but a `<medium>` override MAY omit them, in which case
it inherits the top-level value.

At assembly time, each element's material id is looked up in the resolved
override map: an element whose material id has a `<medium>` override uses
that override's (fully resolved) parameters, and every other element uses
the top-level `<potential_exchange>` block
(`LocalAssemblerInterface::selectPotentialExchangeParameters`).

`OGS_FATAL` if the same `id` is given to two `<medium>` entries within one
`<potential_exchange>` block ("duplicate potential_exchange medium override
for material id ..."), or if `id` names a material that is not present in
`<media>` ("potential_exchange medium override references unknown material
id ...").
