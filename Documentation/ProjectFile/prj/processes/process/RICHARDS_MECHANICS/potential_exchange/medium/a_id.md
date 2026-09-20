The material id (`MaterialIDs`) that this
\ref ogs_file_param__prj__processes__process__RICHARDS_MECHANICS__potential_exchange__medium
"&lt;medium&gt;" override applies to. `int`.

Required — read with the no-default overload of `getConfigAttribute`, so
every `<medium>` must carry an `id`. `OGS_FATAL` if the same `id` value is
used on two `<medium>` entries within one `<potential_exchange>` block, or
if it does not name a material present in `<media>`.
