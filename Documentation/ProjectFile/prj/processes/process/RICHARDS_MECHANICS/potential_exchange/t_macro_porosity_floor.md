A `double`, the macro-porosity floor \f$\phi_{M,\min}\f$ (REV macro
porosity). It MUST be declared in the top-level
\ref ogs_file_param__prj__processes__process__RICHARDS_MECHANICS__potential_exchange
"&lt;potential_exchange&gt;" block (the parser no longer defaults it there); a
`<medium>` override may omit it, in which case it inherits the top-level
value. `OGS_FATAL` unless the value lies in \f$[0, 1)\f$.

It prevents the macro pore from collapsing into the interlayer: the
interlayer water content \f$n_l\f$ is capped at
\f$n_{l,\mathrm{cap}} = (\phi - \phi_{M,\min})/(1 - \phi_{M,\min})\f$, so the
hierarchical split \f$\phi_M \ge \phi_{M,\min}\f$. This is enforced not as a
hard clamp but as a smooth gate that fades the disjoining micro-potential
\f$\mu_{lR}\f$ toward its bulk value as \f$n_l\f$ approaches
\f$n_{l,\mathrm{cap}}\f$ (gate width set by
\ref ogs_file_param__prj__processes__process__RICHARDS_MECHANICS__potential_exchange__macro_floor_cutoff_width
"&lt;macro_floor_cutoff_width&gt;"), so the exchange equilibrates at
\f$n_l \approx n_{l,\mathrm{cap}}\f$ without re-singularising the
pressure-block Jacobian. Beyond the cap the film is treated as saturated and
further water stays in the macro (bulk) pore; the update is porosity- and
water-conserving (\f$\phi = \phi_M + \phi_m\f$ is held).

A value of `0.0` disables the floor: the cutoff is inactive and the
disjoining potential is evaluated unchanged (bit-for-bit identical to the
unfloored form).
