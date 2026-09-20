A `double`, the disjoining-pressure floor \f$n_{l,\min}\f$ (a lower bound on
the micro water content \f$n_l\f$). It MUST be declared in the top-level
\ref ogs_file_param__prj__processes__process__RICHARDS_MECHANICS__potential_exchange
"&lt;potential_exchange&gt;" block (the parser no longer defaults it there); a
`<medium>` override may omit it, in which case it inherits the top-level
value. `OGS_FATAL` unless the value is `>= 0`.

It bounds the water content used **only** inside the van der Waals
disjoining-pressure law (where the potential \f$\Pi \sim 1/n_l^3\f$): the
law's value formulas (\f$\mu_{lR,\mathrm{vdW}}\f$, the augmentation term, and
the film thickness) are evaluated at
\f$n_{l,\mathrm{eff}} = \max(n_l, n_{l,\min})\f$, capping \f$\Pi\f$ at
\f$\Pi(n_{l,\min})\f$ instead of letting it diverge as \f$n_l \to 0\f$. It
does **not** change the global \f$n_l\f$, the mass exchange, or the porosity.
Where the floor clamps (\f$n_l < n_{l,\min}\f$), the value is held constant
in \f$n_l\f$, so its \f$n_l\f$-derivatives are zero there; where it does not
clamp (including the in-struct fallback `0.0`, used only for `<medium>`
inheritance, never as a parse default), every formula and derivative is
byte-identical to the unfloored law. The same non-negativity is re-checked
where the value is consumed, as the disjoining law's `n_l_floor` argument.
