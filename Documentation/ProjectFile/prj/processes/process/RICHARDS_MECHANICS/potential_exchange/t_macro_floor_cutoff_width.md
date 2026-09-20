An optional `double`, the width \f$w\f$ (in \f$n_l\f$, dimensionless) of the
smooth film-to-bulk cutoff gate used by
\ref ogs_file_param__prj__processes__process__RICHARDS_MECHANICS__potential_exchange__macro_porosity_floor
"&lt;macro_porosity_floor&gt;". `OGS_FATAL` unless the value is `>= 0`.
Inherited from the parent block in a `<medium>` override when omitted; the
top-level default is `0.0`.

It only has an effect when `macro_porosity_floor > 0`. If given as `> 0.0`
it is used directly as the gate width. If left at its default `0.0`, the
code instead substitutes \f$\max(10^{-8},\, 0.05\, n_{l,\mathrm{cap}})\f$,
i.e. 5&nbsp;% of the local \f$n_{l,\mathrm{cap}} = (\phi -
\phi_{M,\min})/(1-\phi_{M,\min})\f$.

The gate itself is a smoothstep-like function \f$g = t(2-t)\f$ with
\f$t = (n_{l,\mathrm{cap}} - n_l)/w\f$, chosen so it is \f$C^1\f$ where the
cutoff switches on (\f$t=1\f$) but keeps a nonzero slope at full cutoff
(\f$t \to 0\f$), rather than a smoothstep whose slope also vanishes there.
