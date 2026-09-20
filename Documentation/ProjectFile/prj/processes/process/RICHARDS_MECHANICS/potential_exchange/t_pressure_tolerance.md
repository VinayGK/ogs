Branch tolerance \f$p_{tol}\f$ for the macro (Young-Laplace) potential
\f$\mu_{LR}\f$. `double`. The source states no unit for this tag; it is
compared directly against the liquid pressure \f$p_{LR}\f$, so it carries
whatever pressure unit that variable does.

`computeYoungLaplaceMacroPotential` treats the liquid pressure \f$p_{LR}\f$
as on the *saturated* branch (\f$\mu_{LR} = 0\f$) whenever
\f$p_{LR} > -p_{tol}\f$, and evaluates \f$\mu_{LR} = p_{LR}/\rho_{LR}\f$
otherwise — i.e. it widens the \f$p_{LR}=0\f$ saturated/unsaturated switch
into a band of half-width `pressure_tolerance` around zero.

Optional. Defaults to `0.0` at the top level; a
\ref ogs_file_param__prj__processes__process__RICHARDS_MECHANICS__potential_exchange__medium
"&lt;medium&gt;" override defaults to the enclosing top-level value when
omitted. `OGS_FATAL` unless the resolved value is `>= 0`
("pressure_tolerance must be >= 0"). When potential exchange is not enabled
for a medium, `getPotentialPressureTolerance` returns `0.0` unconditionally,
irrespective of any value configured here.
