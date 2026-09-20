`double`, optional, defaults to `1e-8` (a `<medium>` override that omits it
inherits the enclosing `<potential_exchange>` block's value). Must be
strictly greater than `0`, or the run stops with an
<tt>OGS_FATAL("fd_jacobian_perturbation must be &gt; 0, ...")</tt>.

It only has an effect when
\ref ogs_file_param__prj__processes__process__RICHARDS_MECHANICS__potential_exchange__fd_jacobian_for_exchange "&lt;fd_jacobian_for_exchange&gt;"
is `true`. In that case it scales the central-difference step used to
approximate the macro pressure-pressure "direct macro derivative"
(drho_L_hat/dp_L) Jacobian contribution: the step is
`fd_jacobian_perturbation * max(1, |p_L|)`, evaluated about the current
liquid-pressure integration-point value p_L. It does not affect the local
micro 2x2 nonlinear-solve Jacobian, which — when forced to finite difference
by the same
\ref ogs_file_param__prj__processes__process__RICHARDS_MECHANICS__potential_exchange__fd_jacobian_for_exchange "&lt;fd_jacobian_for_exchange&gt;"
flag — uses a fixed relative step of `1e-8` instead of this value.
