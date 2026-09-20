An optional Jacobian approximation switch for the DSM potential-exchange
contribution only; it does not affect the converged forward solve, only the
Newton tangent. `bool`, defaults to `false` (a `<medium>` override that omits
it inherits the enclosing `<potential_exchange>` block's setting).

When `true` it forces two independent Jacobian terms of the exchange
formulation to be evaluated by central finite difference instead of
analytically:

- the local 2x2 micro nonlinear solve's Jacobian (over micro water content
  and micro liquid density) is taken by a central difference with a fixed
  relative step of `1e-8`, instead of the analytic 2x2 tangent;
- the macro pressure-pressure "direct macro derivative" contribution
  (drho_L_hat/dp_L) is taken by a central difference whose step is scaled by
  \ref ogs_file_param__prj__processes__process__RICHARDS_MECHANICS__potential_exchange__fd_jacobian_perturbation "&lt;fd_jacobian_perturbation&gt;",
  instead of the analytic implicit n_l(p_L) chain-rule tangent.

The two FD switches use different step sizes: the local micro 2x2 one is
fixed at `1e-8` and does not read
\ref ogs_file_param__prj__processes__process__RICHARDS_MECHANICS__potential_exchange__fd_jacobian_perturbation "&lt;fd_jacobian_perturbation&gt;";
only the macro direct-derivative term does.
