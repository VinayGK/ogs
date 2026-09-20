Width of the smooth activation gate that turns on the film-pressure delta as
the confining pressure crosses the gate threshold. `double`, optional,
defaults to `0.0` (a `<medium>` override that omits it inherits the
enclosing `<potential_exchange>` block's setting); must be `>= 0`, or
project-file parsing fails with an error. Unit: Pa (source comment:
"smooth-gate width w [Pa]").

Only used when
\ref ogs_file_param__prj__processes__process__RICHARDS_MECHANICS__potential_exchange__film_pressure_coupling "&lt;film_pressure_coupling&gt;"
is active. Let `x = p_conf - Pi_gate`, with the gate threshold
`Pi_gate = n_S*n_l*Pi`:

- `w == 0` (default): the gate is a sharp Heaviside step in `x` — this
  reproduces the original hard gate exactly.
- `w > 0`: the gate is the C1 smoothstep cubic `g(t) = t^2*(3 - 2t)` with
  `t = x/w` for `0 <= x <= w` (`g = 0` below, `g = 1` above), so the film
  drains continuously rather than snapping on.

As the gate opens (`g > 0`), the film-pressure delta
`+g*b*p_conf/rho_lR` is added to the micro potential `mu_lR`, which raises
`mu_lR` and drives the exchange term `rho_hat_l = alpha*(mu_LR - mu_lR)`
negative, i.e. the micro state drains.
