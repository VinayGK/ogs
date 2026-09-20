Specific surface \f$S_a\f$ of the clay platelets, entering the van der Waals
micro potential. `double`.

Read through the same shared helper as
\ref ogs_file_param__prj__processes__process__RICHARDS_MECHANICS__potential_exchange__hamaker_constant "&lt;hamaker_constant&gt;",
with the same required/optional pattern:

- `enabled == true`: the value must be `> 0` — either given here or inherited
  from the enclosing block's default; otherwise project-file parsing fails
  with an error.
- `enabled == false`: optional; if given it must still be `> 0`, otherwise it
  defaults to `0.0`.

In a `<medium>` override the fallback is the value from the top-level
`<potential_exchange>` block; at the top level the fallback is `0.0`.

It enters the micro potential as `Sa` in
`mu_lR_vdW = (A * Sa^3 / (6*pi)) * (nS^3 * rho_SR^3) / (n_l^3 * rho_lR)` [J/kg]
and in the mean water film thickness `h = n_l / (nS * rho_SR * Sa)` [m] used
by the optional exponential augmentation (see
\ref ogs_file_param__prj__processes__process__RICHARDS_MECHANICS__potential_exchange__potential_augmentation_prefactor "&lt;potential_augmentation_prefactor&gt;").
`computeVanDerWaalsMicroPotential` fails with an error if it is called with
`specific_surface <= 0`.

\attention The dimensional-consistency comment on the van der Waals prefactor
in `PotentialExchange.h` annotates this quantity's unit as \f$m^2/kg\f$
(`A · Sa^3 · nS^3 · rho_SR^3 / (n_l^3 · rho_lR)` = `[J] · [m^2/kg]^3 · [1] ·
[kg/m^3]^3 / ([1] · [kg/m^3])` = J/kg). No other unit for `specific_surface`
is stated in this source tree.
