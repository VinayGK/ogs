Selects which closed-form energy route is used to build the strain-dependent
part of the micro potential `mu_lR` when the strained-film law is active.
`std::string`, optional, defaults to `"operational"` (a `<medium>` override
that omits it inherits the enclosing `<potential_exchange>` block's
setting). Any other string makes project-file parsing fail with an error
listing the supported values.

Possible values:

- `operational` (default): the shipped Derjaguin cut, bit-for-bit the
  original behaviour.
- `exact`: the one-Psi energy pair, i.e. `mu_mech` is derived by
  differentiating a single strain-path energy `Psi_film(n_l, eps_v)`
  (closed-form strain integrals of the disjoining law along the `kinematic`
  spacing law) instead of the operational bolt-on term.

`exact` is only admissible together with
\ref ogs_file_param__prj__processes__process__RICHARDS_MECHANICS__potential_exchange__film_strain_coupling "&lt;film_strain_coupling&gt;"
`= "kinematic"` — the closed-form integrals are derived for that spacing
law only. Requesting `exact` with any other `film_strain_coupling` value
makes project-file parsing fail with an error naming the required
`film_strain_coupling` value. The combination predicate
(`isValidFilmEnergyRouteCombination`) is evaluated at parse time only; the
assembly code does not re-check it.
