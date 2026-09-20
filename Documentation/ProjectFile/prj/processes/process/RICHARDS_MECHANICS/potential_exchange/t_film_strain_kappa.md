Selects the spacing-strain weighting `kappa` used by the `kinematic` variant
of
\ref ogs_file_param__prj__processes__process__RICHARDS_MECHANICS__potential_exchange__film_strain_coupling "&lt;film_strain_coupling&gt;"
(`d h / d eps_v = kappa*h0`). `std::string`, optional, defaults to
`"aggregate"` (a `<medium>` override that omits it inherits the enclosing
`<potential_exchange>` block's setting). Any other string makes
project-file parsing fail with an error.

Possible values:

- `aggregate` (default): the parser's own error message describes this value
  as `kappa = 1 - phi_M`, "the integrable completion of the eigenstress
  scale". In the implementation, `computeStrainedFilmState` sets
  `kappa = active_nS`, i.e. the active micro solid volume fraction supplied by
  \ref ogs_file_param__prj__processes__process__RICHARDS_MECHANICS__potential_exchange__micro_solid_volume_fraction_mode "&lt;micro_solid_volume_fraction_mode&gt;".
  Those two readings of the symbol are not the same quantity in general, and
  the source does not state which one the tag name is meant to denote.
- `unity`: `kappa = 1` — the parser's error message calls this the "naive
  geometric reading", in which the film spacing follows the REV strain
  one-to-one.

`kappa` is frozen at the integration point (no `d(kappa)/d(eps_v)` chain), and
the effective content fed to the disjoining law is `w_eff = n_l*(1 +
kappa*eps_v)`, floored so the spacing factor stays `>= 1e-6`.

This tag only takes effect through the `kinematic` variant of
`film_strain_coupling` (and, correspondingly, through
\ref ogs_file_param__prj__processes__process__RICHARDS_MECHANICS__potential_exchange__film_energy_route "&lt;film_energy_route&gt;"
`= "exact"`, which requires `kinematic`); the `equilibrium` variant does not
read it.
