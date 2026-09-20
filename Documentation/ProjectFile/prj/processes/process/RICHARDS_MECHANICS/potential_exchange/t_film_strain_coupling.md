Selects the strained-film disjoining law variant, i.e. whether and how the
film spacing (interlayer water content) used to evaluate the disjoining law
is displaced by the volumetric strain. `std::string`, optional, defaults to
`"off"` (a `<medium>` override that omits it inherits the enclosing
`<potential_exchange>` block's setting). Any other string than the three
below makes project-file parsing fail with an error citing
`DSM/STRAINED_FILM_IMPLEMENTATION.md`.

Possible values:

- `off` (default): the film geometry is frozen — bit-for-bit the original
  behaviour.
- `kinematic`: the film spacing follows the volumetric strain,
  `h = h0(n_l)*(1 + kappa*eps_v)`, i.e. the bare disjoining law is evaluated
  at the effective water content `w_eff = n_l*(1 + kappa*eps_v)`. The
  weight `kappa` is set by
  \ref ogs_file_param__prj__processes__process__RICHARDS_MECHANICS__potential_exchange__film_strain_kappa "&lt;film_strain_kappa&gt;".
- `equilibrium`: the film spacing tracks the film force balance once the
  load can compress the film: on the loaded branch (confining pressure
  `p_conf` exceeding the disjoining pressure `Pi` at the current `n_l`),
  `w_eff` is the water content solving `Pi(w_eff) = p_conf`; otherwise
  `w_eff = n_l` (the unloaded branch).

When not `off`, this replaces the shipped integrable mechanical partner
(its frozen-geometry, small-strain truncation of the same coupling) rather
than adding to it, to avoid double-counting the same physics.
