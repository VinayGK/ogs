Master switch for the film-pressure coupling of the potential-exchange (DSM)
formulation. `bool`, optional, defaults to `true` (a `<medium>` override that
omits it inherits the enclosing `<potential_exchange>` block's setting).

The bare-Pi film-OFF path is retired: the model is consolidated on the film
coupling (`biot` = `alpha`). If a project file sets this to `false`, the
parser overrides it back to `true` and emits a warning; `false` therefore has
no effect on the run — it can no longer disable the coupling.

When potential exchange is enabled and this flag is (effectively always)
`true`, the disjoining micro potential `mu_lR` is augmented by the
film-pressure term so that `mu_lR = -(Pi - b*p_conf)/rho_lR = -p_film/rho_lR`
in every local micro solve and in the macro exchange, and the swelling stress
takes the eigenstrain form. The eigenstrain Biot coefficient `b` used here is
not a separate film parameter: it is the poroelastic `biot_coefficient`
medium property.

The film-pressure delta is only added when, in addition to this flag, a
finite confining pressure is available to the local context; otherwise the
evaluation is bit-for-bit the plain van-der-Waals potential.
