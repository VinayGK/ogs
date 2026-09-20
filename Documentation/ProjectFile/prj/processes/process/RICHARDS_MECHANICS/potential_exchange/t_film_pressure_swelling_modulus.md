`double`, optional, defaults to `0.0` (a `<medium>` override that omits it
inherits the enclosing `<potential_exchange>` block's setting); must be
`>= 0`, or project-file parsing fails with an error. Source comment gives
the intended unit and meaning: "eigenstrain modulus K_sw [Pa]; 0 -> drained
K".

Deprecated (source comment): the swelling stress is now computed
as `(1 - phi_M) * p_film`, and this modulus is unused. The value is parsed
and range-checked but, in the current tree, not read by any consuming code —
setting it has no effect on a run.
