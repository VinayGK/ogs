An optional string selecting how the macro/micro porosity split
phi = phi_M + phi_m is updated. Implemented values are:

- hierarchical split (macro_porosity_update_mode="algebraic_split", default)
- hierarchical split, kept as a legacy config alias
  (macro_porosity_update_mode="additive_macro_porosity_rate_mode")

Any other string is rejected with `OGS_FATAL`
("unsupported potential_exchange macro_porosity_update_mode ...").

\attention As currently implemented, both values drive the *same*
hierarchical-split law

\f[
  \phi_M = \frac{\phi - n_l}{1 - n_l}, \qquad
  \phi_m = (1 - \phi_M)\, n_l ,
\f]

(each clamped into \f$[0,\phi]\f$). Selecting
`additive_macro_porosity_rate_mode` additionally logs a one-time INFO message
noting that it now evaluates the hierarchical split rather than a distinct
additive-rate update; it does not change the computed porosities.
