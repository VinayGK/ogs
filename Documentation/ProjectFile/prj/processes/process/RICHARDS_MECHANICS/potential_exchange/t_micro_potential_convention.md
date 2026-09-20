Selects the sign convention applied to the micro (van der Waals) potential
\f$\mu_{lR}\f$. `std::string`, converted to the enum
`MicroPotentialConvention`.

Allowed values, exactly as checked by the parser:

- `micro_potential_convention="positive_reduced"` (default) —
  `microPotentialSignFactor` returns `+1.0`.
- `micro_potential_convention="negative_attractive"` —
  `microPotentialSignFactor` returns `-1.0`.

Any other string is rejected with `OGS_FATAL` ("unsupported
potential_exchange micro_potential_convention ..."). Optional; defaults to
`"positive_reduced"` at the top level, and to the enclosing top-level value
in a
\ref ogs_file_param__prj__processes__process__RICHARDS_MECHANICS__potential_exchange__medium
"&lt;medium&gt;" override that omits it.

The resulting sign factor multiplies `mu_lR` in
`computeVanDerWaalsMicroPotential` (both the core van-der-Waals term and its
optional exponential augmentation, additively — see
\ref ogs_file_param__prj__processes__process__RICHARDS_MECHANICS__potential_exchange__hamaker_constant
"&lt;hamaker_constant&gt;"). Per the source comment at its point of use, with
`"negative_attractive"` selected, \f$\mu_{lR} < 0\f$, so the disjoining
pressure \f$\Pi = -\rho\,\mu_{lR} > 0\f$ and the swelling eigenstress
\f$\sigma_{sw} = -\phi_m\,\Pi\f$ comes out compressive (OGS tension-positive
convention), i.e. swelling.
