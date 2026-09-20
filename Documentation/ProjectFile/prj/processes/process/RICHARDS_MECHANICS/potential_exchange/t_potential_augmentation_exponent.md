Characteristic film thickness \f$\lambda\f$ [m] for the optional exponential
augmentation term. `double`, optional, defaults to the value inherited from
the enclosing block (`0.0` at the top level).

Must be `> 0` whenever the resolved
\ref ogs_file_param__prj__processes__process__RICHARDS_MECHANICS__potential_exchange__potential_augmentation_prefactor "&lt;potential_augmentation_prefactor&gt;"
\f$K > 0\f$; otherwise project-file parsing fails with an error. When
`potential_augmentation_prefactor == 0` (the augmentation is inactive) this
tag is unconstrained. `computeVanDerWaalsMicroPotential` enforces the same
condition again at evaluation time.

It enters the augmentation term as `lambda` in
`mu_lR_aug = K * exp(-h / lambda)`, and in the dimensionless ratio
`xi = h / lambda = n_l_eff / (lambda * nS * rho_SR * Sa)`, where
`h = n_l / (nS * rho_SR * Sa)` is the mean water film thickness [m].
