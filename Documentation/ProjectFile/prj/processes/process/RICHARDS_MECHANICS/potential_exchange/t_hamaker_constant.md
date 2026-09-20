Hamaker constant \f$A\f$ [J] for the van der Waals micro potential between
clay platelets. `double`.

Read through a shared helper, so its required/optional status depends on
`enabled`:

- `enabled == true`: the value must be `> 0` — either given here or inherited
  from the enclosing block's default; if neither is `> 0`, project-file
  parsing fails with an error ("`hamaker_constant` must be > 0").
- `enabled == false`: optional; if given it must still be `> 0`, otherwise it
  defaults to `0.0`.

In a `<medium>` override the fallback is the value from the top-level
`<potential_exchange>` block; at the top level the fallback is `0.0`.

The source comment on the potential formula
(`ProcessLib::RichardsMechanics::computeVanDerWaalsMicroPotential`) states
that \f$A\f$ is a material constant, not a fitting parameter, and records a
literature value for montmorillonite-water-montmorillonite of
\f$2.2\times10^{-20}\f$ J (Israelachvili & Adams 1978, SFA mica proxy), with a
literature range of \f$1\f$–\f$5\times10^{-20}\f$ J for smectite (DLVO
literature).

It enters the micro potential as
`mu_lR_vdW = (A * Sa^3 / (6*pi)) * (nS^3 * rho_SR^3) / (n_l^3 * rho_lR)` [J/kg],
where `Sa` is
\ref ogs_file_param__prj__processes__process__RICHARDS_MECHANICS__potential_exchange__specific_surface "&lt;specific_surface&gt;",
and `computeVanDerWaalsMicroPotential` itself fails with an error if it is
called with \f$A \le 0\f$.
