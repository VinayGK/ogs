# Superseded ctest reference VTUs — 2026-08-18

These three files are the pre-Option-C ctest reference VTUs for
`ANCHORS_MS33_ModelI` (dd1400/dd1600/dd1800), superseded 2026-08-18 when
the PRJs were migrated to the Option C (sigma0=0) calibration ratified by
Vinay 2026-08-06.

- Produced under: `potential_augmentation_prefactor` = 45217.0 / 103879.0 /
  266767.3 J/kg, `sigma0` = -1.5e5 Pa (3 of 4 expression slots).
- Superseded by: 46000.0 / 104689.9129 / 265905.06 J/kg at `sigma0` = 0,
  the same K(rho_d) table already live in Models III/IV/VII since
  2026-08-17 (commits 55e39b53b8, c5f50ba34e, 7a13ca874b).
- Reason: the OGS repo's own Model I PRJs had never been migrated to the
  2026-08-06 Option C ruling (campaign `sigma0zero_recal_forma_2026-08-06`,
  executed 2026-08-11) — a gap flagged but left open in
  `ProcessLib/RichardsMechanics/DSM/AGENTS.md` entry #17. This migration
  closes it.
- New references verified: bit-exact match (all 11 declared fields,
  every abs/rel norm via `vtkdiff` at the PRJ's own declared tolerances
  = 0.0) against the doubly-cross-validated campaign record
  (`sigma0zero_recal_forma_2026-08-06/{dd*,verify_dd*}/`, itself
  reproduced a third time in this migration on the same canonical
  binary `maxwell-conjugate-20260602`, md5 c432a156).

Kept per repo convention (§6, never delete tracked VTU/PRJ assets) —
historical record, not live references.
