# Agent instructions for OpenGeoSys-6

**Project**: Scientific THMC simulation framework in modern C++. See [README.md](./README.md) for overview.

## Mandatory code standards

**C++**: C++23 standard (required). Use `ranges-v3` (not `std::ranges`), Eigen (not raw loops). Follow [style guide](https://ufz.github.io/styleguide/cppguide.xml).

**CMake**: 3.31+, target-based commands, no global variables.

**Python**: black-based formatting. ruff check linter.

**Documentation**: British English spelling (colour, behaviour).

## Architecture layers

Foundation:  BaseLib → MathLib → NumLib
Geometry:    GeoLib, MeshLib, MeshGeoToolsLib, MeshToolsLib
Materials:   MaterialLib (MPL), ParameterLib
Processes:   ProcessLib (20+ process implementations)
Apps:        CLI, ApplicationsLib, FileIO, Utils, DataExplorer(Qt)

## Process implementation pattern (REQUIRED)

Every process must have:

1. `{Name}Process.h` - Inherits Process, manages assembly/timestepping
2. `{Name}ProcessData.h` - Material properties, parameters, solver configuration
3. `{Name}LocalAssembler.h` - Element-level assembly (M, K, b matrices)
4. `Create{Name}Process.h` - Factory from XML configuration

## Testing & validation

- Unit tests: `Tests/{LibName}/` (Google Test)
- Integration tests: `Tests/Data/{ProcessName}/` (`.prj` files with reference outputs)
- Always run ctests from release build.
- Check `.clang-format`, `.clang-tidy` for linting rules

## DSM native transition checkpoint

- Active native DSM source tree:
  `/Users/vinaykumar/git/ogs-native-dsm-transition`
- Verified build tree:
  `/Users/vinaykumar/git/build/release-native-transition2`
- As of 2026-05-06, the focused DSM native gate is green:
  `ctest -j 18 --output-on-failure -R "ogs-RichardsMechanics(/DoubleStructureBenchmark/double_porosity_swelling_RM|_.*dsm_micromacro)"`
  passed `32/32`.
- The failing DSM native checks were stale `vtkdiff` references, not solver
  failures and not a tolerance-only issue. Serial and OpenMP regenerated VTUs
  were byte-identical for the failing reference cases.
- The refreshed reference files are:
  - `Tests/Data/RichardsMechanics/beacon_1a01_reference_t_1000.000000.vtu`
  - `Tests/Data/RichardsMechanics/beacon_1a01_dsm_micromacro_inflow_reference_t_100000.000000.vtu`
  - `Tests/Data/RichardsMechanics/beacon_1c_reference_t_1000.000000.vtu`
- DSM source-audit guardrail: keep
  `ProcessLib/RichardsMechanics/ConstitutiveRelations/PotentialExchange.h`
  aligned with the MFront vdW regularisation: use
  `omega3 + omega_min_vdw3` with `h_min_vdw = 5e-11`, and keep the analytic
  derivatives consistent with that denominator.

## 2026-05-11 — Appendix A/B sync protocol (paper-code-presentation)

- Canonical derivation source is Appendix A in the paper:
  `/Users/vinaykumar/tex/dsm-bgr-paper/draft/paper_DSM.tex`
  section *Porosity split and reference-volume consistency*.
- Canonical comment-status source is Appendix B:
  `/Users/vinaykumar/tex/dsm-bgr-paper/draft/appendix_comment_status.tex`.
- Keep DSM-native implementation narrative aligned with Appendix A replacement logic:
  - macro porosity rate follows total-minus-micro coupling,
  - micro evolution is exchange-driven in the local residual,
  - exchange uses `\mu^Macro-\mu^Micro` at consistent REV scale.
- Do not reintroduce ad-hoc `A_H` scaling or ad-hoc `\left(\alpha_\text{Biot}-\phi\right)`
  split statements in docs/scripts when Appendix A equations are the active basis.
- If implementation changes affect Nagel comments #5--#14 scope (reference
  volume, prefactors, exchange scale, missing micro-density-rate term), update
  Appendix B status rows in the same work cycle.
- For Appendix B edits, preserve status taxonomy and validate:
  row count, per-status counts, and `N + M + K = total`.
