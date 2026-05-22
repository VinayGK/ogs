# ANCHORS MS33 status snapshot (streamlined)

Authoritative tracker: `Tests/Data/RichardsMechanics/ANCHORS_MS33_ModelIV/AGENTS.md`.

## Scope
- Model set: I / III / IV / VII (+ V_LE/V_MCC where present).
- Branch: `dsm_native_hierarchical`.
- Runtime binary: `/Users/vinaykumar/git/build/release-omp-mfront/bin/ogs`.

## Roadmap summary (one-line refs)
- Step 1: `0d7a9edd64`
- Step 2: `88d42c98fd`
- Step 3: `c4888b6db4`, `ce9178fa96`
- Step 5: `0d579e8aeb`
- Step 6: `66b782afa1`
- Step 7: `4d47efff55`, `ce9178fa96`
- Step 8: `3ac6b7de1f`

## Operational policy
- Keep baseline runs spec-compliant.
- Keep exploratory/non-spec runs isolated from baseline PRJs and committed status artifacts.
- Treat the 3-gate verification set as mandatory after physics changes:
  1. calibration check,
  2. canonical LE replay,
  3. DSM micro-macro test subset.

## Current summary
- Latest DSM test subset: passing.
- Canonical LE replay: unchanged acceptance/rejection status.
- Remaining open work is benchmark-quality alignment, not current execution stability.
