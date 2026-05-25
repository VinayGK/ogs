# DSM MFront Hierarchical Port Log

## Objective
Port hierarchical DSM behavior from `dsm_native_hierarchical` into MFront-idiomatic implementation on top of `dsm_mfront`.

## Branch set
- `master`
- `dsm_native`
- `dsm_native_hierarchical`
- `dsm_mfront`
- `dsm_mfront_hierarchical`

## Core semantic requirement
Reproduce `dsm_native_hierarchical` functionality via MFront-side implementation in the same architectural style by which `dsm_mfront` implements `dsm_native`.

## Chronology
- Initialized branch `dsm_mfront_hierarchical` from `dsm_mfront`.
- Located and inspected `ProcessLib/RichardsMechanics/DSM/DSM_NATIVE_HIERARCHICAL_PATCH_RECIPE.md` as analytical input.

## Replay seed
From clean `master`:
1. `git checkout dsm_mfront`
2. `git checkout -b dsm_mfront_hierarchical`

Further steps appended as commits are made.
