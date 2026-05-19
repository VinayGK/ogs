# Project Context

Keep this file short. It exists to reduce repeated setup and context gathering in future sessions.

## Current session

- Active task: DSM_native BEACON/ANCHORS calibration follow-up and reproducible run-state capture on branch `dsm_native`.
- Touch points: `WORKLOG.md`, `PROJECT_CONTEXT.md`, `review_dumps/*`, and `Tests/Data/RichardsMechanics/ANCHORS_MS33_ModelI/*`.
- Preferred build tree: `build/release-omp-mfront` (present in this checkout).
- Verified executable: `/Users/vinaykumar/git/build/release-omp-mfront/bin/ogs` with libs in `/Users/vinaykumar/git/build/release-omp-mfront/lib`.
- Verification completed in this session:
  - full pure-vdW ANCHORS calibration (`1400..1800 kg/m³`, `step=25`, 17 points),
  - full augmented-vdW ANCHORS calibration (`lambda=1e-6`, 17 points),
  - BEACON calibrated inflow run (`beacon_1a01_dsm_micromacro_calibrated_inflow.prj`) and stress comparison vs 604 kPa target.

## Stable facts to record

- Current goal or active branch task.
- Relevant process or subsystem names.
- Active build tree, usually `build/release-omp-mfront`.
- Target tests or `ctest` names.
- Any `materialmodels` interface assumptions: parameter order, output names, generated file names, and test data paths.

## Session workflow

1. Read `AGENTS.md`.
2. Read this file.
3. Inspect only the touched files and related tests.
4. Use `git status` and `git diff` for the active work.
5. Update this file and `WORKLOG.md` only when the stable context changes.

## Good habits

- Prefer file paths and line numbers over pasted context.
- Keep changes in the smallest set of files that satisfies the task.
- For long tasks, keep a short `WORKLOG.md` with the current state and next step.
