# Project Context

Keep this file short. It exists to reduce repeated setup and context gathering in future sessions.

## Current session

- Active task: RichardsMechanics transition work tied to `materialmodels/src/TPM/THMDSMRichardsVK_OGS_RM_transition.tex`.
- Touch points: `ProcessLib/RichardsMechanics/RichardsMechanicsFEM-impl.h`, `Tests/ProcessLib/RichardsMechanics/NotebookSingleIntegrationPoint.cpp`, `Tests/Data/RichardsMechanics/beacon_1a01_vk_notebook_roles.prj`, and the generated `Tests/Data/RichardsMechanics/beacon_1a01_vk_stressprobe*` outputs.
- Preferred build tree: `build/release-omp-mfront` when available; it is not present in this checkout right now.
- Verification: run the smallest relevant `ctest` target from that build tree when it exists, and use `git diff` to keep the process file and test data changes aligned.

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
