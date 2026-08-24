# dispatch_nodes

Reusable fan-out helper across Vinay's compute nodes (this MacBook Pro,
`macmini`, `shilpas-mbp`). Built 2026-08-24 per conversation; see project
memory `feedback_remote_compute_nodes` for the standing rule and node
roster this codifies.

## What it does

1. **Probe** — reachability per node, plus an HID-idle-time check for nodes
   flagged `respect_active_use` (e.g. `shilpas-mbp` — it's not exclusively
   Vinay's machine; a node in active interactive use is excluded, not
   silently used anyway).
2. **Atomize** — splits a job manifest into per-node batches weighted by
   each available node's physical core count (§6.9 "one run per processor"),
   one job per core per round; extra jobs automatically spill into further
   rounds.
3. **Dispatch** — local jobs run as ordinary background processes; remote
   jobs launch fully detached (`nohup`, stdio closed — no `setsid`, it's
   Linux-only and absent on macOS) so a dropped SSH session can't kill a
   multi-hour OGS run, then get tracked by polling a completion sentinel
   rather than blocking on the launch session itself.
4. **Assimilate** — `keep` paths (PRJs, logs, metrics: small and durable)
   are always rsync'd back into a local snapshot folder and are **never**
   deleted anywhere by this tool. `mirror_discard` paths (bulk scratch,
   e.g. raw VTU output) are rsync'd with checksum verification into a
   local mirror root; the remote copy is deleted **only** with `--discard`
   **and** a byte-count-confirmed mirror.
5. Writes a mechanical `README.md` run-card (branch/commit per node, job
   list, exit codes, timings) into the snapshot folder — facts only, no
   result interpretation. That stays Vinay's call per CLAUDE.md.

**This tool never commits or pushes.** Snapshotting into git (per the
existing §6.8 per-run-snapshot standard) is a separate, explicit,
user-approved step, same as always.

## Usage

```bash
python3 dispatch.py probe
python3 dispatch.py plan --jobs jobs.json
python3 dispatch.py run  --jobs jobs.json --snapshot-name my_sweep_2026-08-24 --dry-run
python3 dispatch.py run  --jobs jobs.json --snapshot-name my_sweep_2026-08-24
python3 dispatch.py run  --jobs jobs.json --snapshot-name my_sweep_2026-08-24 --discard
```

Always `--dry-run` a new manifest first. `--discard` is opt-in per
invocation on purpose — the default is mirror-only, nothing deleted
anywhere, until you've reviewed a dry run and a real run's output.

## Job manifest (JSON array)

```json
[
  {
    "id": "unique-job-id",
    "cmd": "shell command, run via bash -lc inside cwd",
    "cwd": "optional, relative to the node's ogs_repo (default: repo root)",
    "keep": ["relative paths always copied back; never deleted anywhere"],
    "mirror_discard": ["relative paths of bulk output; mirrored, deleted
                         remotely only with --discard + verified mirror"]
  }
]
```

See `example_jobs.json` for a harmless smoke-test manifest (just runs
`hostname`/`sysctl` on each node, nothing OGS-specific, nothing to keep or
discard) — use it to sanity-check connectivity end to end before pointing
this at a real simulation batch.

## Node roster (`nodes.json`)

Edit this file to add/remove nodes or change core counts. Fields:

- `type`: `"local"` or `"remote"`.
- `ssh_alias`: matches `~/.ssh/config` (remote only).
- `cores`: physical core count, used for the weighted split.
- `ogs_repo`: path to the `ogs` checkout on that node.
- `respect_active_use` / `idle_threshold_s`: gate dispatch on the node
  being idle (used for shared machines, not Vinay's own).

New nodes: get host/username/auth details from Vinay directly rather than
guessing or trying multiple credentials (see the shilpas-mbp setup
incident in `feedback_remote_compute_nodes`).

## Safety notes

- `keep` paths are copy-only, forever. The tool has no code path that
  deletes anything listed under `keep`.
- `mirror_discard` deletion requires: `--discard` on the command line,
  the node being `"remote"` (never the local machine), a successful
  checksummed rsync, and a remote-vs-local byte-count match. Any mismatch
  or rsync failure just skips the delete and logs why — it never guesses.
- Every `keep`/`mirror_discard` entry is resolved to an absolute path
  (`resolve_scoped_path`, via `os.path.normpath`) and **rejected unless it
  stays inside the node's `ogs_repo`** — this is the authoritative check,
  and it catches an escape introduced through `cwd` as well as through the
  path itself (e.g. `cwd: "../../.."` combined with an innocuous-looking
  `mirror_discard` entry). `is_safe_discard_path` is a cheap syntactic
  first-line filter on top of that, not a substitute for it.
- `mirror_discard` deletion additionally requires: `--discard` on the
  command line, the node being `"remote"` (never the local machine), a
  successful checksummed rsync, and a remote-vs-local byte-count match
  against the specific path (no falling back to measuring a parent
  directory). Any mismatch, escape, or rsync failure just skips the delete
  and logs why — it never guesses.
- None of the above is exhaustive — review `mirror_discard` entries in any
  manifest before running with `--discard`, same as you'd review an
  `rm -rf` before running it yourself.
- `cmd` runs with full shell privilege in `cwd` (default: the node's live
  `ogs_repo` checkout, not an isolated worktree) — nothing sandboxes it
  from `git checkout`/`git reset --hard`/etc. against whatever's active on
  that node. Point `cwd` at a worktree if a job needs isolation; see
  CLAUDE.md's "work in a fresh worktree" rule.
- Before dispatching real OGS work to a node, check what branch its
  `~/git/ogs` is on — this tool records it in the README but does not
  check it out or switch it for you (see the deprecated-branch caveat in
  `feedback_remote_compute_nodes` / CLAUDE.md).

## Status

First cut, 2026-08-24, then reviewed the same day by an adversarial
multi-lens pass (destructive-path safety, SSH/dispatch correctness,
CLAUDE.md compliance) before being trusted. The review found and this
build fixes: a critical `cwd`-based path-traversal gap that could
`rm -rf` outside the repo (e.g. onto a shared node's home directory), a
stale-exit-sentinel gap that could misreport a network-killed multi-hour
job as a clean success, and the lack of any protection against a dropped
SSH session killing a long remote job. See the module docstring's "REVIEW
FIXES" section for the full list, including a `setsid`-doesn't-exist-on-
macOS bug introduced by the first attempt at the last fix and caught by
smoke-testing afterward, not by the review itself.

Smoke-tested post-fix, on the live roster: `probe`/`plan`/`dry-run`/a real
run via `example_jobs.json`; the mirror→verify→discard path against real
scratch files on `macmini` (mirrored, byte-verified, then deleted
remotely); the `cwd`-traversal attack from the review (confirmed
rejected before any rsync/rm); an absolute-path `keep` escape (confirmed
rejected); missing-required-field and too-broad-`ogs_repo` config
validation (confirmed both raise clearly). All test artifacts were
cleaned up afterward (locally, on the mac mini, and in `~/ogs-models/`)
— that's why those directories are empty on disk, not evidence nothing
was run.

Not yet exercised: a real OGS simulation batch, `--discard` against real
(large) bulk VTU output, or an actual network-drop during a live
multi-hour job (the detach mechanism is verified to survive the
*launching* session closing cleanly; it hasn't been tested against a
mid-run drop specifically). Do those deliberately, watching the log
output, before trusting this fully unattended on a real campaign.
