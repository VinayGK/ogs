#!/usr/bin/env python3
"""dispatch.py - reusable fan-out helper across Vinay's compute nodes.

Built for Vinay (BGR / EBS Task Force) per discussion 2026-08-24. Reviewed
2026-08-24 by an adversarial multi-lens pass (destructive-path safety,
SSH/dispatch correctness, CLAUDE.md compliance); see "REVIEW FIXES" below
for what that pass found and what changed as a result.

WHAT THIS TOOL DOES (and does NOT do)
--------------------------------------
It does NOT decide what to run, what a job's expected result is, or whether
a run is "correct." It only mechanises the plumbing already agreed in
conversation and in project memory (feedback_remote_compute_nodes):
  1. PROBES each configured node for reachability, and for nodes flagged
     "respect_active_use" (shared machines, e.g. Shilpa's MBP), for whether
     someone is actively using it (HID idle time) -- busy nodes are excluded,
     never silently overridden.
  2. ATOMIZES a job manifest into per-node batches, weighted by each
     available node's physical core count (one job per core per round;
     leftover jobs roll into a further round automatically).
  3. DISPATCHES each round: local jobs run as ordinary background processes;
     remote jobs are launched fully DETACHED (nohup, stdio closed) so a
     dropped SSH session cannot kill a multi-hour OGS run, then tracked by
     polling for a completion sentinel rather than blocking on the launch
     SSH session itself.
  4. ASSIMILATES results: "keep" paths (small, durable -- PRJs, logs,
     metrics) are always rsync'd back into a local snapshot folder and are
     NEVER deleted anywhere by this tool. "mirror_discard" paths (bulk
     scratch, e.g. raw VTU output) are rsync'd (checksummed) into a local
     mirror root; the remote copy is deleted ONLY if invoked with --discard
     AND the mirror is verified present with a matching byte count AND the
     resolved path is confirmed to stay inside the node's repo root.
  5. Writes a mechanical README.md run-card skeleton (branch/commit per
     node, job list, timings, exit codes) into the snapshot folder --
     content only, no interpretation. Per CLAUDE.md, the scientific read of
     a run's results is Vinay's call, not this tool's.

This tool never commits or pushes to git. Snapshotting/committing the
result folder remains an explicit, separate, user-approved step.

REVIEW FIXES (2026-08-24, applied same day as the first cut)
--------------------------------------------------------------
- CRITICAL: job["cwd"] was never checked before being folded into a
  mirror_discard delete target, so a manifest with cwd="../../.." could
  rm -rf outside the repo entirely (worked example in review: a shared
  node's whole home directory). Fixed by resolve_scoped_path(): every
  keep/mirror_discard entry is resolved to an absolute path via
  os.path.normpath and REJECTED unless it is a descendant of the node's
  (already-validated) repo root, regardless of which manifest field
  introduced the escape.
- The mirror "confirmed" byte-count check had a fallback that could
  silently measure a mirror destination's PARENT directory (picking up
  stale bytes from a previous run) instead of the specific path about to
  be deleted. Fixed: confirmation now requires the exact expected copy to
  exist under its own name; any mismatch fails closed (no discard).
- No sentinel-clearing before (re)launch meant a network-killed multi-hour
  job could be misreported as a clean exit=0 by reading a stale .exit file
  left over from an earlier run of the same job id. Fixed: launch always
  removes any prior log/.exit files first.
- Remote jobs ran attached to their launching SSH session, so a transient
  network drop (WiFi blip, laptop sleep, NAT idle-reap -- all realistic
  over an hours-long OGS run, and made more likely by the job's own
  output being redirected to a file, leaving the SSH channel looking idle)
  would SIGHUP-kill the in-progress simulation. Fixed: remote jobs launch
  via nohup with stdio fully closed/redirected away from the SSH channel,
  so the process survives the launching session closing; completion is
  then detected by polling for the sentinel file rather than blocking on
  that session, so a transient poll failure just delays detection instead
  of killing or misreporting the job. (Deliberately not using setsid --
  it's Linux-only and absent on macOS, which is every node in the
  roster; an earlier version of this fix depended on it and silently
  never launched anything as a result -- caught by smoke-testing after
  the fix, not by the review itself.)
- A node config missing "cores"/"ogs_repo" produced a raw KeyError deep in
  probe_all(); load_nodes() now validates required fields up front with a
  clear message.
- ogs_repo is now sanity-checked at resolution time (must not equal "/" or
  the node's own $HOME, must be at least a few path segments deep) so a
  nodes.json typo can't turn the whole discard-scope root into something
  dangerously broad.
- Found on the first real batch (2026-08-24, MS33 I/III/IV/VII): local
  jobs were collected via a sequential blocking loop (collect_local(h)
  per handle, in order), so a fast job ordered right after a slow job on
  the same node had its "end" timestamp captured only once the earlier
  job's .wait() returned -- the job itself finished on time, but its
  reported wall-clock time was wrong (inflated to match the slower job
  ahead of it in the list). Fixed: every handle (local and remote alike)
  is now collected via one thread each, same pattern already used for
  remote polling.

USAGE
-----
  python3 dispatch.py probe [--nodes nodes.json]
  python3 dispatch.py plan  --jobs jobs.json [--nodes nodes.json]
  python3 dispatch.py run   --jobs jobs.json --snapshot-name NAME
                             [--nodes nodes.json] [--mirror-root PATH]
                             [--snapshot-root PATH] [--discard] [--dry-run]
                             [--poll-interval SECONDS]

JOB MANIFEST (JSON array)
--------------------------
  [
    {
      "id": "unique-job-id",
      "cmd": "shell command, run via bash -lc inside cwd",
      "cwd": "optional, relative to the node's ogs_repo (default: repo root)",
      "keep": ["relative paths to always copy back; never deleted anywhere"],
      "mirror_discard": ["relative paths of bulk output; mirrored, then
                           deleted remotely only with --discard + verified"]
    }, ...
  ]

`cmd` runs with full shell privilege in `cwd` (default: the node's live
`ogs_repo` checkout -- NOT an isolated worktree). Treat a job manifest with
the same care as a command you'd type by hand: nothing sandboxes it from
running `git checkout`/`git reset --hard`/etc. against whatever branch is
currently active on that node. See CLAUDE.md's "work in a fresh worktree,
do not disturb active branches" rule -- that's the manifest author's job to
respect, same as it would be by hand; point `cwd` at a worktree if a job
needs isolation from the node's main checkout.
"""
import argparse
import json
import os
import re
import shlex
import subprocess
import sys
import threading
import time
from datetime import datetime, timezone
from pathlib import Path

DEFAULT_NODES_CONFIG = Path(__file__).parent / "nodes.json"
DEFAULT_MIRROR_ROOT = Path.home() / "ogs-models" / "_dispatch_mirror"
DEFAULT_SNAPSHOT_ROOT = Path.home() / "ogs-models" / "_dispatch_snapshots"
SSH_CONNECT_TIMEOUT = 5
DEFAULT_POLL_INTERVAL_S = 20
WARN_AFTER_CONSECUTIVE_UNREACHABLE_POLLS = 30  # ~10 min at the default interval

REQUIRED_NODE_FIELDS = ("type", "cores", "ogs_repo")


# ----------------------------------------------------------------------------
# node config + probing
# ----------------------------------------------------------------------------

def load_nodes(path):
    with open(path) as f:
        cfg = json.load(f)
    cfg.pop("_comment", None)
    for name, node in cfg.items():
        for field in REQUIRED_NODE_FIELDS:
            if field not in node:
                raise ValueError(f"node {name!r} in {path} is missing required field {field!r}")
        if node["type"] not in ("local", "remote"):
            raise ValueError(f"node {name!r}: type must be 'local' or 'remote', got {node['type']!r}")
        if node["type"] == "remote" and "ssh_alias" not in node:
            raise ValueError(f"node {name!r}: type 'remote' requires 'ssh_alias'")
        node.setdefault("respect_active_use", False)
        node.setdefault("idle_threshold_s", 120)
        node["name"] = name
    return cfg


def load_jobs(path):
    with open(path) as f:
        jobs = json.load(f)
    seen = set()
    for j in jobs:
        for req in ("id", "cmd"):
            if req not in j:
                raise ValueError(f"job missing required field {req!r}: {j}")
        if j["id"] in seen:
            raise ValueError(f"duplicate job id: {j['id']}")
        seen.add(j["id"])
        j.setdefault("cwd", "")
        j.setdefault("keep", [])
        j.setdefault("mirror_discard", [])
    return jobs


def ssh_run(alias, remote_cmd, timeout=SSH_CONNECT_TIMEOUT + 15):
    """Run one command on a remote node; remote_cmd is passed as a SINGLE
    argv element to ssh so no ssh-side space-joining/quoting ambiguity."""
    try:
        return subprocess.run(
            ["ssh", "-o", f"ConnectTimeout={SSH_CONNECT_TIMEOUT}", "-o", "BatchMode=yes",
             alias, remote_cmd],
            capture_output=True, text=True, timeout=timeout,
        )
    except subprocess.TimeoutExpired:
        return subprocess.CompletedProcess(args=[], returncode=124, stdout="", stderr="ssh timeout")


def probe_reachable(node):
    if node["type"] == "local":
        return True
    r = ssh_run(node["ssh_alias"], "true")
    return r.returncode == 0


def probe_idle_seconds(node):
    """HID idle time in seconds, or None if unavailable/local."""
    if node["type"] == "local":
        cmd = "ioreg -c IOHIDSystem | grep -o 'HIDIdleTime\" = [0-9]*' | head -1"
        out = subprocess.run(["bash", "-lc", cmd], capture_output=True, text=True)
    else:
        out = ssh_run(node["ssh_alias"], "ioreg -c IOHIDSystem | grep -o 'HIDIdleTime\" = [0-9]*' | head -1")
    m = re.search(r"(\d+)", out.stdout or "")
    if not m:
        return None
    return int(m.group(1)) / 1e9


def probe_all(nodes):
    """Returns {name: {reachable, idle_s, busy, available, reason, cores}}."""
    status = {}
    for name, node in nodes.items():
        reachable = probe_reachable(node)
        idle_s = None
        busy = False
        reason = ""
        if not reachable:
            reason = "unreachable"
        elif node.get("respect_active_use"):
            idle_s = probe_idle_seconds(node)
            if idle_s is None:
                reason = "idle-check failed; treating as busy (fail-safe)"
                busy = True
            elif idle_s < node["idle_threshold_s"]:
                busy = True
                reason = f"active use (idle {idle_s:.0f}s < threshold {node['idle_threshold_s']}s)"
        available = reachable and not busy
        status[name] = dict(reachable=reachable, idle_s=idle_s, busy=busy,
                             available=available, reason=reason, cores=node["cores"])
    return status


def print_probe_table(status):
    print(f"{'node':<14} {'reachable':<10} {'busy':<6} {'cores':<6} verdict")
    for name, s in status.items():
        verdict = "AVAILABLE" if s["available"] else f"EXCLUDED ({s['reason']})"
        print(f"{name:<14} {str(s['reachable']):<10} {str(s['busy']):<6} {s['cores']:<6} {verdict}")


# ----------------------------------------------------------------------------
# atomization
# ----------------------------------------------------------------------------

def plan_rounds(available_node_names, nodes, jobs):
    """Greedy least-loaded-by-core-ratio assignment. Returns list of rounds;
    each round is {node_name: [job, ...]}, jobs assigned at most `cores`
    per node per round; leftovers spill into subsequent rounds."""
    pending = list(jobs)
    rounds = []
    while pending:
        round_assignment = {n: [] for n in available_node_names}
        cores = {n: nodes[n]["cores"] for n in available_node_names}
        while pending:
            candidates = [n for n in available_node_names if len(round_assignment[n]) < cores[n]]
            if not candidates:
                break
            candidates.sort(key=lambda n: (len(round_assignment[n]) / cores[n], n))
            chosen = candidates[0]
            round_assignment[chosen].append(pending.pop(0))
        round_assignment = {n: js for n, js in round_assignment.items() if js}
        if not round_assignment:
            raise RuntimeError("no available node could take any pending job "
                                "(all excluded, or a node has cores <= 0?)")
        rounds.append(round_assignment)
    return rounds


def print_plan(rounds):
    for i, round_assignment in enumerate(rounds, 1):
        print(f"\n-- round {i} --")
        for node, js in round_assignment.items():
            print(f"  {node}: {', '.join(j['id'] for j in js)}")


# ----------------------------------------------------------------------------
# path resolution -- the one place every remote/local path gets built
# ----------------------------------------------------------------------------

def resolve_node_homes(nodes, available_names):
    """Resolve each available node's real $HOME, rewrite a leading '~' in
    ogs_repo to it, and sanity-check ogs_repo isn't dangerously broad.

    The $HOME resolution matters because shlex.quote()-ing a '~/...' path
    disables shell tilde expansion (quoting turns off all expansion), so
    paths must be absolute before they're ever quoted -- this was the
    original bug (cd '~/git/ogs' fails, quoted tilde never expands)."""
    for name in available_names:
        node = nodes[name]
        if node["type"] == "local":
            node["home"] = str(Path.home())
        else:
            r = ssh_run(node["ssh_alias"], "echo $HOME")
            home = r.stdout.strip()
            if not home:
                raise RuntimeError(f"could not resolve $HOME on node {name}: {r.stderr}")
            node["home"] = home
        if node["ogs_repo"].startswith("~"):
            node["ogs_repo"] = node["home"] + node["ogs_repo"][1:]
        node["ogs_repo"] = os.path.normpath(node["ogs_repo"])
        if node["ogs_repo"] in ("/", node["home"]) or len(Path(node["ogs_repo"]).parts) < 3:
            raise RuntimeError(
                f"node {name!r} ogs_repo={node['ogs_repo']!r} looks too broad (equals '/' or "
                f"$HOME, or is too shallow) -- refusing to use it as a discard-scope root. "
                f"Fix nodes.json.")


def resolve_scoped_path(node, cwd, relpath):
    """Resolve (cwd, relpath) against node['ogs_repo'] and require the
    result to stay inside the repo root. Returns the absolute path, or
    raises ValueError if it would escape -- this is the containment check
    that closes the cwd="../../.." traversal the 2026-08-24 review found:
    it catches an escape regardless of whether cwd or relpath introduced
    it, because it checks the final resolved path, not the raw strings."""
    repo_root = node["ogs_repo"]  # already absolute + validated by resolve_node_homes
    joined = os.path.normpath(os.path.join(repo_root, cwd, relpath))
    if joined != repo_root and not joined.startswith(repo_root + os.sep):
        raise ValueError(f"{relpath!r} under cwd {cwd!r} resolves to {joined!r}, "
                          f"outside repo root {repo_root!r}")
    return joined


def remote_repo_path(node, sub=""):
    """Used only for job execution cwd (launch_job) and git_state -- NOT for
    anything that reaches remote_rm(); those go through resolve_scoped_path
    instead. Running a job wherever its manifest points is the same trust
    level as the job's own `cmd` (operator-authored); deleting things is not,
    which is why it gets the stricter containment check."""
    base = node["ogs_repo"]
    return f"{base}/{sub}".rstrip("/") if sub else base


def git_state(node):
    cmd = f"cd {shlex.quote(remote_repo_path(node))} && git branch --show-current && git rev-parse --short HEAD"
    if node["type"] == "local":
        out = subprocess.run(["bash", "-lc", cmd], capture_output=True, text=True)
    else:
        out = ssh_run(node["ssh_alias"], cmd)
    lines = (out.stdout or "").strip().splitlines()
    if len(lines) >= 2:
        return {"branch": lines[0], "commit": lines[1]}
    return {"branch": "?", "commit": "?"}


# ----------------------------------------------------------------------------
# execution -- local jobs use ordinary Popen; remote jobs launch detached
# (nohup+setsid, stdio closed) and are tracked by polling a sentinel file,
# so a dropped SSH session can neither kill the job nor be misread as one
# completing (see REVIEW FIXES in the module docstring).
# ----------------------------------------------------------------------------

def launch_local(node, job, log_dir):
    log_dir.mkdir(parents=True, exist_ok=True)
    local_log = log_dir / f"{job['id']}.log"
    exit_file = Path(str(local_log) + ".exit")
    outer_log = log_dir / f"{job['id']}.outer.log"
    cwd = remote_repo_path(node, job["cwd"])
    full_cmd = (f"rm -f {shlex.quote(str(local_log))} {shlex.quote(str(exit_file))} && "
                f"cd {shlex.quote(cwd)} && ( {job['cmd']} ) > {shlex.quote(str(local_log))} 2>&1 ; "
                f"echo $? > {shlex.quote(str(exit_file))}")
    outer_fh = open(outer_log, "wb")
    proc = subprocess.Popen(["bash", "-lc", full_cmd], stdout=outer_fh, stderr=subprocess.STDOUT)
    return dict(job=job, node=node, mode="popen", proc=proc, outer_fh=outer_fh,
                local_log=local_log, exit_file=exit_file, outer_log=outer_log, start=time.time())


def launch_remote(node, job, log_dir):
    log_dir.mkdir(parents=True, exist_ok=True)
    local_log = log_dir / f"{job['id']}.log"
    outer_log = log_dir / f"{job['id']}.launch.log"
    remote_dispatch_dir = f"{node['home']}/.dispatch_logs"
    remote_log = f"{remote_dispatch_dir}/{job['id']}.log"
    remote_exit = remote_log + ".exit"
    cwd = remote_repo_path(node, job["cwd"])
    inner = (f"cd {shlex.quote(cwd)} && ( {job['cmd']} ) > {shlex.quote(remote_log)} 2>&1 ; "
             f"echo $? > {shlex.quote(remote_exit)}")
    # nohup + fully-closed stdio detaches the job from this SSH session so
    # it survives the session closing -- a transient network drop can no
    # longer SIGHUP-kill an hours-long OGS run. Deliberately NOT using
    # setsid: it's Linux-only (util-linux) and absent on macOS, which is
    # every node in the roster (verified 2026-08-24 on both the mac mini
    # and this machine) -- nohup + closed stdio + backgrounding is the
    # standard macOS-portable idiom for this and needs no setsid. This
    # launch command itself returns almost immediately (it just
    # backgrounds the job and exits).
    launch_cmd = (f"mkdir -p {shlex.quote(remote_dispatch_dir)} && "
                  f"rm -f {shlex.quote(remote_log)} {shlex.quote(remote_exit)} && "
                  f"nohup bash -c {shlex.quote(inner)} </dev/null >/dev/null 2>&1 &")
    r = subprocess.run(
        ["ssh", "-o", f"ConnectTimeout={SSH_CONNECT_TIMEOUT}",
         "-o", "ServerAliveInterval=30", "-o", "ServerAliveCountMax=6",
         node["ssh_alias"], launch_cmd],
        capture_output=True, text=True, timeout=SSH_CONNECT_TIMEOUT + 15,
    )
    Path(outer_log).write_text((r.stdout or "") + (r.stderr or ""))
    return dict(job=job, node=node, mode="poll", remote_log=remote_log, remote_exit=remote_exit,
                local_log=local_log, outer_log=outer_log, launch_rc=r.returncode, start=time.time())


def wait_for_remote(handle, poll_interval):
    node = handle["node"]
    consecutive_unreachable = 0
    while True:
        r = ssh_run(node["ssh_alias"], f"cat {shlex.quote(handle['remote_exit'])} 2>/dev/null")
        content = r.stdout.strip()
        if content:
            handle["exit_code"] = content
            break
        if r.returncode == 255:  # ssh itself couldn't connect, vs. "file not there yet"
            consecutive_unreachable += 1
            if consecutive_unreachable == WARN_AFTER_CONSECUTIVE_UNREACHABLE_POLLS:
                print(f"    warning: {node['name']} unreachable for "
                      f"~{consecutive_unreachable * poll_interval}s while polling "
                      f"{handle['job']['id']} -- job may still be running remotely, "
                      f"still waiting", file=sys.stderr)
        else:
            consecutive_unreachable = 0
        time.sleep(poll_interval)
    handle["end"] = time.time()
    subprocess.run(["scp", "-q", f"{node['ssh_alias']}:{handle['remote_log']}", str(handle["local_log"])],
                    capture_output=True)


def collect_local(handle):
    handle["proc"].wait()
    handle["outer_fh"].close()
    handle["end"] = time.time()
    code = handle["exit_file"].read_text().strip() if handle["exit_file"].exists() else ""
    handle["exit_code"] = code or "?"
    if handle["exit_code"] == "?":
        print(f"    warning: could not read exit code for {handle['job']['id']} on "
              f"{handle['node']['name']}; see {handle['outer_log']}", file=sys.stderr)


def run_round(round_assignment, nodes, log_dir, poll_interval):
    handles = []
    for node_name, js in round_assignment.items():
        node = nodes[node_name]
        for job in js:
            if node["type"] == "local":
                handles.append(launch_local(node, job, log_dir))
            else:
                h = launch_remote(node, job, log_dir)
                handles.append(h)
                if h["launch_rc"] != 0:
                    print(f"  warning: launch failed for {node_name}/{job['id']} "
                          f"(rc={h['launch_rc']}); see {h['outer_log']}", file=sys.stderr)

    # Collect ALL handles (local and remote) concurrently via threads, one
    # per job. A prior version collected local jobs in a sequential
    # blocking loop (collect_local(h) per h in order) -- correctness was
    # fine (each Popen already ran concurrently since launch), but a fast
    # job ordered after a slow job on the same node had its "end" timestamp
    # captured only once the earlier job's .wait() returned, misreporting
    # its wall-clock time as inflated by however long that earlier job
    # took. Found 2026-08-24 on a real batch: modelVII_freeswelling (truly
    # 211.8s per its own log) was reported as 1264.4s, matching
    # modelIV_pellets exactly, because it was collected right after it.
    threads = [threading.Thread(target=wait_for_remote, args=(h, poll_interval), daemon=True)
               for h in handles if h["mode"] == "poll"]
    threads += [threading.Thread(target=collect_local, args=(h,), daemon=True)
                for h in handles if h["mode"] == "popen"]
    for t in threads:
        t.start()
    for t in threads:
        t.join()

    for h in handles:
        print(f"  {h['node']['name']}/{h['job']['id']}: exit={h.get('exit_code', '?')} "
              f"wall={h.get('end', h['start']) - h['start']:.1f}s")
    return handles


# ----------------------------------------------------------------------------
# assimilation: keep (always) + mirror_discard (mirror, then optional discard)
# ----------------------------------------------------------------------------

def rsync_pull(node, abs_src, local_dest_dir, checksum=False):
    """abs_src is an already-resolved absolute path on `node` (local: a
    plain path; remote: prefixed alias:path is added here)."""
    local_dest_dir.mkdir(parents=True, exist_ok=True)
    src = abs_src if node["type"] == "local" else f"{node['ssh_alias']}:{abs_src}"
    args = ["rsync", "-a"]
    if checksum:
        args.append("--checksum")
    args += [src, str(local_dest_dir) + "/"]
    return subprocess.run(args, capture_output=True, text=True)


def remote_du_kb(node, abs_path):
    if node["type"] == "local":
        out = subprocess.run(["bash", "-lc", f"du -sk {shlex.quote(abs_path)} 2>/dev/null | cut -f1"],
                              capture_output=True, text=True)
    else:
        out = ssh_run(node["ssh_alias"], f"du -sk {shlex.quote(abs_path)} 2>/dev/null | cut -f1")
    try:
        return int((out.stdout or "0").strip().splitlines()[0])
    except (ValueError, IndexError):
        return None


def local_du_kb(path):
    if not path.exists():
        return None
    out = subprocess.run(["du", "-sk", str(path)], capture_output=True, text=True)
    try:
        return int(out.stdout.split()[0])
    except (ValueError, IndexError):
        return None


UNSAFE_DISCARD_PATHS = {"", ".", "..", "/", "~"}


def is_safe_discard_path(relpath):
    """Cheap syntactic first-line filter (fails fast, readable message).
    NOT the authoritative check -- resolve_scoped_path() is, because this
    function only ever looks at `relpath` in isolation and can't see a
    traversal introduced via `cwd`."""
    if relpath in UNSAFE_DISCARD_PATHS:
        return False
    if relpath.startswith("/") or relpath.startswith("~") or ".." in Path(relpath).parts:
        return False
    return True


def remote_rm(node, abs_path):
    cmd = f"rm -rf {shlex.quote(abs_path)}"
    if node["type"] == "local":
        return subprocess.run(["bash", "-lc", cmd], capture_output=True, text=True)
    return ssh_run(node["ssh_alias"], cmd)


def assimilate(handle, snapshot_dir, mirror_root, do_discard, dry_run, log_lines):
    job, node = handle["job"], handle["node"]
    node_snap = snapshot_dir / node["name"] / job["id"]
    node_mirror = mirror_root / snapshot_dir.name / node["name"] / job["id"]

    for relpath in job["keep"]:
        try:
            abs_src = resolve_scoped_path(node, job["cwd"], relpath)
        except ValueError as e:
            log_lines.append(f"SKIP keep for {node['name']}:{relpath} -- {e}")
            continue
        dest = node_snap / Path(relpath).parent
        if dry_run:
            log_lines.append(f"[dry-run] would rsync keep {node['name']}:{abs_src} -> {dest}")
            continue
        r = rsync_pull(node, abs_src, dest)
        log_lines.append(f"keep  {node['name']}:{relpath} -> {dest}  rsync_rc={r.returncode}")

    for relpath in job["mirror_discard"]:
        if not is_safe_discard_path(relpath):
            log_lines.append(f"SKIP mirror_discard for {node['name']}:{relpath} -- failed safety guard")
            continue
        try:
            abs_src = resolve_scoped_path(node, job["cwd"], relpath)
        except ValueError as e:
            log_lines.append(f"SKIP mirror_discard for {node['name']}:{relpath} -- {e}")
            continue
        if abs_src == node["ogs_repo"]:
            log_lines.append(f"SKIP mirror_discard for {node['name']}:{relpath} -- "
                              f"resolves to the repo root itself, refusing")
            continue
        dest = node_mirror / Path(relpath).parent
        if dry_run:
            log_lines.append(f"[dry-run] would mirror {node['name']}:{abs_src} -> {dest}"
                              f"{' and discard remote copy' if do_discard else ''}")
            continue
        r = rsync_pull(node, abs_src, dest, checksum=True)
        confirmed = False
        if r.returncode == 0:
            remote_kb = remote_du_kb(node, abs_src)
            copied = dest / Path(abs_src).name
            local_kb = local_du_kb(copied) if copied.exists() else None
            if remote_kb is not None and local_kb is not None and remote_kb > 0:
                confirmed = abs(remote_kb - local_kb) <= max(4, int(0.01 * remote_kb))
        log_lines.append(f"mirror {node['name']}:{relpath} -> {dest}  rsync_rc={r.returncode} confirmed={confirmed}")
        if do_discard and node["type"] == "remote":
            if confirmed:
                rr = remote_rm(node, abs_src)
                log_lines.append(f"discard {node['name']}:{relpath}  rm_rc={rr.returncode}")
            else:
                log_lines.append(f"KEEP {node['name']}:{relpath} on node -- mirror not confirmed, refusing to discard")
        elif do_discard and node["type"] == "local":
            log_lines.append(f"KEEP {node['name']}:{relpath} -- discard only applies to remote nodes")


# ----------------------------------------------------------------------------
# README run-card
# ----------------------------------------------------------------------------

def write_readme(snapshot_dir, nodes, node_git, all_handles, args):
    lines = []
    lines.append(f"# Dispatch run: {snapshot_dir.name}")
    lines.append("")
    lines.append(f"Generated: {datetime.now(timezone.utc).isoformat()}Z by tools/dispatch_nodes/dispatch.py")
    lines.append("")
    lines.append("Mechanical record only -- interpretation, expected values, and")
    lines.append("pass/fail scientific judgement are Vinay's call (CLAUDE.md).")
    lines.append("")
    lines.append("## Nodes used")
    lines.append("")
    lines.append("| node | branch | commit |")
    lines.append("|---|---|---|")
    for name, g in node_git.items():
        lines.append(f"| {name} | {g['branch']} | {g['commit']} |")
    lines.append("")
    lines.append("## Jobs")
    lines.append("")
    lines.append("| id | node | exit | wall (s) | cmd |")
    lines.append("|---|---|---|---|---|")
    for h in all_handles:
        wall = h.get("end", h["start"]) - h["start"]
        cmd_display = h["job"]["cmd"].replace("|", "\\|")
        lines.append(f"| {h['job']['id']} | {h['node']['name']} | {h.get('exit_code', '?')} | "
                      f"{wall:.1f} | `{cmd_display}` |")
    lines.append("")
    lines.append("## Open items")
    lines.append("")
    lines.append("- TODO(Vinay): read the logs under `<node>/<job-id>/`, confirm results.")
    (snapshot_dir / "README.md").write_text("\n".join(lines) + "\n")


# ----------------------------------------------------------------------------
# commands
# ----------------------------------------------------------------------------

def cmd_probe(args):
    nodes = load_nodes(args.nodes)
    status = probe_all(nodes)
    print_probe_table(status)


def cmd_plan(args):
    nodes = load_nodes(args.nodes)
    jobs = load_jobs(args.jobs)
    status = probe_all(nodes)
    print_probe_table(status)
    available = [n for n, s in status.items() if s["available"]]
    if not available:
        print("\nno available nodes -- nothing to plan", file=sys.stderr)
        sys.exit(1)
    resolve_node_homes(nodes, available)  # keeps `plan` from silently drifting out of sync with `run`
    rounds = plan_rounds(available, nodes, jobs)
    print_plan(rounds)


def cmd_run(args):
    nodes = load_nodes(args.nodes)
    jobs = load_jobs(args.jobs)
    status = probe_all(nodes)
    print_probe_table(status)
    available = [n for n, s in status.items() if s["available"]]
    if not available:
        print("\nno available nodes -- ask Vinay to bring one online, not proceeding", file=sys.stderr)
        sys.exit(1)

    resolve_node_homes(nodes, available)
    rounds = plan_rounds(available, nodes, jobs)
    print_plan(rounds)

    mirror_root = Path(args.mirror_root).expanduser()
    snapshot_root = Path(args.snapshot_root).expanduser()
    snapshot_dir = snapshot_root / args.snapshot_name
    log_dir = snapshot_dir / "_dispatch_run_logs"

    if args.dry_run:
        print("\n[dry-run] stopping before any execution/mutation")
        return

    node_git = {n: git_state(nodes[n]) for n in available}
    all_handles = []
    log_lines = []

    for i, round_assignment in enumerate(rounds, 1):
        print(f"\n== executing round {i} ==")
        all_handles.extend(run_round(round_assignment, nodes, log_dir, args.poll_interval))

    print("\n== assimilating ==")
    for h in all_handles:
        assimilate(h, snapshot_dir, mirror_root, args.discard, False, log_lines)
    for line in log_lines:
        print(" ", line)

    write_readme(snapshot_dir, nodes, node_git, all_handles, args)
    print(f"\nsnapshot: {snapshot_dir}")
    print(f"mirror:   {mirror_root / snapshot_dir.name}")
    print("Nothing was committed or pushed -- that stays a separate, explicit step.")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="command", required=True)

    pp = sub.add_parser("probe", help="check reachability/idle-status of all configured nodes")
    pp.add_argument("--nodes", default=str(DEFAULT_NODES_CONFIG))
    pp.set_defaults(func=cmd_probe)

    pl = sub.add_parser("plan", help="preview the atomized job split, no execution")
    pl.add_argument("--jobs", required=True)
    pl.add_argument("--nodes", default=str(DEFAULT_NODES_CONFIG))
    pl.set_defaults(func=cmd_plan)

    pr = sub.add_parser("run", help="probe, atomize, dispatch, assimilate")
    pr.add_argument("--jobs", required=True)
    pr.add_argument("--snapshot-name", required=True)
    pr.add_argument("--nodes", default=str(DEFAULT_NODES_CONFIG))
    pr.add_argument("--mirror-root", default=str(DEFAULT_MIRROR_ROOT))
    pr.add_argument("--snapshot-root", default=str(DEFAULT_SNAPSHOT_ROOT))
    pr.add_argument("--discard", action="store_true",
                     help="after a verified mirror, delete mirror_discard paths on remote nodes")
    pr.add_argument("--dry-run", action="store_true",
                     help="probe + print the plan only, no execution or mutation")
    pr.add_argument("--poll-interval", type=float, default=DEFAULT_POLL_INTERVAL_S,
                     help=f"seconds between remote-job completion checks (default {DEFAULT_POLL_INTERVAL_S}); "
                          f"does not affect how long a job itself takes")
    pr.set_defaults(func=cmd_run)

    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
