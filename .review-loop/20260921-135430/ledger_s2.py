#!/usr/bin/env python3
"""Ledger helper for run 20260921-135430, session 2 (mac mini).

The original `review-fix-loop` skill (and its scripts/ledger.py) exists only on
the laptop and was never pushed, so it is absent on this host.  This script is a
reconstruction that appends to the SAME append-only ledger.jsonl with the SAME
event schema observed in rounds 1-2, and regenerates report.md in the same
layout.  It deliberately has no `init` subcommand: this run must never be
re-initialised.

Usage:
  ledger_s2.py --run <dir> append --json '<one event object>'
  ledger_s2.py --run <dir> report
  ledger_s2.py --run <dir> answer --qid q1 --text "..."
"""
import argparse
import datetime
import json
import os
import sys


def load(run):
    path = os.path.join(run, "ledger.jsonl")
    with open(path) as fh:
        return [json.loads(line) for line in fh if line.strip()]


def append(run, obj):
    obj.setdefault("ts", datetime.datetime.now().replace(microsecond=0).isoformat())
    with open(os.path.join(run, "ledger.jsonl"), "a") as fh:
        fh.write(json.dumps(obj, sort_keys=True) + "\n")
    return obj


def report(run):
    rows = load(run)
    start = next(r for r in rows if r["event"] == "run_start")
    findings = {}
    for r in rows:
        if r["event"] == "findings_added":
            for f in r["findings"]:
                findings[f["id"]] = f
    decisions = {}
    for r in rows:
        if r["event"] == "decision":
            decisions[r["id"]] = r
    resolved = {}
    for r in rows:
        if r["event"] == "resolved":
            resolved[r["id"]] = r
    questions = [r for r in rows if r["event"] == "question"]
    answers = {r["qid"]: r for r in rows if r["event"] == "answer"}
    round_ends = [r for r in rows if r["event"] == "round_end"]
    verifies = [r for r in rows if r["event"] == "verify"]
    notes = [r for r in rows if r["event"] == "note"]

    verify_cmd = verifies[-1]["cmd"] if verifies else ""
    stop = next((r for r in rows if r["event"] == "run_end"), None)

    out = []
    out.append("# Review–fix loop — %s" % start["run_id"])
    out.append("")
    out.append("- repo: `%s`" % start["repo"])
    out.append("- base: `%s`" % start["base"])
    out.append("- verify: `%s`" % verify_cmd)
    out.append("- rounds: %d (budget %d)" % (len(round_ends), start["max_rounds"]))
    if stop:
        out.append("- stop: **%s** — %s" % (stop.get("verdict", "stopped"), stop.get("reason", "")))
    else:
        out.append("- stop: **in progress** — loop still running")
    out.append("")
    out.append("## Rounds")
    out.append("")
    for re_ in round_ends:
        out.append("- round %s: new=%s fixed=%s open=%s verify=%s verdict=%s"
                   % (re_["round"], re_["new"], re_["fixed"], re_["open"],
                      re_["verify"], re_["verdict"]))
        for v in verifies:
            if v["round"] == re_["round"]:
                out.append("  - verify %s: `%s` (%s)" % (v["result"], v["cmd"], v["detail"]))
    out.append("")

    def block(title, ids):
        out.append("## %s (%d)" % (title, len(ids)))
        out.append("")
        for fid in ids:
            f = findings.get(fid, {})
            d = decisions.get(fid, {})
            out.append("- `%s` [%s/%s] %s:%s — %s"
                       % (fid, f.get("severity", "?"), f.get("axis", "?"),
                          f.get("file", "?"), f.get("line", "?"), f.get("summary", "")))
            if d.get("reason"):
                out.append("  - reason: %s" % d["reason"])
            r = resolved.get(fid)
            if r:
                if r.get("commit"):
                    out.append("  - commit: `%s` (fixup of `%s`)" % (r["commit"], (r.get("target") or "")[:12]))
                elif r.get("no_commit_reason"):
                    out.append("  - no commit: %s" % r["no_commit_reason"])
        out.append("")

    sev_rank = {"major": 0, "minor": 1}

    def key(fid):
        f = findings.get(fid, {})
        try:
            line = int(f.get("line") or 0)
        except (TypeError, ValueError):
            line = 0
        return (sev_rank.get(f.get("severity"), 2), f.get("file", ""), line)

    def sel(kind):
        return sorted([i for i, d in decisions.items() if d["decision"] == kind], key=key)
    fixed = sel("fix")
    deferred = sel("defer")
    rejected = sel("reject")
    block("Fixed", fixed)
    block("Deferred", deferred)
    block("Rejected", rejected)

    unanswered = [q for q in questions if q["qid"] not in answers]
    out.append("## Questions for you (%d unanswered)" % len(unanswered))
    out.append("")
    for q in questions:
        a = answers.get(q["qid"])
        state = a["text"] if a else "**unanswered**"
        out.append("- `%s` round %s (%s, %s) \u2014 %s"
                   % (q["qid"], q["round"], q.get("source", "loop"), state, q["text"]))
        if q.get("assumption"):
            out.append("  - the loop proceeded assuming: %s" % q["assumption"])
        out.append("  - finding: `%s`" % q["id"])
    out.append("")
    out.append("## Notes")
    out.append("")
    for n in notes:
        out.append("- round %s: %s" % (n["round"], n["text"]))
    out.append("")

    text = "\n".join(out)
    with open(os.path.join(run, "report.md"), "w") as fh:
        fh.write(text)
    return text


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--run", required=True)
    sub = p.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("append")
    a.add_argument("--json", required=True)
    sub.add_parser("report")
    ans = sub.add_parser("answer")
    ans.add_argument("--qid", required=True)
    ans.add_argument("--text", required=True)
    args = p.parse_args()

    if args.cmd == "append":
        append(args.run, json.loads(args.json))
    elif args.cmd == "answer":
        append(args.run, {"event": "answer", "qid": args.qid, "text": args.text})
        report(args.run)
    else:
        sys.stdout.write(report(args.run))


if __name__ == "__main__":
    main()
