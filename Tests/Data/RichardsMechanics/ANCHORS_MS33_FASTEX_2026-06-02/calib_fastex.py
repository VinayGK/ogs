#!/usr/bin/env python3
"""Calibration probe: run ModelI dd1400 (binding case) with a trial micro-macro
mass_exchange_coefficient (alpha_M), report time-to-95%-steady-state and the
steady-state mean stress vs the Dixon EMDD=rho_d target (4.922 MPa @ dd1400).
Usage: python3 calib_fastex.py <alpha_M>
"""
import os, re, shutil, subprocess, glob, sys
import numpy as np
import meshio

RM = "/Users/vinaykumar/git/ogs-worktrees/dsm_native_hierarchical_wt/Tests/Data/RichardsMechanics"
OGS = "/Users/vinaykumar/git/build/dsm_swelling_full_pdisj-release/bin/ogs"
NEW = RM + "/ANCHORS_MS33_FASTEX_2026-06-02"
calib = NEW + "/calib"
alpha = float(sys.argv[1]) if len(sys.argv) > 1 else 1e-12

os.makedirs(calib, exist_ok=True)
src = RM + "/ANCHORS_MS33_ModelI/ms33_modelI_dd1400.prj"
txt = open(src, encoding="latin-1").read()
srcdir = os.path.dirname(src)
for m in re.findall(r"<mesh[^>]*>([^<]+)</mesh>", txt):
    if not m.endswith(".vtu"):
        continue   # skip <output> mesh-names (referenced by basename, path rewrite leaves them valid)
    shutil.copy(os.path.normpath(os.path.join(srcdir, m)), calib + "/" + os.path.basename(m))
    txt = txt.replace(">" + m + "<", ">" + os.path.basename(m) + "<")
txt = re.sub(r"(<mass_exchange_coefficient>)[^<]+(</mass_exchange_coefficient>)",
             rf"\g<1>{alpha:g}\g<2>", txt)
prj = calib + "/ms33_modelI_dd1400.prj"
open(prj, "w", encoding="latin-1").write(txt)
out = calib + ("/out_a%g" % alpha)
shutil.rmtree(out, ignore_errors=True); os.makedirs(out)

try:
    r = subprocess.run([OGS, "-o", out, prj], cwd=calib, capture_output=True, text=True, timeout=555)
    log = r.stdout + r.stderr
except subprocess.TimeoutExpired as e:
    log = (e.stdout or "") + (e.stderr or "") + "\n[TIMEOUT]"

done = "Simulation completed" in log
fs = sorted(glob.glob(out + "/*_ts_*.vtu"), key=lambda f: int(re.search(r"_ts_(\d+)_", f).group(1)))
if not fs:
    print(f"alpha={alpha:g} NO OUTPUT. log tail:")
    print("\n".join(log.splitlines()[-10:])); sys.exit()
rows = []
for f in fs:
    pd = meshio.read(f).point_data
    p = -np.atleast_2d(pd["sigma"])[:, :3].sum(1).mean() / 3 / 1e6
    t = float(re.search(r"_t_([0-9.]+)\.vtu", f).group(1)) / 86400
    rows.append((t, p))
pend = rows[-1][1]
t95 = next((t for t, p in rows if p >= 0.95 * pend), None)
fail = re.search(r"failed in time step #\d+ at t = ([0-9.]+)", log)
tail = f" CRASH@{float(fail.group(1))/86400:.1f}d" if fail else ""
print(f"alpha={alpha:g} done={done} p_ss={pend:.3f} MPa (Dixon 4.922) t95={t95} d{tail}")
print("  traj: " + ", ".join(f"{t:.1f}d:{p:.2f}" for t, p in rows))
