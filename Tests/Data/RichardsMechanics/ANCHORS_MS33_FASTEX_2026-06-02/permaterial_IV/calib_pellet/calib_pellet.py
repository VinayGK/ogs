#!/usr/bin/env python3
"""Pellet K calibration probe for per-material Model IV (FASTEX 2026-06-02).
Substitutes a trial vdw_augmentation_prefactor K into the single-element pellet
constant-volume probe (pellet phi0/n_s/n_l0 + FASTEX alpha=3e-13, k_rel=Constant 0.1,
min_dt=1e-4), runs OGS (OMP_NUM_THREADS=1), and reports the steady-state mean stress
sigma_sat vs the Dixon (2023) dd0.9 median target 0.350 MPa.
Usage: python3 calib_pellet.py <K1> [K2] [K3] ...
"""
import os, re, shutil, subprocess, glob, sys, math
import numpy as np
import meshio

HERE = os.path.dirname(os.path.abspath(__file__))
OGS = "/Users/vinaykumar/git/build/dsm_swelling_full_pdisj-release/bin/ogs"
TEMPLATE = HERE + "/pellet_probe_template.prj"
TARGET = 0.003 * math.exp(5.2883 * 0.9)   # Dixon (2023) MX-80 Fig.1 median, dd0.9 -> 0.350 MPa


def meanp(s):
    s = np.atleast_2d(s)
    return -s[:, :3].sum(axis=1) / 3.0 / 1e6


def run_K(K):
    tmpl = open(TEMPLATE, encoding="latin-1").read()
    txt = tmpl.replace("__K__", repr(float(K)))
    rundir = HERE + ("/run_K%g" % K)
    shutil.rmtree(rundir, ignore_errors=True)
    os.makedirs(rundir)
    for m in ["square_1e-2_quad_1e0.vtu", "square_1e-2_quad_1e0_left.vtu",
              "square_1e-2_quad_1e0_right.vtu", "square_1e-2_quad_1e0_top.vtu",
              "square_1e-2_quad_1e0_bottom.vtu"]:
        shutil.copy(HERE + "/" + m, rundir + "/" + m)
    prj = rundir + "/pellet_probe.prj"
    open(prj, "w", encoding="latin-1").write(txt)
    out = rundir + "/out"
    os.makedirs(out)
    env = dict(os.environ, OMP_NUM_THREADS="1")
    try:
        r = subprocess.run([OGS, "-o", out, prj], cwd=rundir,
                           capture_output=True, text=True, timeout=900, env=env)
        log = r.stdout + r.stderr
    except subprocess.TimeoutExpired as e:
        so = e.stdout.decode() if isinstance(e.stdout, bytes) else (e.stdout or "")
        se = e.stderr.decode() if isinstance(e.stderr, bytes) else (e.stderr or "")
        log = so + se + "\n[TIMEOUT]"
    done = "Simulation completed" in log
    fs = sorted(glob.glob(out + "/*_ts_*.vtu"),
                key=lambda f: int(re.search(r"_ts_(\d+)_", f).group(1)))
    fail = re.search(r"failed in time step #\d+ at t = ([0-9.]+)", log)
    tail = f" CRASH@{float(fail.group(1))/86400:.2f}d" if fail else ""
    if not fs:
        print(f"K={K:g}  NO OUTPUT{tail}")
        print("\n".join(log.splitlines()[-8:]))
        return None
    last = fs[-1]
    pd = meshio.read(last).point_data
    t_end = float(re.search(r"_t_([0-9.]+)\.vtu", last).group(1)) / 86400
    sig = meanp(pd["sigma"]).mean()
    mp = float(np.mean(np.ravel(pd["micro_pressure"]))) / 1e6 if "micro_pressure" in pd else float("nan")
    S = float(np.mean(np.ravel(pd["saturation"]))) if "saturation" in pd else float("nan")
    err = (sig - TARGET) / TARGET * 100
    print(f"K={K:g}  done={done} t_end={t_end:.0f}d  sigma_sat={sig:.4f} MPa "
          f"(target {TARGET:.4f}, err {err:+.1f}%)  micro_p={mp:.3f} MPa  S={S:.3f}{tail}")
    return sig


if __name__ == "__main__":
    print(f"# Dixon dd0.9 median target sigma_sat = {TARGET:.5f} MPa")
    for a in sys.argv[1:]:
        run_K(float(a))
