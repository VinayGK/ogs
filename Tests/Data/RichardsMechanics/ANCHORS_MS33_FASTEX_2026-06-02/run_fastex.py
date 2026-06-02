#!/usr/bin/env python3
"""General FASTEX model runner: prepare a model into the FASTEX folder with a given
alpha_exchange (and optional K override), run OGS, report centre mean-stress trajectory,
time-to-95%, flatness over the last 25% of the window, and completion/crash.
Usage: python3 run_fastex.py <suite LE|MCC> <key> <alpha> [K]
"""
import os, re, shutil, subprocess, glob, sys
import numpy as np
import meshio

RM = "/Users/vinaykumar/git/ogs-worktrees/dsm_native_hierarchical_wt/Tests/Data/RichardsMechanics"
OGS = "/Users/vinaykumar/git/build/dsm_swelling_full_pdisj-release/bin/ogs"
NEW = RM + "/ANCHORS_MS33_FASTEX_2026-06-02"
KREL = 0.1     # user 2026-06-02: CONSTANT macro k_rel = 0.1 (temporary; reverts to orthotropic power law later)
MIN_DT = 1e-4  # dt-refine: deepen adaptive-retry floor so it steps through the ramp->constant corner
DIXON = {"I_dd1400": 4.922, "I_dd1600": 14.161, "I_dd1800": 40.86}   # EMDD=rho_d Dixon targets [MPa]

PRJ = {
    "LE": {
        "I_dd1400": "ANCHORS_MS33_ModelI/ms33_modelI_dd1400.prj",
        "I_dd1600": "ANCHORS_MS33_ModelI/ms33_modelI_dd1600.prj",
        "I_dd1800": "ANCHORS_MS33_ModelI/ms33_modelI_dd1800.prj",
        "III": "ANCHORS_MS33_ModelIII/ms33_modelIII_gap2mm.prj",
        "IV": "ANCHORS_MS33_ModelIV/ms33_modelIV_pellets.prj",
        "VII": "ANCHORS_MS33_LE_PER_DD/ModelVII/ms33_modelVII_freeswelling_le_perdd.prj",
    },
    "MCC": {
        "I_dd1400": "ANCHORS_MS33_MCC_NATIVE/ModelI_dd1400/ms33_modelI_dd1400_mcc_native.prj",
        "I_dd1600": "ANCHORS_MS33_MCC_NATIVE/ModelI_dd1600/ms33_modelI_dd1600_mcc_native.prj",
        "I_dd1800": "ANCHORS_MS33_MCC_NATIVE/ModelI_dd1800/ms33_modelI_dd1800_mcc_native.prj",
        "III": "ANCHORS_MS33_MCC_NATIVE/ModelIII/ms33_modelIII_gap2mm_mcc_native.prj",
        "IV": "ANCHORS_MS33_MCC_NATIVE/ModelIV/ms33_modelIV_pellets_mcc_native.prj",
        "VII": "ANCHORS_MS33_MCC_NATIVE/ModelVII/ms33_modelVII_freeswelling_mcc_native.prj",
    },
}


def meanp(s):
    return -np.atleast_2d(s)[:, :3].sum(1) / 3 / 1e6


def nidx(p, r, z):
    return int(np.argmin((p[:, 0] - r) ** 2 + (p[:, 1] - z) ** 2))


def main():
    suite, key, alpha = sys.argv[1], sys.argv[2], float(sys.argv[3])
    Kov = sys.argv[4] if len(sys.argv) > 4 else None   # "34368" (set all) or "x0.977" (scale all, keeps ratios)
    src = RM + "/" + PRJ[suite][key]
    dest = NEW + "/" + suite + "/" + key
    os.makedirs(dest, exist_ok=True)
    txt = open(src, encoding="latin-1").read()
    srcdir = os.path.dirname(src)
    for m in re.findall(r"<mesh[^>]*>([^<]+)</mesh>", txt):
        if not m.endswith(".vtu"):
            continue
        ms = os.path.normpath(os.path.join(srcdir, m))
        shutil.copy(ms, dest + "/" + os.path.basename(m))
        txt = txt.replace(">" + m + "<", ">" + os.path.basename(m) + "<")
    txt = re.sub(r"(<mass_exchange_coefficient>)[^<]+(</mass_exchange_coefficient>)",
                 rf"\g<1>{alpha:g}\g<2>", txt)
    # relative permeability -> Constant 0.1 (temporary; replaces the GeneralizedPower block)
    txt = re.sub(r"<name>relative_permeability</name>.*?</property>",
                 f"<name>relative_permeability</name>\n                    <type>Constant</type>\n"
                 f"                    <value>{KREL:g}</value>\n                </property>",
                 txt, count=1, flags=re.DOTALL)
    txt = re.sub(r"(<minimum_dt>)[^<]+(</minimum_dt>)", rf"\g<1>{MIN_DT:g}\g<2>", txt)  # dt-refine corner
    # strip stale tag removed in the full-pdisj branch (some older PRJs, e.g. LE_PER_DD VII, still carry it)
    txt = re.sub(r"\n[ \t]*<accumulate_swelling_contributions>[^<]*</accumulate_swelling_contributions>", "", txt)
    if Kov is not None:
        if Kov.startswith("x"):     # scale all K-tags (preserves IV clay/pellet ratio)
            fac = float(Kov[1:])
            txt = re.sub(r"(<vdw_augmentation_prefactor>)([^<]+)(</vdw_augmentation_prefactor>)",
                         lambda m: f"{m.group(1)}{float(m.group(2)) * fac:g}{m.group(3)}", txt)
        else:                       # set all K-tags to an absolute value
            txt = re.sub(r"(<vdw_augmentation_prefactor>)[^<]+(</vdw_augmentation_prefactor>)",
                         rf"\g<1>{float(Kov):g}\g<2>", txt)
    prj = dest + "/" + os.path.basename(src)
    open(prj, "w", encoding="latin-1").write(txt)
    out = dest + "/out"
    shutil.rmtree(out, ignore_errors=True); os.makedirs(out)
    env = dict(os.environ, OMP_NUM_THREADS="1")
    try:
        r = subprocess.run([OGS, "-o", out, prj], cwd=dest, capture_output=True, text=True, timeout=1800, env=env)
        log = r.stdout + r.stderr
    except subprocess.TimeoutExpired as e:
        so = e.stdout.decode() if isinstance(e.stdout, bytes) else (e.stdout or "")
        se = e.stderr.decode() if isinstance(e.stderr, bytes) else (e.stderr or "")
        log = so + se + "\n[TIMEOUT]"
    done = "Simulation completed" in log
    fs = sorted(glob.glob(out + "/*_ts_*.vtu"), key=lambda f: int(re.search(r"_ts_(\d+)_", f).group(1)))
    fail = re.search(r"failed in time step #\d+ at t = ([0-9.]+)", log)
    tail = f" CRASH@{float(fail.group(1))/86400:.1f}d" if fail else ""
    if not fs:
        print(f"{suite}/{key} a={alpha:g} K={Kov} NO OUTPUT{tail}")
        print("\n".join(log.splitlines()[-8:])); return
    p0 = meshio.read(fs[0]).points
    c = nidx(p0, 0.0125, 0.035)
    rows = []
    for f in fs:
        pd = meshio.read(f).point_data
        t = float(re.search(r"_t_([0-9.]+)\.vtu", f).group(1)) / 86400
        rows.append((t, float(meanp(pd["sigma"])[c])))
    pend = rows[-1][1]
    t95 = next((t for t, p in rows if abs(p) >= 0.95 * abs(pend)), None)
    tmax = rows[-1][0]
    late = [p for t, p in rows if t >= 0.75 * tmax]
    flat = (max(late) - min(late)) / abs(pend) * 100 if pend else 0.0
    dx = f" vs Dixon {DIXON[key]:.3f}({(pend-DIXON[key])/DIXON[key]*100:+.1f}%)" if key in DIXON else ""
    print(f"{suite}/{key} a={alpha:g} K={Kov} krel={KREL} done={done} p_end={pend:.3f}MPa{dx} "
          f"t95={t95}d flat={flat:.1f}%{tail}")
    print("  traj: " + ", ".join(f"{t:.0f}:{p:.2f}" for t, p in rows))


if __name__ == "__main__":
    main()
