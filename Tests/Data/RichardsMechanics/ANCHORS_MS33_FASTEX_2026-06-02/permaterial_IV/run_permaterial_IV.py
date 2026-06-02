#!/usr/bin/env python3
"""Per-material Model IV runner (FASTEX 2026-06-02 settings, PER-MEDIUM DSM).

Fixes the single-K FASTEX IV bug: the prior run applied ONE process-level
<potential_exchange> (K=83375.2, n_s=0.5755, n_l0=1.1699e-3 -- all CLAY values) to BOTH
media, i.e. the dd-0.9 pellet was treated as dd-1.6 clay. Here each medium carries its
OWN DSM parameters via the per-medium <medium id="1"> override inside <potential_exchange>
(parsed by CreateRichardsMechanicsProcess.cpp: process-level = defaults/clay, <medium id>
= per-material override inheriting unspecified tags).

  CLAY  (medium id 0, process-level): K=83377 J/kg (FASTEX dd-1.6), n_s=0.5755395683453237,
        n_l0=1.1699e-3.
  PELLET(medium id 1, override):       K=6600 J/kg (CALIBRATED to Dixon dd-0.9 median
        sigma_sat=0.350 MPa: probe gives 0.3481 MPa, -0.6%), n_s=0.3237410071942446,
        n_l0=6.59e-4.

FASTEX run settings (both media): mass_exchange_coefficient alpha=3e-13,
macro relative_permeability=Constant 0.1, minimum_dt=1e-4, stale
<accumulate_swelling_contributions> stripped. OMP_NUM_THREADS=1. 200 d.

Usage: python3 run_permaterial_IV.py <LE|MCC>
"""
import os, re, shutil, subprocess, glob, sys
import numpy as np
import meshio

RM = "/Users/vinaykumar/git/ogs-worktrees/dsm_native_hierarchical_wt/Tests/Data/RichardsMechanics"
OGS = "/Users/vinaykumar/git/build/dsm_swelling_full_pdisj-release/bin/ogs"
NEW = RM + "/ANCHORS_MS33_FASTEX_2026-06-02/permaterial_IV"
KREL = 0.1
MIN_DT = 1e-4
ALPHA = 3e-13
K_CLAY = 83377.0    # FASTEX dd-1.6 fitted K (CALIBRATION_PROVENANCE.md table)
K_PELLET = 6600.0   # CALIBRATED this task: Dixon dd-0.9 median sigma_sat 0.350 MPa -> probe 0.3481 (-0.6%)
NS_PELLET = "0.3237410071942446"
NL0_PELLET = "6.59e-4"

BASE = {
    "LE": "ANCHORS_MS33_ModelIV/ms33_modelIV_pellets.prj",
    "MCC": "ANCHORS_MS33_MCC_NATIVE/ModelIV/ms33_modelIV_pellets_mcc_native.prj",
}

PELLET_OVERRIDE = (
    "                <!-- PER-MATERIAL pellet override (material id 1, rho_d=0.9 g/cm3),\n"
    "                     added 2026-06-02 to fix the single-K FASTEX IV bug. n_s=0.32374\n"
    "                     (=1-phi0_pellet), n_l0=6.59e-4 (dd900 vdW eq vs clay 1.17e-3), and a\n"
    "                     K CALIBRATED to the Dixon (2023) dd-0.9 median sigma_sat=0.350 MPa at the\n"
    "                     production FASTEX settings (alpha=3e-13, k_rel=Constant 0.1): single-element\n"
    "                     pellet probe gives sigma_sat=0.3481 MPa (-0.6%) at K=6600. Inherits other\n"
    "                     DSM params from the process-level (clay) block. -->\n"
    f"                <medium id=\"1\">\n"
    f"                    <micro_solid_volume_fraction_reference>{NS_PELLET}</micro_solid_volume_fraction_reference>\n"
    f"                    <initial_micro_water_content>{NL0_PELLET}</initial_micro_water_content>\n"
    f"                    <vdw_augmentation_prefactor>{K_PELLET:g}</vdw_augmentation_prefactor>\n"
    f"                </medium>\n"
)


def set_constant_krel_all(txt):
    """Replace EVERY relative_permeability property block with Constant KREL (both media)."""
    return re.sub(
        r"<name>relative_permeability</name>.*?</property>",
        ("<name>relative_permeability</name>\n"
         "                    <type>Constant</type>\n"
         f"                    <value>{KREL:g}</value>\n"
         "                </property>"),
        txt, flags=re.DOTALL)


def inject_pellet_override(txt):
    """Insert/replace the <medium id="1"> override inside the process-level <potential_exchange>."""
    # If a per-medium override already exists (MCC base has one), replace it wholesale.
    if re.search(r"<medium id=\"1\">.*?</medium>", txt, flags=re.DOTALL):
        # only the FIRST <medium ...> after <potential_exchange> and before </potential_exchange>
        pe = re.search(r"(<potential_exchange>)(.*?)(</potential_exchange>)", txt, flags=re.DOTALL)
        block = pe.group(2)
        block2 = re.sub(r"[ \t]*<!--[^>]*?Per-medium DSM override.*?-->\n", "", block, flags=re.DOTALL)
        block2 = re.sub(r"[ \t]*<medium id=\"1\">.*?</medium>\n", "", block2, flags=re.DOTALL)
        block2 = block2.rstrip("\n") + "\n" + PELLET_OVERRIDE
        txt = txt[:pe.start()] + pe.group(1) + block2 + pe.group(3) + txt[pe.end():]
        return txt
    # LE base: no override yet -> insert just before </potential_exchange>
    txt = txt.replace("</potential_exchange>", PELLET_OVERRIDE + "            </potential_exchange>", 1)
    return txt


def set_clay_K(txt):
    """Set the process-level (clay) vdw_augmentation_prefactor to K_CLAY.
    The FIRST vdw_augmentation_prefactor in the file is the process-level (clay) one."""
    return re.sub(r"(<vdw_augmentation_prefactor>)[^<]+(</vdw_augmentation_prefactor>)",
                  rf"\g<1>{K_CLAY:g}\g<2>", txt, count=1)


def meanp(s):
    s = np.atleast_2d(s)
    return -s[:, :3].sum(axis=1) / 3.0 / 1e6


def main():
    suite = sys.argv[1]
    src = RM + "/" + BASE[suite]
    dest = NEW + "/" + suite
    os.makedirs(dest, exist_ok=True)
    txt = open(src, encoding="latin-1").read()
    srcdir = os.path.dirname(src)
    # copy meshes local, rewrite refs to basename
    for m in re.findall(r"<mesh[^>]*>([^<]+)</mesh>", txt):
        if not m.endswith(".vtu"):
            continue
        ms = os.path.normpath(os.path.join(srcdir, m))
        shutil.copy(ms, dest + "/" + os.path.basename(m))
        txt = txt.replace(">" + m + "<", ">" + os.path.basename(m) + "<")
    # FASTEX settings
    txt = re.sub(r"(<mass_exchange_coefficient>)[^<]+(</mass_exchange_coefficient>)",
                 rf"\g<1>{ALPHA:g}\g<2>", txt)
    txt = set_constant_krel_all(txt)
    txt = re.sub(r"(<minimum_dt>)[^<]+(</minimum_dt>)", rf"\g<1>{MIN_DT:g}\g<2>", txt)
    txt = re.sub(r"\n[ \t]*<accumulate_swelling_contributions>[^<]*</accumulate_swelling_contributions>", "", txt)
    # PER-MATERIAL DSM
    txt = set_clay_K(txt)
    txt = inject_pellet_override(txt)

    prj = dest + "/" + os.path.basename(src)
    open(prj, "w", encoding="latin-1").write(txt)
    out = dest + "/out"
    shutil.rmtree(out, ignore_errors=True); os.makedirs(out)

    # sanity echo of the DSM block
    pe = re.search(r"<potential_exchange>.*?</potential_exchange>", txt, flags=re.DOTALL).group(0)
    print(f"# {suite} K tags found:", re.findall(r"<vdw_augmentation_prefactor>([^<]+)</vdw_augmentation_prefactor>", pe))
    print(f"# {suite} n_s tags found:", re.findall(r"<micro_solid_volume_fraction_reference>([^<]+)</micro_solid_volume_fraction_reference>", pe))
    print(f"# {suite} n_l0 tags found:", re.findall(r"<initial_micro_water_content>([^<]+)</initial_micro_water_content>", pe))
    print(f"# {suite} alpha:", re.findall(r"<mass_exchange_coefficient>([^<]+)</mass_exchange_coefficient>", txt))
    print(f"# {suite} krel Constant count:", txt.count("<value>0.1</value>"))

    env = dict(os.environ, OMP_NUM_THREADS="1")
    try:
        r = subprocess.run([OGS, "-o", out, prj], cwd=dest, capture_output=True, text=True, timeout=3000, env=env)
        log = r.stdout + r.stderr
    except subprocess.TimeoutExpired as e:
        so = e.stdout.decode() if isinstance(e.stdout, bytes) else (e.stdout or "")
        se = e.stderr.decode() if isinstance(e.stderr, bytes) else (e.stderr or "")
        log = so + se + "\n[TIMEOUT]"
    open(dest + "/run.log", "w").write(log)
    done = "Simulation completed" in log
    fail = re.search(r"failed in time step #(\d+) at t = ([0-9.eE+-]+)", log)
    rej = len(re.findall(r"Time step rejected", log)) + len(re.findall(r"reduce.*time step|repeat.*time step", log))
    fs = sorted(glob.glob(out + "/*_ts_*.vtu"), key=lambda f: int(re.search(r"_ts_(\d+)_", f).group(1)))
    tail = ""
    if fail:
        tail = f"  CRASH step#{fail.group(1)} t={float(fail.group(2))/86400:.2f}d"
    if not fs:
        print(f"{suite} NO OUTPUT done={done}{tail}")
        print("\n".join(log.splitlines()[-12:])); return
    pts = meshio.read(fs[0]).points
    clay = pts[:, 1] >= 0.035; pel = pts[:, 1] < 0.035
    last = fs[-1]; pd = meshio.read(last).point_data
    t_end = float(re.search(r"_t_([0-9.]+)\.vtu", last).group(1)) / 86400
    sig = meanp(pd["sigma"]); dd = np.ravel(pd["dry_density_solid"])
    print(f"{suite} done={done} t_end={t_end:.1f}d steps={len(fs)} rej~{rej}{tail}")
    print(f"  clay-zone:  mean_p={sig[clay].mean():.3f} MPa  rho_d_mean={dd[clay].mean():.1f}  (centre node 0.0125,0.0525)")
    print(f"  pellet-zone:mean_p={sig[pel].mean():.3f} MPa  rho_d_mean={dd[pel].mean():.1f}  (centre node 0.0125,0.0175)")


if __name__ == "__main__":
    main()
