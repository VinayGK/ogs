#!/usr/bin/env python3
"""Run the full LE + MCC MS33 suites with the calibrated FASTEX parameters
(krel=0.1 constant, alpha=3e-13, K per dd / scaled by the dd1600 factor),
4 OGS runs in parallel (each single-threaded), writing a live + final summary.

Calibration (2026-06-02): krel=0.1 (constant, temporary), alpha_exchange=3e-13
(transient: Model I t95<20 d), K from Dixon EMDD=rho_d (vdw_augmentation_prefactor):
  dd1400=34368, dd1600=83377, dd1800=224610 J/kg; III/IV/VII inherit dd1.6 (x0.97729).
Known: dd1800 saturated-corner linear-solver singularity (path 0-20 d + Dixon captured);
       MCC III/IV/VII expected to fail at the cam-clay tension apex (documented).
Run (background): python3 run_suite.py
"""
import subprocess
import concurrent.futures
import os
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ALPHA = "3e-13"
# (suite, key, Kspec): Model I = absolute calibrated K; III/IV/VII = scale by dd1600 factor
SUITE = [
    ("LE", "I_dd1400", "34368"), ("LE", "I_dd1600", "83377"), ("LE", "I_dd1800", "224610"),
    ("LE", "III", "x0.97729"), ("LE", "IV", "x0.97729"), ("LE", "VII", "x0.97729"),
    ("MCC", "I_dd1400", "34368"), ("MCC", "I_dd1600", "83377"), ("MCC", "I_dd1800", "224610"),
    ("MCC", "III", "x0.97729"), ("MCC", "IV", "x0.97729"), ("MCC", "VII", "x0.97729"),
]


def run(job):
    suite, key, K = job
    t0 = time.time()
    try:
        r = subprocess.run(["python3", HERE + "/run_fastex.py", suite, key, ALPHA, K],
                           capture_output=True, text=True, timeout=2000)
        out = r.stdout
    except subprocess.TimeoutExpired:
        out = "[suite-timeout]"
    line = [l for l in out.splitlines() if l.startswith(suite + "/")]
    res = (line[0] if line else f"{suite}/{key} NO RESULT ({out.strip()[-150:]})") + f"  [{time.time()-t0:.0f}s]"
    print(res, flush=True)
    return res


if __name__ == "__main__":
    print(f"== FASTEX suite start ({len(SUITE)} models, 4-parallel) ==", flush=True)
    results = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:
        for res in ex.map(run, SUITE):
            results.append(res)
    open(HERE + "/SUITE_SUMMARY.txt", "w").write("\n".join(results) + "\n")
    print("== FASTEX suite DONE -> SUITE_SUMMARY.txt ==", flush=True)
