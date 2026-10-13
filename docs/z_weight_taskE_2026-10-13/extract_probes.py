#!/usr/bin/env python3
"""Task E probe extraction (CLOUD_PROMPT_Z_weight_taskE_2026-10-13.md, item 4).

Usage: extract_probes.py <run output dir> <prefix> <csv out dir>

Probes: nearest mesh node (reference/undeformed coordinates, as written in the VTU) to
r = 0 and z = 0.070, 0.04025, 0.0105 m. Quantities per output frame (all from written point data):
  rho_d      = dry_density_solid                                  [kg/m3]
  e          = 2780/rho_d - 1   (rho_s = 2780 kg/m3, task prompt) [-]
  phi_M      = transport_porosity                                 [-]
  nl_over_phi= micro_water_content / porosity                     [-]
  p_sw_mean  = -tr(swelling_stress)/3  (compression positive)     [Pa]
  p_mean     = -tr(sigma)/3            (compression positive)     [Pa]
Tensors are written with 4 Kelvin components (rr, zz, thetatheta, rz); the trace is the sum of the first three.
"""
import csv, json, os, sys
import xml.etree.ElementTree as ET
import numpy as np
import meshio

RHO_S = 2780.0
PROBES = {"top": 0.070, "centre": 0.04025, "bottom": 0.0105}
PHI_ZERO = 1e-6
REPORT_TIMES = {"200d": 17280000.0, "240d": 20736000.0, "550d": 47520000.0}


def frames(run_dir, prefix):
    pvd = os.path.join(run_dir, prefix + ".pvd")
    root = ET.parse(pvd).getroot()
    out = []
    for ds in root.iter("DataSet"):
        out.append((float(ds.get("timestep")), os.path.join(run_dir, ds.get("file"))))
    return out


def main():
    run_dir, prefix, csv_dir = sys.argv[1:4]
    os.makedirs(csv_dir, exist_ok=True)
    fr = frames(run_dir, prefix)
    m0 = meshio.read(fr[0][1])
    pts = m0.points
    idx = {k: int(np.argmin(np.hypot(pts[:, 0] - 0.0, pts[:, 1] - z))) for k, z in PROBES.items()}
    rows = {k: [] for k in PROBES}
    for t, f in fr:
        m = meshio.read(f)
        pd = m.point_data
        for k, i in idx.items():
            rho = float(pd["dry_density_solid"][i])
            sw = pd["swelling_stress"][i]
            sg = pd["sigma"][i]
            rows[k].append({
                "time_s": t, "time_d": t / 86400.0,
                "dry_density_solid": rho,
                "void_ratio": RHO_S / rho - 1.0,
                "transport_porosity": float(pd["transport_porosity"][i]),
                "micro_water_content": float(pd["micro_water_content"][i]),
                "porosity": float(pd["porosity"][i]),
                "nl_over_phi": float(pd["micro_water_content"][i]) / float(pd["porosity"][i]),
                "swelling_stress_mean_comp_pos_Pa": -float(sw[0] + sw[1] + sw[2]) / 3.0,
                "sigma_mean_comp_pos_Pa": -float(sg[0] + sg[1] + sg[2]) / 3.0,
            })
    summary = {"run_dir": run_dir, "prefix": prefix, "n_frames": len(fr), "probes": {}}
    for k, i in idx.items():
        fn = os.path.join(csv_dir, f"{prefix}_{k}.csv")
        with open(fn, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[k][0].keys()))
            w.writeheader()
            for r in rows[k]:
                w.writerow({kk: (f"{v:.10g}" if isinstance(v, float) else v) for kk, v in r.items()})
        zero = next((r["time_s"] for r in rows[k] if r["transport_porosity"] <= PHI_ZERO), None)
        at = {}
        for lab, tt in REPORT_TIMES.items():
            hit = [r for r in rows[k] if abs(r["time_s"] - tt) < 1e-3]
            if hit:
                at[lab] = hit[0]
        summary["probes"][k] = {
            "node": i, "r": float(pts[i, 0]), "z": float(pts[i, 1]),
            "end": rows[k][-1],
            "phi_M_zero_first_frame_time_s": zero,
            "at": at,
        }
    with open(os.path.join(csv_dir, f"{prefix}_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1)


if __name__ == "__main__":
    main()
