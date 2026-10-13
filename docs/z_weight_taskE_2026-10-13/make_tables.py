#!/usr/bin/env python3
"""Markdown tables for RESULTS.md from the *_summary.json files written by extract_probes.py.
Usage: make_tables.py <csv dir>"""
import json, os, sys

d = sys.argv[1]
RUNS = [("OFF", "ms33_modelVII_freeswelling"), ("level", "ms33_modelVII_freeswelling_level"),
        ("level+Z", "ms33_modelVII_freeswelling_levelz"),
        ("OFF", "ms33_modelVII_freeswelling_ladder550"), ("level", "ms33_modelVII_freeswelling_ladder550_level"),
        ("level+Z", "ms33_modelVII_freeswelling_ladder550_levelz")]
PROBES = ["top", "centre", "bottom"]
S = {p: json.load(open(os.path.join(d, p + "_summary.json"))) for _, p in RUNS}


def row(lab, probe, r):
    return (f"| {lab} | {probe} | {r['dry_density_solid']:.2f} | {r['void_ratio']:.4f} | {r['transport_porosity']:.4g} | "
            f"{r['nl_over_phi']:.4f} | {r['swelling_stress_mean_comp_pos_Pa']/1e6:.4f} | {r['sigma_mean_comp_pos_Pa']/1e6:.4f} |")


HDR = ("| run | probe | rho_d [kg/m3] | e [-] | phi_M [-] | n_l/phi [-] | -tr(sigma_sw)/3 [MPa] | -tr(sigma)/3 [MPa] |\n"
       "|---|---|---|---|---|---|---|---|")
for group, sel, times in (("240-d decks", RUNS[:3], ["200d", "240d"]),
                          ("ladder550 decks", RUNS[3:], ["200d", "240d", "550d"])):
    for t in times:
        print(f"\n#### {group}, t = {t[:-1]} d\n")
        print(HDR)
        for lab, p in sel:
            for pr in PROBES:
                at = S[p]["probes"][pr]["at"]
                print(row(lab, pr, at[t]) if t in at else f"| {lab} | {pr} | (no frame) | | | | | |")
    print(f"\n#### {group}, end of run and first frame with phi_M <= 1e-6\n")
    print("| run | probe | node (r, z) [m] | t_end [d] | rho_d end [kg/m3] | e end [-] | phi_M end [-] | phi_M <= 1e-6 first at [d] |\n|---|---|---|---|---|---|---|---|")
    for lab, p in sel:
        for pr in PROBES:
            q = S[p]["probes"][pr]; e = q["end"]; z = q["phi_M_zero_first_frame_time_s"]
            print(f"| {lab} | {pr} | {q['node']} ({q['r']:g}, {q['z']:g}) | {e['time_d']:g} | {e['dry_density_solid']:.2f} | "
                  f"{e['void_ratio']:.4f} | {e['transport_porosity']:.4g} | {('%g' % (z/86400)) if z is not None else 'not reached'} |")

# Differences (derived by subtraction of the rows above): level - OFF and (level+Z) - level
Q = [("dry_density_solid", "rho_d [kg/m3]", 1.0, "%.3f"), ("void_ratio", "e [-]", 1.0, "%.5f"),
     ("nl_over_phi", "n_l/phi [-]", 1.0, "%.2e"),
     ("swelling_stress_mean_comp_pos_Pa", "-tr(sigma_sw)/3 [kPa]", 1e-3, "%.2f"),
     ("sigma_mean_comp_pos_Pa", "-tr(sigma)/3 [kPa]", 1e-3, "%.3f")]
print("\n#### Differences (derived: subtraction of the computed rows)\n")
print("| deck set | t [d] | probe | pair | " + " | ".join(q[1] for q in Q) + " |\n|---|---|---|---|" + "---|" * len(Q))
for sel, times in ((RUNS[:3], ["200d", "240d"]), (RUNS[3:], ["200d", "240d", "550d"])):
    off, lev, lz = [S[p] for _, p in sel]
    name = "ladder550" if "ladder" in sel[0][1] else "240-d"
    for t in times:
        for pr in PROBES:
            a, b, c = (x["probes"][pr]["at"][t] for x in (off, lev, lz))
            for pair, (u, v) in (("level - OFF", (b, a)), ("level+Z - level", (c, b))):
                print(f"| {name} | {t[:-1]} | {pr} | {pair} | " + " | ".join(q[3] % ((u[q[0]] - v[q[0]]) * q[2]) for q in Q) + " |")
