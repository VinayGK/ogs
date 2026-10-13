"""Post-processing of the DSM element Jacobian FD test (DSMElementJacobianFD).

Reads the CSV files written by
    OGS_DSM_JACFD_CSV=<dir>/results.csv OGS_DSM_JACFD_MATRICES=<dir>/mat \
        bin/testrunner --gtest_filter='RichardsMechanics.DSMElementJacobianFD'
and writes the markdown tables of docs/dsm_jacobian_fd_verification_2026-10.md
plus the attribution fits:

* P = int B^T m N_p w and Q = int B^T m m^T B w of the axisymmetric 1 cm Quad4
  (2x2 Gauss, w = 2 pi r detJ w_q), m = (1, 1, 1, 0) in Kelvin order
  (rr, zz, tt, rz). P is checked against the analytic up block at the wet
  state (alpha chi = 1, p_cap = 0, no body force: J_up = -P).
* err_pu, err_up = J_an - J_fd and err_uu are fitted as c P^T, c P and c Q (least
  squares); c is the missing isotropic stress derivative [Pa per Pa] and
  [Pa per unit eps_v]; the fit residual says whether the error is purely
  volumetric-isotropic.

Usage: python3 dsm_jacobian_fd_analysis.py <dir>
"""

import csv
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

N_P, N_U = 4, 8
NODES = np.array([[0.0, 0.0], [0.01, 0.0], [0.01, 0.01], [0.0, 0.01]])


def shape(xi, eta):
    n = 0.25 * np.array(
        [
            (1 - xi) * (1 - eta),
            (1 + xi) * (1 - eta),
            (1 + xi) * (1 + eta),
            (1 - xi) * (1 + eta),
        ]
    )
    dn = 0.25 * np.array(
        [
            [-(1 - eta), (1 - eta), (1 + eta), -(1 + eta)],
            [-(1 - xi), -(1 + xi), (1 + xi), (1 - xi)],
        ]
    )
    return n, dn


def patterns():
    g = 1 / np.sqrt(3)
    P = np.zeros((N_U, N_P))
    Q = np.zeros((N_U, N_U))
    m = np.array([1.0, 1.0, 1.0, 0.0])
    for xi in (-g, g):
        for eta in (-g, g):
            n, dn = shape(xi, eta)
            jac = dn @ NODES
            dndx = np.linalg.solve(jac, dn)
            r = n @ NODES[:, 0]
            w = 2 * np.pi * r * np.linalg.det(jac)
            B = np.zeros((4, N_U))
            B[0, 0:4] = dndx[0]
            B[1, 4:8] = dndx[1]
            B[2, 0:4] = n / r
            B[3, 0:4] = dndx[1] / np.sqrt(2)
            B[3, 4:8] = dndx[0] / np.sqrt(2)
            P += np.outer(B.T @ m, n) * w
            Q += np.outer(B.T @ m, B.T @ m) * w
    return P, Q


def fit(err, pat):
    c = float(np.sum(err * pat) / np.sum(pat * pat))
    res = np.linalg.norm(err - c * pat) / max(np.linalg.norm(err), 1e-300)
    return c, res


def main(d):
    d = Path(d)
    P, Q = patterns()
    with (d / "results.csv").open() as f:
        rows = list(csv.DictReader(f))
    with (d / "results.csv.state.csv").open() as f:
        state = {(r["switch_set"], r["state"]): r for r in csv.DictReader(f)}
    # check P against the wet analytic up block
    A = np.loadtxt(
        d / "mat" / "deck_step_klevel_wet_ceiling_analytic.csv", delimiter=","
    )
    chk = np.linalg.norm(A[N_P:, :N_P] + P) / np.linalg.norm(P)
    print(f"P check at wet state: ||J_up + P|| / ||P|| = {chk:.2e}\n")

    sets = []
    states = []
    for r in rows:
        if r["switch_set"] not in sets:
            sets.append(r["switch_set"])
        if r["state"] not in states:
            states.append(r["state"])
    tab = defaultdict(dict)
    for r in rows:
        key = (
            r["switch_set"],
            r["state"],
            r["h_mode"] + ("+relag5" if r["relag"] != "0" else ""),
        )
        tab[key][r["block"]] = r

    def e(x):
        return f"{float(x):.1e}"

    for hm in ["abs", "rel1e-2", "rel1e-7", "abs+relag5"]:
        print(
            f"\n### h = {hm}: max abs error / max rel error (block-normalised) per block\n"
        )
        print(
            "| switch set | state | kkt IPs | pp abs | pp rel | pu abs | pu rel | up abs | up rel | uu abs | uu rel | ‖·‖ rel |"
        )
        print("|---|---|---|---|---|---|---|---|---|---|---|---|")
        for s in sets:
            for st in states:
                b = tab.get((s, st, hm))
                if not b:
                    continue
                cells = []
                for blk in ["pp", "pu", "up", "uu"]:
                    cells += [e(b[blk]["max_abs_err"]), e(b[blk]["max_rel_err"])]
                print(
                    f"| {s} | {st} | {b['pp']['kkt_active_ip']} | "
                    + " | ".join(cells)
                    + f" | {e(b['pp']['fro_rel'])} |"
                )

    print("\n### block magnitudes (max |J_fd|, h = abs)\n")
    print("| switch set | state | pp | pu | up | uu |")
    print("|---|---|---|---|---|---|")
    for s in sets:
        for st in states:
            b = tab.get((s, st, "abs"))
            if b:
                print(
                    f"| {s} | {st} | "
                    + " | ".join(
                        e(b[k]["max_abs_fd"]) for k in ["pp", "pu", "up", "uu"]
                    )
                    + " |"
                )

    print("\n### attribution fits (h = abs)\n")
    print(
        "| switch set | state | c_pu = fit err_pu / P^T | resid | c_up = fit err_up / P [Pa/Pa] | resid | -dsigma_sw/dp (IP FD) | c_uu = fit err_uu / Q [Pa] | resid | dsigma_sw/deps_v (IP FD) | dn_l/deps_v | (dsigma_sw/dp)/(dn_l/dp) dn_l/deps_v |"
    )
    print("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for s in sets:
        for st in states:
            fa = d / "mat" / f"{s}_{st}_analytic.csv"
            ff = d / "mat" / f"{s}_{st}_fd_abs.csv"
            if not fa.exists():
                continue
            A = np.loadtxt(fa, delimiter=",")
            F = np.loadtxt(ff, delimiter=",")
            err = A - F
            cpu, rpu = fit(err[:N_P, N_P:], P.T)
            cup, rup = fit(err[N_P:, :N_P], P)
            cuu, ruu = fit(err[N_P:, N_P:], Q)
            sr = state[(s, st)]
            dsdp = float(sr["dsigma_sw_dp"])
            dndp = float(sr["dn_l_dp"])
            dnde = float(sr["dn_l_deps_v"])
            nl_ch = dsdp / dndp * dnde if dndp != 0 else 0.0
            print(
                f"| {s} | {st} | {cpu:.4e} | {rpu:.1e} | {cup:.4e} | {rup:.1e} | {-dsdp:.4e} | {cuu:.4e} | {ruu:.1e} | "
                f"{float(sr['dsigma_sw_deps_v']):.4e} | {dnde:.4e} | {nl_ch:.4e} |"
            )


if __name__ == "__main__":
    main(sys.argv[1])
