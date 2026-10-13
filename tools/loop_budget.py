#!/usr/bin/env python3
"""Stress-work budget of a one-element DSM closed loop (2026-10).

Reads the VTU frames of one run (via its .pvd), averages every nodal field over
the nodes of the mesh, and evaluates over the window [t_start, t_end] (days):

  eps_v      from the porosity, d eps_v = d phi / (1 - phi)  (and, as a check,
             from the written epsilon, eps_rr + eps_zz + eps_tt)
  sigma_sw   mean of the normal components of `swelling_stress` (OGS sign,
             tension positive; the eigenstress is isotropic)
  rho_lR     `micro_liquid_density` if written, else --rho-micro
  mu_lR      = -micro_pressure / rho_lR          (micro_pressure = -rho_lR mu_lR)
  mu_LR      = p_L / rho_LR                      (rho_LR = --rho-macro, 1000)
  m_mic      = rho_lR (1 - phi_M) n_l            (phi_M = transport_porosity)
  m_tot      = rho_LR S_L phi_M + m_mic
  W_sigma    = closed int sigma_sw d eps_v
  W_mu       = closed int mu_lR d m_mic          (= closed int (rho mu_lR) d m_mic / rho)
  W_mu_nl    = closed int (1 - phi_M) rho_lR mu_lR d n_l
  R_pair     = W_sigma + W_mu      (residual of the pair (sigma_sw, mu_lR))
  R_pair_nl  = W_sigma + W_mu_nl
  D          = closed int (mu_LR - mu_lR) d m_mic   (exchange dissipation)
  W_chem     = closed int mu_LR d m_tot
  budget     = W_sigma + W_chem - D
All integrals: trapezoidal rule over consecutive frames. Units J/m3 (= Pa).
Also: start/end differences, per-leg integrals, and the stationarity at the
end of every hold (relative change over the last output interval).

Usage: loop_budget.py RUN.pvd --t-start D [--t-end D] [--legs t0,t1,...]
       [--json out.json]
"""

import argparse
import json
import xml.etree.ElementTree as ET
from pathlib import Path

import meshio
import numpy as np

DAY = 86400.0


def read_frames(pvd):
    pvd = Path(pvd)
    root = ET.parse(pvd).getroot()
    rows = []
    for ds in root.iter("DataSet"):
        rows.append((float(ds.get("timestep")), pvd.parent / ds.get("file")))
    rows.sort(key=lambda r: r[0])
    # keep the last frame written for a given time (OGS writes once per time)
    out = {}
    for t, f in rows:
        out[t] = f
    times = np.array(sorted(out))
    fields = {}
    for i, t in enumerate(times):
        m = meshio.read(out[t])
        for k, v in m.point_data.items():
            v = np.asarray(v, dtype=float)
            mean = v.mean(axis=0)
            fields.setdefault(k, []).append(np.atleast_1d(mean))
    return times / DAY, {k: np.array(v) for k, v in fields.items()}


def trap(f, g):
    return float(np.sum(0.5 * (f[1:] + f[:-1]) * np.diff(g)))


def evaluate(times, F, t_start, t_end, rho_micro, rho_macro, legs=None):
    sel = (times >= t_start - 1e-9) & (times <= t_end + 1e-9)
    t = times[sel]
    g = {k: v[sel] for k, v in F.items()}
    phi = g["porosity"][:, 0]
    phiM = g["transport_porosity"][:, 0]
    nl = g["micro_water_content"][:, 0]
    pL = g["pressure"][:, 0]
    SL = g["saturation"][:, 0]
    pm = g["micro_pressure"][:, 0]
    sw = g["swelling_stress"][:, :3].mean(axis=1)
    sig = g["sigma"][:, :3].mean(axis=1)
    if "micro_liquid_density" in g:
        rho = g["micro_liquid_density"][:, 0]
        rho_src = "micro_liquid_density (written)"
    else:
        rho = np.full_like(phi, rho_micro)
        rho_src = f"constant {rho_micro} kg/m3 (micro_liquid_density not written)"
    deps = np.diff(phi) / (1.0 - 0.5 * (phi[1:] + phi[:-1]))
    eps = np.concatenate([[0.0], np.cumsum(deps)])
    eps_w = g["epsilon"][:, :3].sum(axis=1) if "epsilon" in g else None
    mu_l = -pm / rho
    mu_L = pL / rho_macro
    m_mic = rho * (1.0 - phiM) * nl
    m_tot = rho_macro * SL * phiM + m_mic

    def integrals(a, b):
        s = slice(a, b + 1)
        Ws = trap(sw[s], eps[s])
        Wm = trap(mu_l[s], m_mic[s])
        Wmn = trap((1.0 - phiM[s]) * rho[s] * mu_l[s], nl[s])
        D = trap(mu_L[s] - mu_l[s], m_mic[s])
        Wc = trap(mu_L[s], m_tot[s])
        return dict(W_sigma=Ws, W_mu=Wm, W_mu_nl=Wmn, R_pair=Ws + Wm,
                    R_pair_nl=Ws + Wmn, D=D, W_chem=Wc, budget=Ws + Wc - D)

    res = dict(frames=int(len(t)), t_start=float(t[0]), t_end=float(t[-1]),
               rho_micro_source=rho_src, rho_macro=rho_macro)
    res.update(integrals(0, len(t) - 1))
    # D increments per frame, to check the sign along the path
    dD = (0.5 * ((mu_L - mu_l)[1:] + (mu_L - mu_l)[:-1])) * np.diff(m_mic)
    res["D_increment_min"] = float(dD.min())
    res["D_increment_negative_sum"] = float(dD[dD < 0].sum())

    def rel(a, b, scale):
        return float(abs(b - a) / max(abs(a), abs(b), scale))

    st = dict(n_l=(nl, 1e-12), phi=(phi, 1e-12), phi_M=(phiM, 1e-12),
              p_L=(pL, 1.0), sigma_sw=(sw, 1.0), sigma_mean=(sig, 1.0),
              micro_pressure=(pm, 1.0), rho_lR=(rho, 1e-12), m_mic=(m_mic, 1e-12))
    res["start"] = {k: float(v[0][0]) for k, v in st.items()}
    res["end"] = {k: float(v[0][-1]) for k, v in st.items()}
    res["start_end_rel"] = {k: rel(v[0][0], v[0][-1], v[1]) for k, v in st.items()}
    res["start_end_abs"] = {k: float(v[0][-1] - v[0][0]) for k, v in st.items()}
    res["eps_v_end_from_phi"] = float(eps[-1])
    res["eps_v_max_from_phi"] = float(eps.max())
    if eps_w is not None:
        res["eps_v_written_start_end"] = [float(eps_w[0]), float(eps_w[-1])]
        res["eps_v_written_max"] = float(eps_w.max())
        res["W_sigma_written_eps"] = trap(sw, eps_w)
    if legs:
        idx = [int(np.argmin(abs(t - x))) for x in legs]
        res["legs"] = []
        for a, b in zip(idx[:-1], idx[1:]):
            d = integrals(a, b)
            d["t0"], d["t1"] = float(t[a]), float(t[b])
            res["legs"].append(d)
    return res, t, st


def stationarity(times, F, hold_ends):
    """Relative change over the last output interval before each hold end."""
    keys = dict(micro_water_content=1e-12, porosity=1e-12,
                transport_porosity=1e-12, micro_pressure=1.0,
                swelling_stress=1.0, sigma=1.0, micro_liquid_density=1e-12,
                saturation=1e-12, pressure=1.0)
    out = []
    for te in hold_ends:
        i = int(np.argmin(abs(times - te)))
        row = dict(t=float(times[i]), dt=float(times[i] - times[i - 1]))
        worst = 0.0
        for k, sc in keys.items():
            if k not in F:
                continue
            a, b = F[k][i - 1], F[k][i]
            r = float(np.max(np.abs(b - a) / np.maximum(np.maximum(abs(a), abs(b)), sc)))
            row[k] = r
            worst = max(worst, r)
        row["max"] = worst
        out.append(row)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pvd")
    ap.add_argument("--t-start", type=float, required=True)
    ap.add_argument("--t-end", type=float, default=None)
    ap.add_argument("--legs", default=None, help="comma list of leg boundaries [d]")
    ap.add_argument("--hold-ends", default=None, help="comma list [d]")
    ap.add_argument("--rho-micro", type=float, default=1000.0)
    ap.add_argument("--rho-macro", type=float, default=1000.0)
    ap.add_argument("--json", default=None)
    a = ap.parse_args()
    times, F = read_frames(a.pvd)
    t_end = a.t_end if a.t_end is not None else times[-1]
    legs = [float(x) for x in a.legs.split(",")] if a.legs else None
    res, _, _ = evaluate(times, F, a.t_start, t_end, a.rho_micro, a.rho_macro, legs)
    if a.hold_ends:
        res["stationarity"] = stationarity(times, F, [float(x) for x in a.hold_ends.split(",")])
    s = json.dumps(res, indent=1)
    print(s)
    if a.json:
        Path(a.json).write_text(s)


if __name__ == "__main__":
    main()
