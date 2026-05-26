#!/usr/bin/env python3
"""
Run DSM saturation/desaturation verification ramps as OGS single-element tests.

Replaces the Gauss-point Jupyter notebook
  /Users/vinaykumar/tex/cc2024/examples/DSM/Potentials_penalty.ipynb
by running the corresponding ramp on the OGS native DSM hierarchical
implementation, at the three Model-I calibrated dry densities (1400, 1600,
1800 kg/m^3). The single-element mesh is reused from ANCHORS_MS33_ModelI.

Pipeline:
  1. Generate three PRJ files in this directory from the Model-I templates,
     modifying:
         - initial pressure  : -100 MPa -> -1 MPa  (match notebook IC)
         - time loop         : t_end = 4e8 s, dt = 4e6 s (100 steps)
         - pressure curve    : triangle ramp
                                 t=0     p=-1 MPa  (initial suction)
                                 t=2e8s  p=+1e5 Pa (saturation peak)
                                 t=4e8s  p=-1 MPa  (return)
         - output prefix     : dsm_ramp_dd{dd}
         - output variables  : add DSM secondary variables (micro_water_content,
                               micro_porosity, micro_liquid_density,
                               micro_saturation, micro_pressure,
                               micro_exchange_source)
  2. Run each PRJ with the native DSM binary (dsm_native_hierarchical worktree)
  3. Post-process the resulting PVD/VTU series: extract cell-averaged values
     of pressure, n_l, phi_m, total porosity, S_L, sigma into a CSV per density
  4. Plot three verification figures per density:
         (a) water-content components vs time
         (b) water-content components vs macro potential
         (c) saturation/desaturation hysteresis loop (S_L vs |psi|)

Run:
    python3 run_verification_ramps.py             # generate, run, plot
    python3 run_verification_ramps.py --skip-run  # plot only (re-process VTUs)
    python3 run_verification_ramps.py --only-gen  # only regenerate PRJ files
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
MODEL_I_DIR = HERE.parent / "ANCHORS_MS33_ModelI"
PARITY_DIR = HERE.parent / "ANCHORS_MS33_StrictParity"
TEMPLATE_PRJ = PARITY_DIR / "ms33_dsm_parity_native.prj"   # DSM hierarchical active
OGS_BIN = Path("/Users/vinaykumar/git/build/native-release-omp-sharedcache/bin/ogs")
DENSITIES = (1400, 1600, 1800)

# Per-density Model-I calibrated values (phi0, IntrinsicPermeability0)
MODEL_I_PARAMS = {
    1400: dict(phi0="0.49640287769784175", perm="1.2264011098493141e-20"),
    1600: dict(phi0="0.4244604316546763",  perm="5.870260900123441e-21"),
    1800: dict(phi0="0.3525179856115108",  perm="2.6569572325678573e-21"),
}

# Triangle pressure ramp: saturation half (P_LOW -> P_HIGH) then desaturation
# half (P_HIGH -> P_LOW). Range matches the suction regime over which the
# Model-I calibration is active and the DSM constitutive update varies
# meaningfully (cf. notebook shape, but spanning the model's working range).
P_LOW  = -100.0e6  # Pa, dry-side suction (matches MS33 baseline IC)
P_HIGH = 1.0e5     # Pa, near-atmospheric (saturation peak)
T_HALF = 5.184e6   # s, half-cycle = 60 days (matches MS33 baseline length)
T_END  = 1.0368e7  # s, full saturation+desaturation cycle = 120 days
DT     = 86400.0   # s, 1 day (120 steps), matches MS33 baseline

# Free-top variant flag. When True, removes the top-boundary displacement
# Dirichlet BC (the column is then traction-free at the top and may expand
# vertically) and replaces sigma0 with the Biot-balancing initial effective
# stress sigma0_iso = alpha * chi(p_init) * p_init. The lateral and bottom
# BCs are retained (axisymmetric + bottom support), so the column free-swells
# in y only.
FREE_TOP = True
# Macro saturation at p_LR = P_LOW used to build sigma0; taken from the
# confined verification CSV at t=0 (consistent across the three Model-I
# densities at the same p_L, since the macro retention closure is
# density-invariant in this model).
SAT_AT_PINIT = 0.329

DSM_OUTPUT_VARIABLES = [
    "micro_water_content",
    "micro_porosity",
    "micro_liquid_density",
    "micro_saturation",
    "micro_pressure",
    "micro_exchange_source",
]


# ----------------------------------------------------------------------------- generation
def generate_prj(dd: int) -> Path:
    """Generate verification-ramp PRJ from the DSM-active native parity template,
    overriding phi0 / permeability with Model-I-calibrated values for this density.
    The DSM hierarchical path is active (potential_exchange, vdW Hamaker, mass
    exchange) and the micro_* secondary variables are already in the output list."""
    target = HERE / f"dsm_ramp_dd{dd}.prj"
    text = TEMPLATE_PRJ.read_text()
    params = MODEL_I_PARAMS[dd]

    # 1. output prefix
    text = text.replace("ms33_dsm_parity_native", f"dsm_ramp_dd{dd}")

    # 2. initial pressure -100e6 -> -1e6
    text = re.sub(
        r"(<name>pressure_ic</name>\s*<type>Constant</type>\s*<value>)[^<]+(</value>)",
        rf"\g<1>{P_LOW}\g<2>",
        text,
    )

    # 3. phi0 -> Model-I value at this density
    text = re.sub(
        r"(<name>phi0</name>\s*<type>Constant</type>\s*<value>)[^<]+(</value>)",
        rf"\g<1>{params['phi0']}\g<2>",
        text,
    )

    # 4. IntrinsicPermeability0 -> Model-I value at this density
    text = re.sub(
        r"(<name>IntrinsicPermeability0</name>\s*<type>Constant</type>\s*<value>)[^<]+(</value>)",
        rf"\g<1>{params['perm']}\g<2>",
        text,
    )

    # 5. time-stepping: t_end + delta_t (FixedTimeStepping single pair)
    n_steps = int(T_END / DT)
    text = re.sub(
        r"(<time_stepping>.*?<t_initial>0</t_initial>\s*<t_end>)[^<]+(</t_end>.*?<pair>\s*<repeat>)[0-9]+(</repeat>\s*<delta_t>)[^<]+(</delta_t>\s*</pair>)",
        rf"\g<1>{T_END:.0f}\g<2>{n_steps}\g<3>{DT:.0f}\g<4>",
        text,
        flags=re.DOTALL,
    )

    # 6. output frequency: every 5 steps
    each = 5
    repeat_out = n_steps // each
    text = re.sub(
        r"(<output>.*?<timesteps>\s*<pair>\s*<repeat>)[0-9]+(</repeat>\s*<each_steps>)[0-9]+(</each_steps>)",
        rf"\g<1>{repeat_out}\g<2>{each}\g<3>",
        text,
        flags=re.DOTALL,
    )

    # 7. curve: triangle ramp (saturation + desaturation)
    new_curve = (
        "        <curve>\n"
        "            <name>ms33_pressure_release</name>\n"
        f"            <coords>0 {T_HALF:.0f} {T_END:.0f}</coords>\n"
        f"            <values>{P_LOW} {P_HIGH} {P_LOW}</values>\n"
        "        </curve>"
    )
    text = re.sub(
        r"<curve>\s*<name>ms33_pressure_release</name>.*?</curve>",
        new_curve,
        text,
        flags=re.DOTALL,
    )

    # 8. Free-top variant: remove the top-boundary displacement Dirichlet BC so
    #    the column can expand upward. Lateral and bottom BCs (axisymmetric +
    #    bottom support) are retained. Set sigma0 to balance the initial Biot
    #    pressure force so the column starts in mechanical equilibrium under the
    #    high initial suction (memory rule: sigma0_eff = alpha * chi * p_L_init).
    if FREE_TOP:
        # Remove the entire top-boundary block from the displacement process_variable.
        top_bc_pat = (
            r"\s*<boundary_condition>\s*"
            r"<mesh>square_1x1_quad_1e0_top</mesh>\s*"
            r"<type>Dirichlet</type>\s*"
            r"<component>1</component>\s*"
            r"<parameter>zero</parameter>\s*"
            r"</boundary_condition>"
        )
        text, nsub = re.subn(top_bc_pat, "", text, count=1)
        assert nsub == 1, "top displacement BC pattern not found"

        # Replace sigma0 with a Biot-balancing isotropic initial effective stress.
        # At p_L_init = -100 MPa the macro saturation chi ~ 0.33 (parity case);
        # alpha=1, so sigma0_eff_iso = 0.33 * (-100e6) = -33e6 Pa.
        sigma0_iso = 1.0 * SAT_AT_PINIT * P_LOW   # negative (compressive)
        text = re.sub(
            r"(<name>sigma0</name>\s*<type>Function</type>)\s*"
            r"<expression>0</expression>\s*"
            r"<expression>0</expression>\s*"
            r"<expression>0</expression>\s*"
            r"<expression>0</expression>",
            (
                r"\g<1>" + "\n"
                f"            <expression>{sigma0_iso:.6e}</expression>\n"
                f"            <expression>{sigma0_iso:.6e}</expression>\n"
                f"            <expression>{sigma0_iso:.6e}</expression>\n"
                f"            <expression>0</expression>"
            ),
            text,
            count=1,
        )

    target.write_text(text)
    return target


# ----------------------------------------------------------------------------- run
def run_one(prj: Path) -> bool:
    print(f"[run] {prj.name}", flush=True)
    res = subprocess.run(
        [str(OGS_BIN), prj.name, "-o", "."],
        cwd=prj.parent,
        capture_output=True,
        text=True,
    )
    log = prj.with_suffix(".log")
    log.write_text(res.stdout + "\n--- STDERR ---\n" + res.stderr)
    if res.returncode != 0:
        print(f"  FAILED -> see {log}")
        return False
    print(f"  ok ({len([l for l in res.stdout.splitlines() if 'Time step' in l])//2} steps)")
    return True


# ----------------------------------------------------------------------------- post
def extract_series(dd: int) -> Path:
    """Read the PVD for this density, extract scalars per timestep, write CSV
    next to it. For a single-element axisymmetric quad, point_data is uniform
    across the 4 corner nodes, so np.mean is exact."""
    import xml.etree.ElementTree as ET
    import meshio
    import numpy as np

    pvd = HERE / f"dsm_ramp_dd{dd}.pvd"
    if not pvd.exists():
        raise FileNotFoundError(pvd)
    tree = ET.parse(pvd)
    rows = []
    for ds in tree.getroot().iter("DataSet"):
        t = float(ds.attrib["timestep"])
        vtu_path = HERE / ds.attrib["file"]
        m = meshio.read(vtu_path)

        def pmean(name):
            return float(np.mean(m.point_data[name])) if name in m.point_data else float("nan")

        row = {
            "t": t,
            "pressure": pmean("pressure"),
            "saturation": pmean("saturation"),
            "porosity": pmean("porosity"),
            "micro_water_content": pmean("micro_water_content"),
            "micro_porosity": pmean("micro_porosity"),
            "micro_pressure": pmean("micro_pressure"),
            "micro_exchange_source": pmean("micro_exchange_source"),
            "dry_density_solid": pmean("dry_density_solid"),
        }
        # macro porosity phi_M = porosity - micro_porosity (when both present)
        if not (np.isnan(row["porosity"]) or np.isnan(row["micro_porosity"])):
            row["phi_M"] = row["porosity"] - row["micro_porosity"]
        else:
            row["phi_M"] = float("nan")

        # mean stress = trace(sigma)/3 (sigma is Kelvin-like 4-vector for axisym 2D)
        if "sigma" in m.point_data:
            sig = np.asarray(m.point_data["sigma"]).mean(axis=0)
            # Kelvin vector order: sigma_xx, sigma_yy, sigma_zz, sqrt(2)*sigma_xy
            row["sigma_mean"] = float((sig[0] + sig[1] + sig[2]) / 3.0)
            row["sigma_xx"] = float(sig[0]); row["sigma_yy"] = float(sig[1])
            row["sigma_zz"] = float(sig[2])

        # Vertical displacement of the top boundary (free-top runs only):
        # take the maximum y-displacement across the four nodes. The bottom is
        # Dirichlet-clamped at u_y=0, so this is the linear axial swelling
        # strain of the 1 m x 1 m element.
        if "displacement" in m.point_data:
            disp = np.asarray(m.point_data["displacement"])
            row["u_top"] = float(disp[:, 1].max())
        rows.append(row)

    csv_path = HERE / f"dsm_ramp_dd{dd}.csv"
    fieldnames = list(rows[0].keys())
    with csv_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader(); w.writerows(rows)
    print(f"[post] wrote {csv_path.name} ({len(rows)} rows)")
    return csv_path


def plot_density(dd: int) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd

    csv_path = HERE / f"dsm_ramp_dd{dd}.csv"
    df = pd.read_csv(csv_path).sort_values("t").reset_index(drop=True)
    seconds_per_day = 86400.0
    t_days = df["t"] / seconds_per_day

    # (a) Water-content components vs time, with porosities on twin axis
    fig, ax = plt.subplots(figsize=(8, 5))
    ax2 = ax.twinx()
    omega_macro = df["phi_M"] * df["saturation"]
    ax.plot(t_days, df["micro_water_content"], color="tab:orange", lw=2, label=r"$n_l$ (micro)")
    ax.plot(t_days, omega_macro, color="tab:blue", lw=2, label=r"$\phi_M S_L$ (macro)")
    ax.plot(t_days, df["micro_water_content"] + omega_macro, color="tab:red", lw=2, label=r"$\omega^{\rm Total}$")
    ax2.plot(t_days, df["porosity"], color="tab:green", linestyle="--", lw=1.5, label=r"$\phi$ (total)")
    ax2.plot(t_days, df["phi_M"], color="tab:purple", linestyle=":", lw=1.5, label=r"$\phi_M$")
    ax2.plot(t_days, df["micro_porosity"], color="tab:gray", linestyle="-.", lw=1.5, label=r"$\phi_m$")
    ax.set_xlabel("time / days")
    ax.set_ylabel(r"water content $\omega$ / -")
    ax2.set_ylabel(r"porosity $\phi$ / -")
    ax.set_title(rf"DSM verification ramp, $\rho_d = {dd}$ kg/m$^3$ -- water content vs time")
    ax.set_ylim(0, max(0.6, df["porosity"].max() * 1.1))
    ax2.set_ylim(0, max(0.6, df["porosity"].max() * 1.1))
    ax.grid(True, alpha=0.4)
    ax.legend(loc="upper left", fontsize=9, framealpha=0.85)
    ax2.legend(loc="upper right", fontsize=9, framealpha=0.85)
    fig.tight_layout()
    fig.savefig(HERE / f"dsm_ramp_dd{dd}_content_vs_time.png", dpi=140)
    plt.close(fig)

    # (b) Water content vs macro suction |p_macro|
    fig, ax = plt.subplots(figsize=(8, 5))
    psi_MPa = -df["pressure"] * 1e-6        # suction in MPa (>0 in unsaturated)
    psi_plot = np.where(psi_MPa > 1e-3, psi_MPa, np.nan)
    ax.plot(psi_plot, df["micro_water_content"], color="tab:orange", lw=2, marker="o", ms=4, label=r"$n_l$ (micro)")
    ax.plot(psi_plot, df["phi_M"] * df["saturation"], color="tab:blue", lw=2, marker="s", ms=4, label=r"$\phi_M S_L$ (macro)")
    ax.plot(psi_plot, df["micro_water_content"] + df["phi_M"] * df["saturation"], color="tab:red", lw=2, marker="^", ms=4, label=r"$\omega^{\rm Total}$")
    ax.set_xscale("log")
    ax.set_xlabel(r"$|p^{LR}|$ / MPa  (suction)")
    ax.set_ylabel(r"water content $\omega$ / -")
    ax.set_title(rf"DSM verification ramp, $\rho_d = {dd}$ kg/m$^3$ -- water content vs suction")
    ax.grid(True, which="both", alpha=0.4)
    ax.legend(loc="best", fontsize=10, framealpha=0.85)
    fig.tight_layout()
    fig.savefig(HERE / f"dsm_ramp_dd{dd}_content_vs_potential.png", dpi=140)
    plt.close(fig)

    # (c) Hysteresis: macro saturation S_L vs suction; saturation/desaturation branches
    fig, ax = plt.subplots(figsize=(8, 5))
    n_half = int(np.argmax(df["t"] >= T_HALF))
    if n_half == 0:
        n_half = len(df) // 2
    sat_branch = df.iloc[: n_half + 1]
    des_branch = df.iloc[n_half:]
    ax.plot(psi_MPa.iloc[: n_half + 1], sat_branch["saturation"], color="tab:blue", lw=2, marker="o", ms=4, label="saturation branch")
    ax.plot(psi_MPa.iloc[n_half:], des_branch["saturation"], color="tab:red", lw=2, marker="s", ms=4, label="desaturation branch")
    ax.set_xscale("symlog", linthresh=1e-2)
    ax.set_xlabel(r"$|p^{LR}|$ / MPa")
    ax.set_ylabel(r"macro saturation $S_L$ / -")
    ax.set_title(rf"DSM verification ramp, $\rho_d = {dd}$ kg/m$^3$ -- macro saturation hysteresis")
    ax.grid(True, alpha=0.4)
    ax.legend(loc="best", fontsize=10, framealpha=0.85)
    fig.tight_layout()
    fig.savefig(HERE / f"dsm_ramp_dd{dd}_hysteresis.png", dpi=140)
    plt.close(fig)

    # (d') Free-top: axial swelling u_top vs time -- the headline kinematic
    # signature absent in the confined runs.
    if "u_top" in df.columns:
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot(t_days, df["u_top"], color="tab:purple", lw=2, label=r"$u_y$ (top)")
        ax.set_xlabel("time / days")
        ax.set_ylabel(r"axial swelling $u_y^{\,\rm top}$ / m")
        ax.set_title(rf"DSM free-top ramp, $\rho_d = {dd}$ kg/m$^3$ -- axial swelling vs time")
        ax.grid(True, alpha=0.4)
        ax.legend(loc="best", fontsize=10, framealpha=0.85)
        fig.tight_layout()
        fig.savefig(HERE / f"dsm_ramp_dd{dd}_uTop_vs_time.png", dpi=140)
        plt.close(fig)

    # (d) Bonus: mean stress vs time -- swelling-pressure-style trace
    if "sigma_mean" in df.columns:
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot(t_days, -df["sigma_mean"] * 1e-6, color="tab:red", lw=2, label=r"$-\sigma_{\rm mean}$")
        ax.plot(t_days, -df["pressure"] * 1e-6, color="tab:blue", lw=1.5, linestyle="--", label=r"$|p^{LR}|$ (BC)")
        ax.set_xlabel("time / days")
        ax.set_ylabel("stress / MPa")
        ax.set_title(rf"DSM verification ramp, $\rho_d = {dd}$ kg/m$^3$ -- mean stress vs time")
        ax.grid(True, alpha=0.4)
        ax.legend(loc="best", fontsize=10, framealpha=0.85)
        fig.tight_layout()
        fig.savefig(HERE / f"dsm_ramp_dd{dd}_stress_vs_time.png", dpi=140)
        plt.close(fig)

    print(f"[plot] dd{dd}: PNGs written")


# ----------------------------------------------------------------------------- main
def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--skip-run", action="store_true", help="skip OGS sim, just post-process")
    ap.add_argument("--only-gen", action="store_true", help="only regenerate PRJ files")
    ap.add_argument("--densities", nargs="*", type=int, default=list(DENSITIES))
    args = ap.parse_args(argv)

    print(f"[setup] working dir: {HERE}")
    if not OGS_BIN.exists() and not args.only_gen and not args.skip_run:
        print(f"  ERROR: OGS binary not found at {OGS_BIN}", file=sys.stderr)
        return 1

    for dd in args.densities:
        prj = generate_prj(dd)
        print(f"[gen] wrote {prj.name}")

    if args.only_gen:
        return 0

    if not args.skip_run:
        ok_all = True
        for dd in args.densities:
            prj = HERE / f"dsm_ramp_dd{dd}.prj"
            if not run_one(prj):
                ok_all = False
        if not ok_all:
            print("[run] one or more sims failed; continuing to post-process what is available", file=sys.stderr)

    for dd in args.densities:
        try:
            extract_series(dd)
            plot_density(dd)
        except Exception as exc:
            print(f"[post] dd{dd}: skipped due to {exc!r}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
