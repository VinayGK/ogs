#!/usr/bin/env python3
"""Model IV figures for the PER-MATERIAL run (clay K=83377, pellet K=6600), both suites.
Adapted from plot_ms33_fastex.py model_IV (same CIMNE house style + spec axes + BGR tag),
reading permaterial_IV/<LE|MCC>/out and writing to permaterial_IV/figs/.
Outputs: ms33_<LE|MCC>_modelIV_fig1_stress.pdf, fig2_density.pdf, fig3_ps_path.pdf .
Run: python3 plot_permaterial_IV.py
"""
import glob, re, os
import numpy as np
import meshio
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

PM = "/Users/vinaykumar/git/ogs-worktrees/dsm_native_hierarchical_wt/Tests/Data/RichardsMechanics/ANCHORS_MS33_FASTEX_2026-06-02/permaterial_IV"
OUT = PM + "/figs"
QC = "/tmp/ms33png_permaterial_IV"
os.makedirs(OUT, exist_ok=True)
os.makedirs(QC, exist_ok=True)

GRID = "#D9D9D9"
MAT = {"clay": ("bentonite (clay)", "#000000"), "pellet": ("pellets", "#C00000")}


def series(suite):
    return sorted(glob.glob(f"{PM}/{suite}/out/*_ts_*.vtu"),
                  key=lambda f: int(re.search(r"_ts_(\d+)_", f).group(1)))


def tdays(f):
    return float(re.search(r"_t_([0-9.]+)\.vtu", f).group(1)) / 86400.0


def meanp(s):
    s = np.atleast_2d(s)
    return -s[:, :3].sum(axis=1) / 3.0


def style(ax, xl, yl, xlim=None, ylim=None, xlog=False, ylog=False):
    if xlog:
        ax.set_xscale("log")
    if ylog:
        ax.set_yscale("log")
    if xlim:
        ax.set_xlim(*xlim)
    if ylim:
        ax.set_ylim(*ylim)
    ax.set_xlabel(xl, fontsize=11, fontweight="bold")
    ax.set_ylabel(yl, fontsize=11, fontweight="bold")
    ax.grid(True, which="major", color=GRID, lw=0.8, zorder=0)
    ax.grid(True, which="minor", color=GRID, lw=0.4, zorder=0, alpha=0.7)
    ax.tick_params(labelsize=9)
    for sp in ax.spines.values():
        sp.set_color("black"); sp.set_linewidth(0.9)
    ax.set_axisbelow(True)


def tag(ax, suite):
    ax.text(0.5, 1.02, f"BGR · {suite}", transform=ax.transAxes, ha="center", va="bottom",
            fontsize=11, fontweight="bold",
            bbox=dict(boxstyle="square,pad=0.3", fc="white", ec="black", lw=1.1))


def arrow(ax, x, y, frac=0.5, color="k"):
    x, y = np.asarray(x), np.asarray(y)
    n = len(x)
    if n < 2:
        return
    i = min(n - 1, max(1, int(frac * n)))
    ax.annotate("", xy=(x[i], y[i]), xytext=(x[i - 1], y[i - 1]),
                arrowprops=dict(arrowstyle="-|>", color=color, lw=1.5), zorder=6)


def boxleg(ax, handles=None, loc="upper right"):
    ax.legend(handles=handles, fontsize=8.5, loc=loc, framealpha=1.0, edgecolor="#BFBFBF", handlelength=1.8)


def save(fig, name):
    fig.tight_layout()
    fig.savefig(f"{OUT}/{name}", dpi=200)
    fig.savefig(f"{QC}/{name.replace('.pdf', '.png')}", dpi=130)
    plt.close(fig)
    print("  wrote", name)


def model_IV(suite):
    fs = series(suite)
    if not fs:
        print(f"  {suite} IV: NO DATA"); return
    pts = meshio.read(fs[0]).points
    clay = pts[:, 1] >= 0.035
    pel = pts[:, 1] < 0.035
    t, cp, pp, crho, prho, cs, ps = [], [], [], [], [], [], []
    for f in fs:
        pd = meshio.read(f).point_data
        sig = meanp(np.atleast_2d(pd["sigma"]))
        dd = np.ravel(pd["dry_density_solid"])
        pr = -np.ravel(pd["pressure"]) / 1e6
        t.append(tdays(f))
        cp.append(sig[clay].mean() / 1e6); pp.append(sig[pel].mean() / 1e6)
        crho.append(dd[clay].mean()); prho.append(dd[pel].mean())
        cs.append(pr[clay].mean()); ps.append(pr[pel].mean())
    t = np.array(t); m = t <= 200.0
    t = t[m]; cp = np.array(cp)[m]; pp = np.array(pp)[m]
    crho = np.array(crho)[m]; prho = np.array(prho)[m]
    cs = np.array(cs)[m]; ps = np.array(ps)[m]

    # fig1: stress vs time
    fig, ax = plt.subplots(figsize=(5.3, 4.0))
    ax.plot(t, cp, "-", color=MAT["clay"][1], lw=2.0, label="bentonite (clay)", zorder=3)
    ax.plot(t, pp, "-", color=MAT["pellet"][1], lw=2.0, label="pellets", zorder=3)
    style(ax, "time (days)", r"mean stress $\bar p$ (MPa)", (0, 200), (0, 15))
    boxleg(ax, loc="lower right"); tag(ax, suite)
    save(fig, f"ms33_{suite}_modelIV_fig1_stress.pdf")

    # fig2: dry density vs time
    fig, ax = plt.subplots(figsize=(5.3, 4.0))
    ax.plot(t, crho, "-", color=MAT["clay"][1], lw=2.0, label="bentonite (clay)", zorder=3)
    ax.plot(t, prho, "-", color=MAT["pellet"][1], lw=2.0, label="pellets", zorder=3)
    style(ax, "time (days)", r"dry density (kg/m$^3$)", (0, 200), (900, 1700))
    boxleg(ax, loc="center right"); tag(ax, suite)
    save(fig, f"ms33_{suite}_modelIV_fig2_density.pdf")

    # fig3: p-s path
    fig, ax = plt.subplots(figsize=(5.3, 4.2))
    for vp, vs, col, lab in [(cp, cs, MAT["clay"][1], "bentonite (clay)"),
                             (pp, ps, MAT["pellet"][1], "pellets")]:
        spc = np.clip(vs, 0.1, None)
        ax.plot(vp, spc, "-", color=col, lw=2.0, label=lab, zorder=3); arrow(ax, vp, spc, 0.4, col)
    style(ax, r"mean stress $\bar p$ (MPa)", r"suction $s$ (MPa)", (0, 15), (0.1, 100), ylog=True)
    boxleg(ax); tag(ax, suite)
    save(fig, f"ms33_{suite}_modelIV_fig3_ps_path.pdf")
    print(f"  {suite} IV: clay p={cp[-1]:.2f}/pellet p={pp[-1]:.2f}; rho {crho[-1]:.0f}/{prho[-1]:.0f}")


if __name__ == "__main__":
    for suite in ("LE", "MCC"):
        print(f"== {suite} ==")
        model_IV(suite)
    print("== done -> figs/ ==")
