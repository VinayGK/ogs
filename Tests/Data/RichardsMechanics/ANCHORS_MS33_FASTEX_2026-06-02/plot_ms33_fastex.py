#!/usr/bin/env python3
"""Standardized MS33 figures from the FASTEX runs, for BOTH suites (LE + MCC).

Same house style + spec axes as plot_ms33_standard.py (CIMNE convention + Required_figures.tex),
but reads the recalibrated FASTEX runs (krel=0.1, alpha=3e-13, K->Dixon) and names every output
by suite:  ms33_<LE|MCC>_modelX_figN_*.pdf , tag "BGR . <LE|MCC>".

Gaps are handled honestly (MCC dd1800 + MCC VII have no data -> skipped/annotated).
Run: python3 plot_ms33_fastex.py
"""
import glob
import re
import os
import numpy as np
import meshio
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

FX = "/Users/vinaykumar/git/ogs-worktrees/dsm_native_hierarchical_wt/Tests/Data/RichardsMechanics/ANCHORS_MS33_FASTEX_2026-06-02"
OUT = FX + "/figures"
QC = "/tmp/ms33png_fastex"
os.makedirs(OUT, exist_ok=True)
os.makedirs(QC, exist_ok=True)

GRID = "#D9D9D9"
DD_COL = {1400: ("1.4", "#E0A100"), 1600: ("1.6", "#000000"), 1800: ("1.8", "#2E74B5")}
COMP = {"mean": (r"$\bar p$ (mean)", "#000000"), "v": (r"$\sigma_v$ (axial)", "#C00000"),
        "h": (r"$\sigma_h$ (radial)", "#2E74B5")}
MAT = {"clay": ("bentonite (clay)", "#000000"), "pellet": ("pellets", "#C00000")}
DIXON_T = {1400: 4.922, 1600: 14.161, 1800: 40.860}


def series(suite, key):
    # LE dd1800: use the power-law-k_rel variant (completes 200 d; constant-k_rel corner-crashed)
    d = "I_dd1800_pl" if (suite == "LE" and key == "I_dd1800") else key
    return sorted(glob.glob(f"{FX}/{suite}/{d}/out/*_ts_*.vtu"),
                  key=lambda f: int(re.search(r"_ts_(\d+)_", f).group(1)))


def tdays(f):
    return float(re.search(r"_t_([0-9.]+)\.vtu", f).group(1)) / 86400.0


def meanp(s):
    s = np.atleast_2d(s)
    return -s[:, :3].sum(axis=1) / 3.0


def nidx(p, r, z):
    return int(np.argmin((p[:, 0] - r) ** 2 + (p[:, 1] - z) ** 2))


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


# ============================================================== Model I
def model_I(suite):
    figp, axp = plt.subplots(figsize=(5.3, 4.3))
    figk, axk = plt.subplots(figsize=(5.3, 4.3))
    got = False
    for dd, (lab, col) in DD_COL.items():
        fs = series(suite, f"I_dd{dd}")
        if len(fs) < 3:
            print(f"  {suite} I dd{dd}: no real run ({len(fs)} steps) -- skip"); continue
        got = True
        s, p, k = [], [], []
        for f in fs:
            pd = meshio.read(f).point_data
            s.append(-np.mean(pd["pressure"]) / 1e6)
            p.append(np.mean(meanp(pd["sigma"])) / 1e6)
            k.append(np.mean(np.atleast_2d(pd["intrinsic_permeability"])[:, 0]))
        s, p, k = np.array(s), np.array(p), np.array(k)
        sp = np.clip(s, 0.1, None)
        axp.plot(p, sp, "-", color=col, lw=2.0, zorder=3); arrow(axp, p, sp, 0.4, col)
        axk.plot(k, sp, "-", color=col, lw=2.0, zorder=3)
        axk.plot(k[0], sp[0], "o", color=col, ms=7, zorder=4); arrow(axk, k, sp, 0.4, col)
        print(f"  {suite} I dd{dd}: p_sw={p[-1]:.2f} k={k[-1]:.2e}")
    if not got:
        plt.close(figp); plt.close(figk); return
    for dd in DD_COL:
        axp.plot(DIXON_T[dd], 0.105, "s", mfc="#9A9A9A", mec="#5A5A5A", ms=9, mew=0.7, clip_on=False, zorder=5)
    style(axp, r"mean effective stress $p'$ (MPa)", r"suction $s$ (MPa)", (0, 60), (0.1, 100), ylog=True)
    h = [Line2D([], [], color=DD_COL[d][1], lw=2.0, label=fr"$\rho_d$={DD_COL[d][0]} g/cm$^3$") for d in DD_COL]
    h.append(Patch(fc="#9A9A9A", ec="#5A5A5A", label="Dixon (2023)"))
    boxleg(axp, h); tag(axp, suite); save(figp, f"ms33_{suite}_modelI_fig1_ps_path.pdf")
    style(axk, r"intrinsic permeability $k$ (m$^2$)", r"suction $s$ (MPa)", (1e-23, 1e-19), (0.1, 100), xlog=True, ylog=True)
    boxleg(axk, h[:-1]); tag(axk, suite); save(figk, f"ms33_{suite}_modelI_fig2_k_suction.pdf")


# ============================================================== Model III
def model_III(suite):
    fs = series(suite, "III")
    if not fs:
        print(f"  {suite} III: NO DATA"); return
    pts = meshio.read(fs[0]).points
    c = nidx(pts, 0.0125, 0.035); gi = nidx(pts, 0.025, 0.035)
    t, sm, sv, sh, su, ap = [], [], [], [], [], []
    for f in fs:
        pd = meshio.read(f).point_data; sig = np.atleast_2d(pd["sigma"])
        t.append(tdays(f)); sm.append(meanp(sig)[c] / 1e6); sv.append(-sig[c, 1] / 1e6); sh.append(-sig[c, 0] / 1e6)
        su.append(-np.ravel(pd["pressure"])[c] / 1e6)
        ap.append(2.0 - np.atleast_2d(pd["displacement"])[gi, 0] * 1e3)
    t = np.array(t); m = t <= 200.0
    t, sm, sv, sh, su, ap = (t[m], np.array(sm)[m], np.array(sv)[m], np.array(sh)[m], np.array(su)[m], np.array(ap)[m])
    fig, ax = plt.subplots(figsize=(5.3, 4.0)); ax.plot(t, ap, "-", color="#000000", lw=2.0, zorder=3)
    style(ax, "time (days)", "gap aperture (mm)", (0, 200), (0.0, 5.0)); tag(ax, suite); save(fig, f"ms33_{suite}_modelIII_fig1_gap.pdf")
    fig, ax = plt.subplots(figsize=(5.3, 4.0))
    ax.plot(t, sm, "-", color=COMP["mean"][1], lw=2.4, label=COMP["mean"][0], zorder=4)
    ax.plot(t, sv, "-", color=COMP["v"][1], lw=1.8, label=COMP["v"][0], zorder=3)
    ax.plot(t, sh, "-", color=COMP["h"][1], lw=1.8, label=COMP["h"][0], zorder=3)
    style(ax, "time (days)", "stress (MPa)", (0, 200), (0, 15)); boxleg(ax, loc="lower right"); tag(ax, suite)
    save(fig, f"ms33_{suite}_modelIII_fig2_stress.pdf")
    fig, ax = plt.subplots(figsize=(5.3, 4.2)); spc = np.clip(su, 0.1, None)
    ax.plot(sm, spc, "-", color="#000000", lw=2.0, zorder=3); arrow(ax, sm, spc, 0.4, "#000000")
    style(ax, r"mean stress $\bar p$ (MPa)", r"suction $s$ (MPa)", (0, 15), (0.1, 100), ylog=True); tag(ax, suite)
    save(fig, f"ms33_{suite}_modelIII_fig3_ps_path.pdf")
    print(f"  {suite} III: gap {ap[0]:.2f}->{ap[-1]:.2f} mm; centre p={sm[-1]:.2f} MPa")


# =============================================================== Model IV
def model_IV(suite):
    fs = series(suite, "IV")
    if not fs:
        print(f"  {suite} IV: NO DATA"); return
    pts = meshio.read(fs[0]).points; clay = pts[:, 1] >= 0.035; pel = pts[:, 1] < 0.035
    t, cp, pp, crho, prho, cs, ps = [], [], [], [], [], [], []
    for f in fs:
        pd = meshio.read(f).point_data; sig = meanp(np.atleast_2d(pd["sigma"])); dd = np.ravel(pd["dry_density_solid"])
        pr = -np.ravel(pd["pressure"]) / 1e6
        t.append(tdays(f)); cp.append(sig[clay].mean() / 1e6); pp.append(sig[pel].mean() / 1e6)
        crho.append(dd[clay].mean()); prho.append(dd[pel].mean()); cs.append(pr[clay].mean()); ps.append(pr[pel].mean())
    t = np.array(t); m = t <= 200.0
    t = t[m]; cp = np.array(cp)[m]; pp = np.array(pp)[m]; crho = np.array(crho)[m]; prho = np.array(prho)[m]
    cs = np.array(cs)[m]; ps = np.array(ps)[m]
    fig, ax = plt.subplots(figsize=(5.3, 4.0))
    ax.plot(t, cp, "-", color=MAT["clay"][1], lw=2.0, label="bentonite (clay)", zorder=3)
    ax.plot(t, pp, "-", color=MAT["pellet"][1], lw=2.0, label="pellets", zorder=3)
    style(ax, "time (days)", r"mean stress $\bar p$ (MPa)", (0, 200), (0, 15)); boxleg(ax, loc="lower right"); tag(ax, suite)
    save(fig, f"ms33_{suite}_modelIV_fig1_stress.pdf")
    fig, ax = plt.subplots(figsize=(5.3, 4.0))
    ax.plot(t, crho, "-", color=MAT["clay"][1], lw=2.0, label="bentonite (clay)", zorder=3)
    ax.plot(t, prho, "-", color=MAT["pellet"][1], lw=2.0, label="pellets", zorder=3)
    style(ax, "time (days)", r"dry density (kg/m$^3$)", (0, 200), (900, 1700)); boxleg(ax, loc="center right"); tag(ax, suite)
    save(fig, f"ms33_{suite}_modelIV_fig2_density.pdf")
    fig, ax = plt.subplots(figsize=(5.3, 4.2))
    for vp, vs, col, lab in [(cp, cs, MAT["clay"][1], "bentonite (clay)"), (pp, ps, MAT["pellet"][1], "pellets")]:
        spc = np.clip(vs, 0.1, None); ax.plot(vp, spc, "-", color=col, lw=2.0, label=lab, zorder=3); arrow(ax, vp, spc, 0.4, col)
    style(ax, r"mean stress $\bar p$ (MPa)", r"suction $s$ (MPa)", (0, 15), (0.1, 100), ylog=True); boxleg(ax); tag(ax, suite)
    save(fig, f"ms33_{suite}_modelIV_fig3_ps_path.pdf")
    print(f"  {suite} IV: clay p={cp[-1]:.2f}/pellet p={pp[-1]:.2f}; rho {crho[-1]:.0f}/{prho[-1]:.0f}")


# =============================================================== Model VII
def model_VII(suite):
    fs = series(suite, "VII")
    if len(fs) < 3:
        print(f"  {suite} VII: no real run ({len(fs)} steps) -- skip"); return
    pts = meshio.read(fs[0]).points; c = nidx(pts, 0.0125, 0.035)
    t, e, p, kk = [], [], [], []
    for f in fs:
        pd = meshio.read(f).point_data; phi = float(np.ravel(pd["porosity"])[c])
        t.append(tdays(f)); e.append(phi / (1 - phi)); p.append(meanp(np.atleast_2d(pd["sigma"]))[c] / 1e6)
        kk.append(np.atleast_2d(pd["intrinsic_permeability"])[c, 0])
    t, e, p, kk = np.array(t), np.array(e), np.array(p), np.array(kk)
    fig, ax = plt.subplots(figsize=(5.6, 4.3))
    for lab, col, lo, hi in [("hydration (0--200 d)", "#2E74B5", 0, 200), ("loading (200--220 d)", "#C00000", 200, 220),
                             ("unloading (220--240 d)", "#1A9850", 220, 245)]:
        mk = (t >= lo) & (t <= hi)
        if mk.sum():
            ax.plot(np.clip(p[mk], 0.1, None), e[mk], "-o", color=col, ms=3, lw=1.8, label=lab, zorder=3)
    style(ax, r"mean net stress $\bar p$ (MPa)", "void ratio $e$ (-)", (0.1, 10), (0.4, 1.2), xlog=True)
    boxleg(ax, loc="lower left"); tag(ax, suite); save(fig, f"ms33_{suite}_modelVII_fig1_void_ratio.pdf")
    fig, ax = plt.subplots(figsize=(5.3, 4.0)); ax.plot(t, kk, "-", color="#000000", lw=2.0, zorder=3)
    style(ax, "time (days)", r"intrinsic permeability $k$ (m$^2$)", (0, float(t.max())), (1e-22, 1e-18), ylog=True)
    tag(ax, suite); save(fig, f"ms33_{suite}_modelVII_fig2_k.pdf")
    print(f"  {suite} VII: e {e[0]:.3f}->{e[-1]:.3f}, p_peak={np.nanmax(p):.2f} MPa")


if __name__ == "__main__":
    for suite in ("LE", "MCC"):
        print(f"== {suite} ==")
        model_I(suite); model_III(suite); model_IV(suite); model_VII(suite)
    print("== done -> figures/ ==")
