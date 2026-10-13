#!/usr/bin/env python3
"""Generate the one-element closed-loop decks of the DSM loop test (2026-10).

Parent deck: ../ANCHORS_MS33_ModelI/ms33_modelI_dd1600.prj (read, never changed).
Changes against the parent (all else, materials and numerics, as the parent):
  * <swelling_stress_form>level</swelling_stress_form> under <potential_exchange>
    (scalar K kept), optionally <micro_potential_partner_z_weight>true</...>;
  * the top face (y = 0.01 m) is driven in y by the curve `strain` scaled by
    1e-4 m (= +1 % of the 0.01 m height; left/right/bottom stay fixed, so the
    volumetric strain is the axial strain);
  * the suction is driven by the curve `suction` (pressure on all four nodes);
  * fixed output times every `hold_dt` d in holds, every `ramp_dt` d in ramps;
  * micro_liquid_density and epsilon are written in addition;
  * no <test_definition> unless `--test` (the registered deck).

Schedule (days): a dry hold of H0 d at 100 MPa suction (the loop starts at its
end), then four legs, each a 200 d linear ramp followed by a hold of H d:
  cycle A: wet (p_L -1e8 -> 0), strain (0 -> 1 %), dry (0 -> -1e8), unstrain
  cycle B: strain, wet, unstrain, dry
Usage: make_loop_decks.py OUTDIR [--H0 d] [--H d] [--ramp-dt d] [--hold-dt d]
"""

import argparse
import re
from pathlib import Path

DAY = 86400.0
HERE = Path(__file__).resolve().parent
PARENT = HERE.parent / "ANCHORS_MS33_ModelI" / "ms33_modelI_dd1600.prj"

LEGS = {
    "A": ["wet", "strain", "dry", "unstrain"],
    "B": ["strain", "wet", "unstrain", "dry"],
}


def schedule(cycle, H0, H, ramp=200.0):
    """Return (t, p_L, strain) knots in days and the ramp intervals."""
    p, e = -1e8, 0.0
    t = [0.0, H0]
    P = [p, p]
    E = [e, e]
    ramps = []
    for leg in LEGS[cycle]:
        t0 = t[-1]
        if leg == "wet":
            p = 0.0
        elif leg == "dry":
            p = -1e8
        elif leg == "strain":
            e = 1.0
        elif leg == "unstrain":
            e = 0.0
        t += [t0 + ramp, t0 + ramp + H]
        P += [p, p]
        E += [e, e]
        ramps.append((t0, t0 + ramp))
    return t, P, E, ramps


def output_times(t, ramps, ramp_dt, hold_dt):
    """Output times in days: every ramp_dt inside ramps, hold_dt elsewhere."""
    out = set()

    def fill(a, b, d):
        n = int(round((b - a) / d))
        for i in range(n + 1):
            out.add(round(a + i * (b - a) / n, 9))

    edges = sorted(set(t))
    for a, b in zip(edges[:-1], edges[1:]):
        is_ramp = any(abs(a - r0) < 1e-9 and abs(b - r1) < 1e-9 for r0, r1 in ramps)
        fill(a, b, ramp_dt if is_ramp else hold_dt)
    out.discard(0.0)
    return sorted(out)


def fmt(x):
    return repr(float(x))


def make(cycle, z, H0, H, ramp_dt, hold_dt, prefix, test_frame=None):
    s = PARENT.read_text(encoding="utf-8")
    s = re.sub(r"<!--.*?-->", "", s, flags=re.S)
    s = re.sub(r"\n\s*\n+", "\n", s)
    t, P, E, ramps = schedule(cycle, H0, H)

    extra = "                <swelling_stress_form>level</swelling_stress_form>\n"
    if z:
        extra += (
            "                <micro_potential_partner_z_weight>true"
            "</micro_potential_partner_z_weight>\n"
        )
    s = s.replace("            </potential_exchange>", extra + "            </potential_exchange>", 1)

    top_old = (
        "<boundary_condition><mesh>square_1e-2_quad_1e0_top</mesh><type>Dirichlet</type>"
        "<component>1</component><parameter>zero</parameter></boundary_condition>"
    )
    assert top_old in s
    s = s.replace(top_old, top_old.replace("<parameter>zero</parameter>", "<parameter>u_top</parameter>"))

    par_old = (
        "<parameter><name>pressure_bc</name><type>CurveScaled</type>"
        "<curve>suction_ramp</curve><parameter>p_bc_scale</parameter></parameter>"
    )
    assert par_old in s
    s = s.replace(
        par_old,
        "<parameter><name>pressure_bc</name><type>CurveScaled</type>"
        "<curve>suction</curve><parameter>p_bc_scale</parameter></parameter>\n"
        "        <parameter><name>u_top_scale</name><type>Constant</type><value>1e-4</value></parameter>\n"
        "        <parameter><name>u_top</name><type>CurveScaled</type>"
        "<curve>strain</curve><parameter>u_top_scale</parameter></parameter>",
    )
    # p_bc_scale stays 1: the suction curve carries p_L in Pa.
    curves = (
        "    <curves>\n"
        "        <curve>\n            <name>suction</name>\n"
        f"            <coords>{' '.join(fmt(x * DAY) for x in t)}</coords>\n"
        f"            <values>{' '.join(fmt(x) for x in P)}</values>\n"
        "        </curve>\n"
        "        <curve>\n            <name>strain</name>\n"
        f"            <coords>{' '.join(fmt(x * DAY) for x in t)}</coords>\n"
        f"            <values>{' '.join(fmt(x) for x in E)}</values>\n"
        "        </curve>\n"
        "    </curves>"
    )
    s = re.sub(r"    <curves>.*?</curves>", curves, s, flags=re.S)
    t_end = t[-1] * DAY
    s = re.sub(r"<t_end>[^<]*</t_end>", f"<t_end>{fmt(t_end)}</t_end>", s)
    s = re.sub(r"<prefix>[^<]*</prefix>", f"<prefix>{prefix}</prefix>", s)
    outs = output_times(t, ramps, ramp_dt, hold_dt)
    s = re.sub(
        r"<fixed_output_times>.*?</fixed_output_times>",
        "<fixed_output_times>\n                "
        + "\n                ".join(fmt(x * DAY) for x in outs)
        + "\n            </fixed_output_times>",
        s,
        flags=re.S,
    )
    s = s.replace(
        "                <variable>micro_exchange_source</variable>\n",
        "                <variable>micro_exchange_source</variable>\n"
        "                <variable>micro_liquid_density</variable>\n"
        "                <variable>epsilon</variable>\n",
    )
    s = s.replace(
        '                <secondary_variable name="micro_exchange_source"/>\n',
        '                <secondary_variable name="micro_exchange_source"/>\n'
        '                <secondary_variable name="micro_liquid_density"/>\n',
    )
    s = s.replace("../square_1e-2_quad_1e0", "square_1e-2_quad_1e0")
    td = re.search(r"    <test_definition>.*?</test_definition>\n", s, flags=re.S)
    s = s[: td.start()] + s[td.end() :]
    if test_frame is not None:
        fields = [
            ("displacement", "1e-9", "1e-7"),
            ("porosity", "1e-8", "1e-7"),
            ("micro_water_content", "1e-8", "1e-7"),
            ("sigma", "1e3", "1e-2"),
            ("swelling_stress", "1e3", "1e-2"),
            ("pressure", "1e3", "1e-2"),
        ]
        rx = re.escape(prefix) + r"_ts_.*_t_" + f"{t_end:.6f}".replace(".", r"\.") + r"\.vtu"
        td = "    <test_definition>\n" + "".join(
            f"        <vtkdiff><regex>{rx}</regex><field>{f}</field>"
            f"<absolute_tolerance>{a}</absolute_tolerance>"
            f"<relative_tolerance>{r}</relative_tolerance></vtkdiff>\n"
            for f, a, r in fields
        ) + "    </test_definition>\n"
        s = s.replace("</OpenGeoSysProject>", td + "</OpenGeoSysProject>")
    header = (
        "<!--\nDSM one-element CLOSED LOOP, cycle {c}, level-form swelling stress, partner: {p}.\n"
        "Generated by make_loop_decks.py (this directory) from ../ANCHORS_MS33_ModelI/ms33_modelI_dd1600.prj;\n"
        "dry hold H0 = {h0} d, ramps 200 d, holds H = {h} d, output every {rd} d in ramps and {hd} d in holds.\n"
        "See docs/dsm_closed_loop_2026-10.md.\n-->\n"
    ).format(c=cycle, p="Z-weighted" if z else "drained-line (default)", h0=H0, h=H, rd=ramp_dt, hd=hold_dt)
    s = s.replace("<OpenGeoSysProject>", header + "<OpenGeoSysProject>", 1)
    return s, outs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("outdir")
    ap.add_argument("--H0", type=float, default=400.0)
    ap.add_argument("--H", type=float, default=400.0)
    ap.add_argument("--ramp-dt", type=float, default=0.2)
    ap.add_argument("--hold-dt", type=float, default=2.0)
    ap.add_argument("--tag", default="")
    ap.add_argument("--test", action="store_true")
    a = ap.parse_args()
    od = Path(a.outdir)
    od.mkdir(parents=True, exist_ok=True)
    for cycle in "AB":
        for z in (False, True):
            prefix = f"dsm_loop_{cycle}_{'zon' if z else 'off'}{a.tag}"
            s, outs = make(cycle, z, a.H0, a.H, a.ramp_dt, a.hold_dt, prefix,
                           test_frame=True if a.test else None)
            (od / f"{prefix}.prj").write_text(s, encoding="utf-8")
            print(prefix, len(outs), "output times, t_end", outs[-1], "d")


if __name__ == "__main__":
    main()
