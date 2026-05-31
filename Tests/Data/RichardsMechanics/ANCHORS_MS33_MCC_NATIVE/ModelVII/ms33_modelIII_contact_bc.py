# SPDX-License-Identifier: BSD-3-Clause
# Model III technological-gap CONTACT surrogate via a solution-dependent
# Python Neumann (traction) BC on the bentonite outer radial boundary.
#
# Physics: the bentonite cylinder (r=25 mm) swells radially into a 2 mm gap.
# While u_r < gap the boundary is free (gap open). Once u_r >= gap the rigid
# host pushes back with a one-sided penalty reaction traction (inward).
#
# REGULARISED (C1) penalty. A hard kink t_n = -K_PEN*max(0, u_r-gap) makes the
# tangent jump 0 -> -K_PEN at first contact; under the adaptive time stepping the
# boundary crosses the gap in one large step and the abrupt tangent makes the
# global Newton stall (the |dx| floor rises until no tolerance is reachable;
# observed failures at ~16 d). Instead the contact stiffness ramps LINEARLY from
# 0 to K_PEN over a thin band DELTA, so the reaction is C1 and engagement is
# gentle:
#     g = u_r - gap
#     g <= 0     : t_n = 0,                       dt/du_r = 0          (gap open)
#     0 < g < DELTA : t_n = -0.5*K_PEN*g^2/DELTA, dt/du_r = -K_PEN*g/DELTA
#     g >= DELTA : t_n = -K_PEN*(g - DELTA/2),    dt/du_r = -K_PEN     (full penalty)
# Value and slope are continuous at g=0 and g=DELTA; the open gap carries exactly
# zero spurious force. cf. OGS Tests/Data/Mechanics/Linear/PythonHertzContact.
try:
    import ogs.callbacks as OpenGeoSys
except ModuleNotFoundError:
    import OpenGeoSys

GAP = 0.002        # technological gap [m] (2 mm)
K_PEN = 5.0e10     # full penalty stiffness [Pa/m] (t=K_PEN*penetration at full contact)
DELTA = 5.0e-5     # engagement smoothing band [m] (0.05 mm): stiffness ramps 0 -> K_PEN

_DBG = {"n": 0, "eng": 0}


class RadialContactBC(OpenGeoSys.BoundaryCondition):
    def getFlux(self, t, coords, primary_vars):
        # RM primary_vars has 3 components; verify order on first calls.
        if _DBG["n"] < 4:
            print(f"[contactBC] t={t:.3g} coords={coords} primary_vars={list(primary_vars)}")
            _DBG["n"] += 1
        # RM primary_vars order is [pressure, u_x(=u_r, radial), u_y] (verified via debug).
        u_r = primary_vars[1]          # radial displacement (axisymmetric x-component)
        g = u_r - GAP
        if g <= 0.0:
            return (True, 0.0, [0.0, 0.0, 0.0])               # gap open -> free
        if g < DELTA:
            flux = -0.5 * K_PEN * g * g / DELTA               # quadratic ramp
            dflux_dur = -K_PEN * g / DELTA
        else:
            flux = -K_PEN * (g - 0.5 * DELTA)                 # full linear penalty
            dflux_dur = -K_PEN
        if _DBG["eng"] < 12:
            print(f"[contactBC ENGAGE] t={t:.4g} u_r={u_r * 1e3:.4f}mm "
                  f"penetration={g * 1e3:.4f}mm t_n={flux:.4g}Pa")
            _DBG["eng"] += 1
        return (True, flux, [0.0, dflux_dur, 0.0])            # d(flux)/d(p, u_x, u_y)


contact_bc = RadialContactBC()
