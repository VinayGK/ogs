#!/usr/bin/env python3
"""Symbolic check of the Z-weighted strain partner of the level-form swelling stress (2026-10-13).

Stress (level form, linear skeleton, closed form of RichardsMechanicsFEM-impl.h `L (1 + c) = ...`):
    sigma_sw = -(phi_m/Z) [Pi + b (s0 - L0) + b K eps] - L0,   Z = 1 + b phi_m,  phi_m = nS n_l.
Energy by integration over eps at fixed n_l; partner rho mu = (1/nS) dPsi/dn_l (nS, rho, K, s0, L0 held).
Checks: (1) dPsi/deps = sigma_sw + L0 (the constant L0 does not come from Psi), (2) Maxwell identity
d sigma_sw/d n_l = nS * d(rho mu)/d eps, (3) the partner equals the closed form coded in
computeZWeightedMechanicalMicroPotential, (4) at Z = 1 (b = 0) it is the drained-line partner,
(5) the difference Z-partner minus drained-line partner, as an explicit expression (sign re-derived here).
Run: python3 Tests/Python/dsm_z_weight_partner_sympy.py  -> prints PASS/FAIL per check, exit 1 on FAIL."""
import sys
import sympy as sp

nl, eps, nS, b, K, s0, L0, rho = sp.symbols("n_l eps_v n_S b K s0 L0 rho", real=True)
Pi = sp.Function("Pi")(nl)
phim = nS * nl
Z = 1 + b * phim
A2 = Pi + b * (s0 - L0)
sigma_sw = -(phim / Z) * (A2 + b * K * eps) - L0
Psi = -(phim / Z) * (A2 * eps + sp.Rational(1, 2) * b * K * eps**2)
rho_mu = sp.diff(Psi, nl) / nS

ok = True
def check(name, expr):
    global ok
    z = sp.simplify(expr)
    print(("PASS " if z == 0 else "FAIL ") + name + ("" if z == 0 else f": residual {z}"))
    ok = ok and z == 0

check("(1) dPsi/deps = sigma_sw + L0", sp.diff(Psi, eps) - (sigma_sw + L0))
check("(2) Maxwell: d sigma_sw/d n_l = nS d(rho mu)/d eps", sp.diff(sigma_sw, nl) - nS * sp.diff(rho_mu, eps))
coded = -(A2 * eps + sp.Rational(1, 2) * b * K * eps**2) / Z**2 - nl * sp.diff(Pi, nl) * eps / Z
check("(3) partner = coded closed form", rho_mu - coded)
drained = -((Pi + nl * sp.diff(Pi, nl)) * eps + sp.Rational(1, 2) * b * K * eps**2)
check("(4) b = 0 -> drained-line partner", (rho_mu - drained).subs(b, 0))
# coded derivative in n_l (total; Pi, Pi' move; nS held) against sympy
Pip, Pipp = sp.symbols("Pi1 Pi2")
dZ = b * nS
Q = A2 * eps + sp.Rational(1, 2) * b * K * eps**2
coded_dnl = -(-sp.diff(Pi, nl) * eps / Z**2 + 2 * Q * dZ / Z**3 - (sp.diff(Pi, nl) + nl * sp.diff(Pi, nl, 2)) * eps / Z + nl * sp.diff(Pi, nl) * eps * dZ / Z**2) * (-1)
check("(3b) coded d(rho mu)/d n_l", sp.diff(coded, nl) - coded_dnl)
diff = sp.simplify(rho_mu - drained)
print("(5) Z-partner minus drained-line partner, rho (mu_Z - mu_code) =")
print("   ", sp.collect(sp.expand(diff), eps))
print("    at Z -> 1 this vanishes when s0 = L0; the Pi eps term carries the factor (1 - 1/Z^2) with the SIGN of -Pi eps:")
print("   ", sp.simplify(diff.subs({s0: L0, K: 0}).subs(sp.diff(Pi, nl), 0)))
sys.exit(0 if ok else 1)
