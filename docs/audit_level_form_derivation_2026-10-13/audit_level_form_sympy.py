"""Symbolic checks of the appendix 'The level form: what is derived, what is
posited, what is missing' (audit 2026-10-13). Run: python3 audit_level_form_sympy.py
Every check prints OK (difference simplifies to 0) or the non-zero residual."""

import sympy as sp

nl, eps, e, nS, b, K, s0, L0, rho = sp.symbols(
    "n_l eps_v e n_S b K_sk sigma_0 L0 rho_lR", real=True
)
Pi = sp.Function("Pi")(nl)
dPi = sp.diff(Pi, nl)
phim = nS * nl  # (1 - phi_M) n_l
Z = 1 + b * phim


def check(name, expr):
    r = sp.simplify(sp.expand(expr))
    print(f"{'OK  ' if r == 0 else 'FAIL'} {name}" + ("" if r == 0 else f"\n     residual = {r}"))
    return r == 0


# BGR_dv_film: integrate the posited stress at frozen film (nS, n_l, Pi fixed)
P = sp.Symbol("P")  # Pi frozen
Psi_film_int = -sp.integrate(nS * nl * (P + b * K * e), (e, 0, eps))
Psi_film_txt = -nS * nl * (P * eps + sp.Rational(1, 2) * b * K * eps**2)
check("film energy: integral = closed text", Psi_film_int - Psi_film_txt)
check("film energy: dPsi/deps = -(1-phiM) n_l (Pi + b K eps) (eq BGR_dv_sw on drained line)",
      sp.diff(Psi_film_txt, eps) + nS * nl * (P + b * K * eps))

# BGR_dv_partner: n_l derivative / ((1-phiM) rho), Pi = Pi(n_l)
Psi_film = Psi_film_txt.subs(P, Pi)
mu_mech = sp.diff(Psi_film, nl) / (nS * rho)
mu_mech_txt = -((Pi + nl * dPi) * eps + sp.Rational(1, 2) * b * K * eps**2) / rho
check("partner mu_mech (BGR_dv_partner)", mu_mech - mu_mech_txt)
# code PotentialExchange.h l.584-587
A = Pi + nl * dPi
mu_code = -(A * eps + sp.Rational(1, 2) * b * K * eps**2) / rho
check("partner = code PotentialExchange.h l.584-587", mu_mech_txt - mu_code)

# BGR_dv_ident
sig_film = sp.diff(Psi_film, eps)
lhs = sp.diff(sig_film, nl)
mid = nS * rho * sp.diff(mu_mech_txt, eps)
rhs = -nS * (Pi + nl * dPi + b * K * eps)
check("Maxwell ident: dsigma_sw/dn_l = (1-phiM) rho dmu/deps", lhs - mid)
check("Maxwell ident: = -(1-phiM)(Pi + n_l Pi' + b K eps)", lhs - rhs)

# Code closed form (FEM l.3581-3585): L = -nS nl (Pi + b m), m = m_hat + s,
# s = s_prev + L - L_prev
L, mh, sp_, Lp = sp.symbols("L m_hat s_prev L_prev")
sol = sp.solve(sp.Eq(L, -phim * (Pi + b * (mh + sp_ + L - Lp))), L)[0]
code = (-phim * Pi - phim * b * (mh + sp_ - Lp)) / (1 + phim * b)
check("code closed form L(1+c) = -nS nl Pi - c(m_hat + s_prev - L_prev)", sol - code)

# BGR_dv_closed: telescoped, s_prev - L_prev = -L0, m_hat = sigma_0 + K eps
s = sp.Symbol("s")
sol_s = sp.solve(sp.Eq(s + L0, -phim * (Pi + b * (s0 + K * eps + s))), s)[0]
closed_txt = -(phim * (Pi + b * (s0 + K * eps)) + L0) / (1 + b * phim)
check("closed form BGR_dv_closed from L = -phim[Pi + b m], m = s0+K eps+s, s = L-L0", sol_s - closed_txt)
check("code closed form with m_hat = s0 + K eps, s_prev - L_prev = -L0 gives BGR_dv_closed",
      (code.subs({mh: s0 + K * eps, sp_: Lp - L0}) - L0) - closed_txt)

# BGR_dv_energy
Psi_cpl = -(phim / Z) * ((Pi + b * (s0 - L0)) * eps + sp.Rational(1, 2) * b * K * eps**2)
check("dPsi_cpl/deps - L0 = sigma_sw (closed form)", sp.diff(Psi_cpl, eps) - L0 - closed_txt)

# BGR_dv_need
mu_need = sp.diff(Psi_cpl, nl) / (nS * rho)
mu_need_txt = -(1 / rho) * (
    ((Pi + b * (s0 - L0)) * eps + sp.Rational(1, 2) * b * K * eps**2) / Z**2
    + nl * dPi * eps / Z
)
check("mu_need (BGR_dv_need)", mu_need - mu_need_txt)
# branch-tip docs/dsm_z_weight_design_2026-10.md l.35 (rho mu_m)
mu_tip = (-((Pi + b * (s0 - L0)) * eps + sp.Rational(1, 2) * b * K * eps**2) / Z**2
          - nl * dPi * eps / Z) / rho
check("mu_need == branch-tip partner (design doc l.35 / sympy test l.33)", mu_need_txt - mu_tip)
check("Maxwell for Psi_cpl: dsigma_sw/dn_l = (1-phiM) rho dmu_need/deps",
      sp.diff(closed_txt, nl) - nS * rho * sp.diff(mu_need_txt, eps))

# BGR_dv_miss
miss = sp.diff(closed_txt, nl) - nS * rho * sp.diff(mu_mech_txt, eps)
miss_txt = nS * (Pi + nl * dPi + b * K * eps
                 - (Pi + b * (s0 - L0) + b * K * eps) / Z**2 - nl * dPi / Z)
check("miss BGR_dv_miss", miss - miss_txt)

# vanish at eps = 0
check("mu_mech(eps=0) = 0", mu_mech_txt.subs(eps, 0))
check("mu_need(eps=0) = 0", mu_need_txt.subs(eps, 0))
# miss at b = 0 (Z = 1): should be 0
check("miss vanishes at b = 0", miss_txt.subs(b, 0))
print("miss at eps = 0:", sp.simplify(miss_txt.subs(eps, 0)))
