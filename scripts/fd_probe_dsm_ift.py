"""
FD probe of the analytical IFT formula for the DSM microstate tangent.

What this verifies
------------------
The companion derivation in
  /Users/vinaykumar/tex/cc2024/VK_B35_Pinion_May_2026/dsm_microstate_ift_units.tex
gives two analytical implicit-function-theorem formulas at a converged
microstate (n_l*, rho_lR*):

    d n_l*/d(eps_v) = phi_m(n_l*) * rho_lR* * J22 / det(J)
    d n_l*/d(p)     = dt * alpha_eff * J22 / (rho_macro_safe * det(J))

where J is the analytical 2x2 Jacobian of the residual system

    r1 = phi_m(n_l)*rho_lR - prev_mass - dt*alpha_eff*(mu_LR - mu_lR(n_l, rho_lR))
                            - phi_m(n_l)*rho_lR*eps_v
    r2 = rho_lR - rho_micro_eos(n_l, rho_lR)

This script:
  1. solves the inner system at a chosen probe point (eps_v, p)
     using scipy.optimize.fsolve;
  2. evaluates the analytical IFT formulas at the converged state;
  3. re-solves at (eps_v ± h_eps, p) and (eps_v, p ± h_p), forms
     centered finite differences;
  4. compares (1)+(2) vs (3) — relative error should be < 1e-6.

What this does NOT verify
-------------------------
This is a check of the *derivation*, not of the mfront/native C++
implementation. A bug in the mfront translation (sign flip, missing
chain term, wrong variable capture) would NOT be caught here. To
verify the implementation, the same probe must be done by running
OGS itself — that requires a build, which is outside the scope of
this Python prototype.

Parameters / citations
----------------------
All material constants below come from
  Tests/Data/RichardsMechanics/ANCHORS_MS33_ModelI/ms33_model_i_dd1600.prj
which carries inline citations:
  - HamakerConstant 2.2e-20 J: Israelachvili & Adams 1978
    (J. Chem. Soc. Faraday Trans. I, vol. 74, p. 975, Table 2)
  - SpecificSurface 5.23e5 m^2/kg: FEBEX montmorillonite
  - VdwAugmentationPrefactor 29985.0 J/kg, VdwAugmentationDecayLength
    7.5e-7 m: paper_DSM.tex Tab. msIII_vs_villar (rho_d=1600)
  - phi0 0.424460: dd1600 reference porosity
  - MassExchangeCoefficient, FluidViscosity, densities: from the
    same PRJ (uncomment block citation in script header).

The probe point (eps_v, p, n_l_prev, rho_lR_prev) is a CHOICE made
to exercise an "interior" state where neither the EOS clamp nor the
saturation cap is active. Different choices test the same math at
different operating points; agreement should hold at any reasonable
state.

Usage
-----
    python3 scripts/fd_probe_dsm_ift.py

Exit code 0 if both relative errors < 1e-5; non-zero otherwise.
"""

from __future__ import annotations

import sys

import numpy as np
from scipy.optimize import fsolve

# ============================================================================
# Material parameters from ms33_model_i_dd1600.prj (block citation per CLAUDE.md
# §1.1: all values from one source, single citation at top of block).
# ============================================================================
HAMAKER = 2.2e-20            # J     (Israelachvili & Adams 1978; mica-water SFA)
SPECIFIC_SURFACE = 5.23e5    # m^2/kg (FEBEX montmorillonite, 523 m^2/g)
RHO_LR_REF = 1000.0          # kg/m^3 (ReferenceLiquidDensityMacro)
RHO_L0 = 1.0                 # kg/m^3 (ReferenceLiquidDensityMicro)
RHO_SR = 2780.0              # kg/m^3 (ReferenceDensitySolid)
DENSITY_A = 50.0             # MicroLiquidDensityA (a in EOS)
DENSITY_B = 1.0              # MicroLiquidDensityB (b in EOS)
PHI0 = 0.424460              # InitialPorosity (dd1600 anchor)
K_AUG = 29985.0              # J/kg  (paper_DSM.tex Tab. msIII_vs_villar, rho_d=1600)
LAMBDA_AUG = 7.5e-7          # m     (same source)
ALPHA_BAR = 1e-10            # MassExchangeCoefficient
MU_VISC = 1e-3               # FluidViscosity (Pa s)

PI = np.pi
MIN_NL = 1e-14
MIN_DENS_GAP = 1e-12

ALPHA_EFF = ALPHA_BAR * RHO_LR_REF / MU_VISC    # kg s / m^5

# ============================================================================
# Probe point: a representative interior state.
# eps_v small (linear-elastic regime); p in mid-suction band; n_l hydrated
# enough that no clamps engage.
# ============================================================================
P_TEST = -5e7              # Pa  (50 MPa suction)
EPS_V_TEST = -0.005        # 1   (small compression)
DT = 1.0                   # s   (arbitrary; only enters via dt*alpha_eff)
N_L_PREV = 0.05            # 1   (previous-step micro water content)
RHO_LR_PREV = 1000.05      # kg/m^3 (close to macro density; EOS is mild here)

# ============================================================================
# Core functions (mirror the mfront mu_micro_from_state etc.).
# ============================================================================

def active_nS(n_l: float) -> float:
    return max(MIN_NL, 1.0 - n_l)


def phi_M_from_nl(n_l: float) -> float:
    return float(np.clip((PHI0 - n_l) / (1.0 - n_l), 0.0, PHI0))


def phi_m_from_nl(n_l: float) -> float:
    return n_l * (1.0 - phi_M_from_nl(n_l))


def mu_micro(n_l: float, rho_lR: float) -> float:
    n_l_s = max(n_l, MIN_NL)
    rho_lR_s = max(rho_lR, RHO_LR_REF + MIN_DENS_GAP)
    nS = active_nS(n_l_s)
    mu_vdW = -(HAMAKER * SPECIFIC_SURFACE**3 / (6.0 * PI)) * \
              (nS**3 * RHO_SR**3) / (n_l_s**3 * rho_lR_s)
    if K_AUG > 0.0 and LAMBDA_AUG > 0.0:
        xi = n_l_s / (LAMBDA_AUG * nS * RHO_SR * SPECIFIC_SURFACE)
        mu_aug = -K_AUG * np.exp(-xi)
    else:
        mu_aug = 0.0
    return mu_vdW + mu_aug


def mu_macro(p: float) -> float:
    if p > -1e-12:
        return 0.0
    rho_macro_safe = max(RHO_LR_PREV, RHO_LR_REF + MIN_DENS_GAP)
    return p / rho_macro_safe


def rho_micro_eos(n_l: float, rho_lR: float) -> float:
    n_l_s = max(n_l, MIN_NL)
    nS = active_nS(n_l_s)
    omega = max(n_l_s * rho_lR / (nS * RHO_SR), MIN_NL)
    return RHO_L0 * np.exp(-DENSITY_A * omega**DENSITY_B) + RHO_LR_REF


def residual(state, eps_v: float, p: float):
    n_l, rho_lR = state
    phi_m = phi_m_from_nl(n_l)
    mu_LR = mu_macro(p)
    mu_lR_value = mu_micro(n_l, rho_lR)
    rho_l_hat = ALPHA_EFF * (mu_LR - mu_lR_value)
    prev_mass = phi_m_from_nl(N_L_PREV) * RHO_LR_PREV
    r1 = (phi_m * rho_lR - prev_mass - DT * rho_l_hat -
          phi_m * rho_lR * eps_v)
    r2 = rho_lR - rho_micro_eos(n_l, rho_lR)
    return [r1, r2]


def solve_microstate(eps_v: float, p: float,
                     initial_guess: tuple = (0.05, 1000.05)):
    sol, info, ier, msg = fsolve(
        residual, initial_guess, args=(eps_v, p),
        full_output=True, xtol=1e-13)
    return float(sol[0]), float(sol[1]), int(ier)


# ============================================================================
# Analytical Jacobian + IFT outputs.
# ============================================================================

def analytical_J(n_l: float, rho_lR: float, eps_v: float) -> dict:
    n_l_s = max(n_l, MIN_NL)
    rho_lR_s = max(rho_lR, RHO_LR_REF + MIN_DENS_GAP)
    nS = active_nS(n_l_s)
    one_minus_n_l = max(1e-12, 1.0 - n_l_s)
    phi_m = phi_m_from_nl(n_l_s)
    dphi_m_dnl = (1.0 - PHI0) / (one_minus_n_l ** 2)

    mu_vdW = -(HAMAKER * SPECIFIC_SURFACE**3 / (6.0 * PI)) * \
              (nS**3 * RHO_SR**3) / (n_l_s**3 * rho_lR_s)
    aug_active = (K_AUG > 0.0 and LAMBDA_AUG > 0.0)
    if aug_active:
        xi = n_l_s / (LAMBDA_AUG * nS * RHO_SR * SPECIFIC_SURFACE)
        mu_aug = -K_AUG * np.exp(-xi)
        dmu_aug_dnl = -mu_aug / (LAMBDA_AUG * nS**2 * RHO_SR * SPECIFIC_SURFACE)
    else:
        mu_aug = 0.0
        dmu_aug_dnl = 0.0

    dmu_lR_dnl = -3.0 * mu_vdW / (n_l_s * nS) + dmu_aug_dnl
    dmu_lR_drho_lR = -mu_vdW / rho_lR_s    # augmentation has no rho_lR dependence

    omega_raw = n_l_s * rho_lR_s / (nS * RHO_SR)
    omega_clamp = max(omega_raw, MIN_NL)
    EOS_attr = RHO_L0 * np.exp(-DENSITY_A * omega_clamp ** DENSITY_B)
    domega_dnl = (rho_lR_s / (RHO_SR * nS**2)) if omega_raw > MIN_NL else 0.0
    domega_drho_lR = (n_l_s / (nS * RHO_SR)) if omega_raw > MIN_NL else 0.0
    deos_pref = -EOS_attr * DENSITY_A * DENSITY_B * omega_clamp ** (DENSITY_B - 1.0)
    drho_eos_dnl = deos_pref * domega_dnl
    drho_eos_drho_lR = deos_pref * domega_drho_lR

    one_minus_dEv = 1.0 - eps_v
    J11 = dphi_m_dnl * rho_lR_s * one_minus_dEv + DT * ALPHA_EFF * dmu_lR_dnl
    J12 = phi_m * one_minus_dEv + DT * ALPHA_EFF * dmu_lR_drho_lR
    J21 = -drho_eos_dnl
    J22 = 1.0 - drho_eos_drho_lR
    det = J11 * J22 - J12 * J21

    return dict(J11=J11, J12=J12, J21=J21, J22=J22, det=det,
                phi_m=phi_m, mu_vdW=mu_vdW, mu_aug=mu_aug,
                dphi_m_dnl=dphi_m_dnl,
                dmu_lR_dnl=dmu_lR_dnl, dmu_lR_drho_lR=dmu_lR_drho_lR,
                drho_eos_dnl=drho_eos_dnl, drho_eos_drho_lR=drho_eos_drho_lR)


def analytical_dnl_deps_v(n_l: float, rho_lR: float, eps_v: float) -> float:
    J = analytical_J(n_l, rho_lR, eps_v)
    return J['phi_m'] * rho_lR * J['J22'] / J['det']


def analytical_dnl_dp(n_l: float, rho_lR: float, eps_v: float) -> float:
    J = analytical_J(n_l, rho_lR, eps_v)
    rho_macro_safe = max(RHO_LR_PREV, RHO_LR_REF + MIN_DENS_GAP)
    return DT * ALPHA_EFF * J['J22'] / (rho_macro_safe * J['det'])


# ============================================================================
# Probe driver.
# ============================================================================

def main() -> int:
    print(f"FD probe: DSM microstate analytical IFT vs. centered FD")
    print(f"  Probe point: eps_v={EPS_V_TEST}, p={P_TEST/1e6:.2f} MPa")
    print(f"  K_aug={K_AUG:.1f} J/kg, lambda={LAMBDA_AUG:.2e} m "
          f"(dd1600 calibrated)")
    print()

    # 1. Converged state.
    n_l_star, rho_lR_star, ier = solve_microstate(EPS_V_TEST, P_TEST)
    print(f"Converged microstate:")
    print(f"  n_l*    = {n_l_star:.10e}")
    print(f"  rho_lR* = {rho_lR_star:.10e} kg/m^3")
    print(f"  fsolve ier = {ier}")
    if ier != 1:
        print(f"  WARN: solver did not converge cleanly")

    J = analytical_J(n_l_star, rho_lR_star, EPS_V_TEST)
    print()
    print(f"Analytical Jacobian at converged state:")
    print(f"  J11 = {J['J11']:.6e}, J12 = {J['J12']:.6e}")
    print(f"  J21 = {J['J21']:.6e}, J22 = {J['J22']:.6e}")
    print(f"  det = {J['det']:.6e}")
    print(f"  phi_m = {J['phi_m']:.6e}, dphi_m/dnl = {J['dphi_m_dnl']:.6e}")
    print(f"  mu_vdW = {J['mu_vdW']:.6e} J/kg, "
          f"mu_aug = {J['mu_aug']:.6e} J/kg")
    print(f"  dmu_lR/dnl = {J['dmu_lR_dnl']:.6e}, "
          f"dmu_lR/drho_lR = {J['dmu_lR_drho_lR']:.6e}")

    # 2. Analytical IFT.
    dnl_deps_ana = analytical_dnl_deps_v(n_l_star, rho_lR_star, EPS_V_TEST)
    dnl_dp_ana = analytical_dnl_dp(n_l_star, rho_lR_star, EPS_V_TEST)
    print()
    print(f"Analytical IFT:")
    print(f"  dn_l/d(eps_v) = {dnl_deps_ana:.10e}     [1]")
    print(f"  dn_l/d(p)     = {dnl_dp_ana:.10e}     [1/Pa]")

    # 3. Numerical centered FD. Step sizes chosen O(sqrt(eps_machine)) at scale.
    h_eps = 1e-7                # eps_v scale ~ 1e-3 -> h/x ~ 1e-4
    h_p = 1e2                   # p scale ~ 5e7 -> h/p ~ 2e-6

    init = (n_l_star, rho_lR_star)
    n_l_eps_p, _, _ = solve_microstate(EPS_V_TEST + h_eps, P_TEST,
                                       initial_guess=init)
    n_l_eps_m, _, _ = solve_microstate(EPS_V_TEST - h_eps, P_TEST,
                                       initial_guess=init)
    dnl_deps_fd = (n_l_eps_p - n_l_eps_m) / (2.0 * h_eps)

    n_l_p_p, _, _ = solve_microstate(EPS_V_TEST, P_TEST + h_p,
                                     initial_guess=init)
    n_l_p_m, _, _ = solve_microstate(EPS_V_TEST, P_TEST - h_p,
                                     initial_guess=init)
    dnl_dp_fd = (n_l_p_p - n_l_p_m) / (2.0 * h_p)

    print()
    print(f"Centered FD (h_eps={h_eps:.0e}, h_p={h_p:.0e} Pa):")
    print(f"  dn_l/d(eps_v) = {dnl_deps_fd:.10e}     [1]")
    print(f"  dn_l/d(p)     = {dnl_dp_fd:.10e}     [1/Pa]")

    # 4. Relative error.
    def rel(a, f):
        denom = max(abs(a), abs(f), 1e-30)
        return abs(a - f) / denom

    rel_eps = rel(dnl_deps_ana, dnl_deps_fd)
    rel_p = rel(dnl_dp_ana, dnl_dp_fd)
    print()
    print(f"Relative agreement (analytical vs. FD):")
    print(f"  eps_v derivative: {rel_eps:.3e}  "
          f"{'PASS' if rel_eps < 1e-5 else 'FAIL'}")
    print(f"  p     derivative: {rel_p:.3e}  "
          f"{'PASS' if rel_p < 1e-5 else 'FAIL'}")

    # 5. Second probe at a different point (sanity).
    print()
    print("Second probe at (eps_v=-0.001, p=-1.0e7 Pa)...")
    eps2, p2 = -0.001, -1.0e7
    n2, r2, _ = solve_microstate(eps2, p2)
    d_eps2_a = analytical_dnl_deps_v(n2, r2, eps2)
    d_p2_a = analytical_dnl_dp(n2, r2, eps2)
    n_eps_p, _, _ = solve_microstate(eps2 + h_eps, p2, initial_guess=(n2, r2))
    n_eps_m, _, _ = solve_microstate(eps2 - h_eps, p2, initial_guess=(n2, r2))
    d_eps2_f = (n_eps_p - n_eps_m) / (2 * h_eps)
    n_p_p, _, _ = solve_microstate(eps2, p2 + h_p, initial_guess=(n2, r2))
    n_p_m, _, _ = solve_microstate(eps2, p2 - h_p, initial_guess=(n2, r2))
    d_p2_f = (n_p_p - n_p_m) / (2 * h_p)
    re_eps2 = rel(d_eps2_a, d_eps2_f)
    re_p2 = rel(d_p2_a, d_p2_f)
    print(f"  eps_v derivative: ana={d_eps2_a:.4e}, fd={d_eps2_f:.4e}, "
          f"rel={re_eps2:.3e}  {'PASS' if re_eps2 < 1e-5 else 'FAIL'}")
    print(f"  p     derivative: ana={d_p2_a:.4e}, fd={d_p2_f:.4e}, "
          f"rel={re_p2:.3e}  {'PASS' if re_p2 < 1e-5 else 'FAIL'}")

    all_pass = (rel_eps < 1e-5 and rel_p < 1e-5 and
                re_eps2 < 1e-5 and re_p2 < 1e-5)
    print()
    print("=" * 64)
    print(f"OVERALL: {'PASS' if all_pass else 'FAIL'}")
    print("=" * 64)
    return 0 if all_pass else 1


if __name__ == '__main__':
    sys.exit(main())
