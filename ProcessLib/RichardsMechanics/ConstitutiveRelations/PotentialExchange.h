// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <algorithm>
#include <cmath>
#include <limits>

#include "BaseLib/Error.h"
#include "ProcessLib/RichardsMechanics/PotentialExchangeParameters.h"

namespace ProcessLib::RichardsMechanics
{
struct YoungLaplaceMacroPotentialData
{
    double mu_LR = 0.0;
    double dmu_LR_dpLR = 0.0;
    double dmu_LR_drho_LR = 0.0;
    bool saturated_branch = true;
};

// DSM dsm_micromacro Phase-2 macro potential helper (Young-Laplace side):
// mu_LR = 0            for pLR > -ptol
// mu_LR = pLR / rho_LR otherwise
inline YoungLaplaceMacroPotentialData computeYoungLaplaceMacroPotential(
    double const p_LR, double const rho_LR, double const pressure_tolerance = 0.0)
{
    if (!(rho_LR > 0.0))
    {
        OGS_FATAL(
            "computeYoungLaplaceMacroPotential requires rho_LR > 0, got {:g}.",
            rho_LR);
    }
    if (!(pressure_tolerance >= 0.0))
    {
        OGS_FATAL(
            "computeYoungLaplaceMacroPotential requires pressure_tolerance >= "
            "0, got {:g}.",
            pressure_tolerance);
    }

    YoungLaplaceMacroPotentialData out;
    out.saturated_branch = (p_LR > -pressure_tolerance);
    if (out.saturated_branch)
    {
        return out;
    }

    out.mu_LR = p_LR / rho_LR;
    out.dmu_LR_dpLR = 1.0 / rho_LR;
    out.dmu_LR_drho_LR = -p_LR / (rho_LR * rho_LR);
    return out;
}

struct VanDerWaalsMicroPotentialData
{
    double omega_l = 0.0;
    double mu_lR = 0.0;

    double domega_l_dnl = 0.0;
    double domega_l_drho_lR = 0.0;
    double domega_l_dnS = 0.0;
    double domega_l_drho_SR = 0.0;

    double dmu_lR_dnl = 0.0;
    double dmu_lR_drho_lR = 0.0;   // NON-zero in the current form: computed as
                                   // -mu_lR/rho_lR (see line ~181, "/rho_lR fix"),
                                   // because mu_lR carries an explicit 1/rho_lR.
                                   // (Was exactly zero in the earlier reduced form.)
    double dmu_lR_dnS = 0.0;
    double dmu_lR_drho_SR = 0.0;
    // PARTIAL second derivative d^2 mu_lR/d n_l^2 (nS held fixed), needed by the
    // INTEGRABLE Maxwell mechanical partner (mu_lR_mech, below) for its analytic
    // n_l-derivative via Pi'' = -rho_lR*d^2 mu_lR_vdw/d n_l^2. Analytic, derived
    // in this file from the cubic core (12*mu/n_l^2) plus augmentation
    // (mu_aug*(xi/n_l)^2). The live-nS chain (F2) is NOT threaded into this
    // second partial (tangent-only path; the integrable partner's n_l-tangent
    // is used in the analytic predictor/scalar micro-solve Jacobians only).
    double d2mu_lR_dnl2 = 0.0;
    // PARTIAL derivatives w.r.t. the augmentation prefactor K (all other
    // state held fixed), needed by the live-K(rho_d) Jacobian chain
    // (K_OF_RHO_D_LIVE.md; Vinay 2026-06-12): every K-channel below is
    // LINEAR in K, so these are exact:
    //   mu_aug = sign*K*exp(-xi)        -> dmu_lR/dK        = sign*exp(-xi)
    //   d mu_aug/d n_l = -mu_aug*xi/n_l -> d(dmu_lR/dnl)/dK = -dmu_lR_dK*xi/n_l
    // Both are 0 when the augmentation is inactive (K == 0: the residual then
    // carries no aug term; the tangent at that isolated point is left 0) and
    // ddmu_lR_dnl_dK is 0 when the disjoining floor clamps (flat in n_l).
    double dmu_lR_dK = 0.0;        // [-] = (J/kg) per (J/kg)
    double ddmu_lR_dnl_dK = 0.0;   // [1/n_l] = (J/kg per n_l) per (J/kg)
};

// DSM dsm_micromacro microscale vdW potential helper:
// omega_l = n_l * rho_lR / (nS * rho_SR)
// mu_lR_vdW = (A * Sa^3 / (6*pi)) * (nS^3 * rho_SR^3) / (n_l^3 * rho_lR)   [J/kg]
//
//   A      = hamaker_constant   [J]    Hamaker constant for clay-water-clay
//                                      Literature: montmorillonite-water-montmorillonite
//                                      A = 2.2e-20 J (Israelachvili & Adams 1978, SFA mica proxy)
//                                      Range: 1–5e-20 J (smectite, DLVO literature)
//                                      DO NOT calibrate A — it is a material constant.
//
// Dimensional derivation:
//   Film thickness h = n_l / (nS * rho_SR * Sa)                         [m]
//   Surface energy/area: E = -A/(12*pi*h^2)                             [J/m^2]
//   Specific free energy: mu_lR_vdW = E * (nS*rho_SR*Sa) / (nS*rho_lR)  [J/kg]
//   = -A*Sa^3*nS^3*rho_SR^3 / (12*pi*n_l^3*rho_lR) (adsorption sign)
//   Disjoining pressure Pi=-dE/dh gives factor 2: A*Sa^3*nS^3*rho_SR^3 / (6*pi*n_l^3*rho_lR)
//   Consistent with p_L_m = -rho_lR * mu_lR  [Pa]  (impl.h lines 276, 1044)
//
// Optional lumped exponential augmentation (activated when
// potential_augmentation_prefactor > 0):
// h = n_l / (nS * rho_SR * Sa)               mean water film thickness [m]
// mu_lR_aug = K * exp(-h / lambda)
//   K      = potential_augmentation_prefactor    [J/kg]  augmentation amplitude (calibrate to Villar)
//   lambda = potential_augmentation_exponent [m]     characteristic film thickness (calibrate)
// Total: mu_lR = sign * (mu_lR_vdW + mu_lR_aug)   — ADDITIVE, augmentation never replaces vdW
// Setting K = 0 (default) reduces exactly to the pure vdW form.
//
// Optional disjoining-pressure FLOOR (n_l_floor > 0): the value formulas
// (mu_lR_vdW, mu_aug, omega_l, film thickness h) are evaluated at
//   n_l_eff = max(n_l, n_l_floor),
// so Pi ~ 1/n_l^3 is CAPPED at Pi(n_l_floor) instead of diverging as n_l -> 0.
// When clamped (n_l < n_l_floor) the value is CONSTANT in n_l, so the
// n_l-derivatives (dmu_lR_dnl, d2mu_lR_dnl2 and the augmentation chain) are 0;
// when n_l >= n_l_floor (incl. the default n_l_floor = 0) every formula and
// derivative is byte-identical to the unfloored form.
inline VanDerWaalsMicroPotentialData computeVanDerWaalsMicroPotential(
    double const n_l, double const rho_lR, double const nS, double const rho_SR,
    double const hamaker_constant, double const specific_surface,
    double const potential_sign_factor = 1.0,
    double const potential_augmentation_prefactor = 0.0,
    double const potential_augmentation_exponent = 0.0,
    double const dnS_dnl = 0.0, double const n_l_floor = 0.0)
{
    if (!(n_l > 0.0))
    {
        OGS_FATAL("computeVanDerWaalsMicroPotential requires n_l > 0, got {:g}.",
                  n_l);
    }
    if (!(rho_lR > 0.0))
    {
        OGS_FATAL(
            "computeVanDerWaalsMicroPotential requires rho_lR > 0, got {:g}.",
            rho_lR);
    }
    if (!(nS > 0.0))
    {
        OGS_FATAL("computeVanDerWaalsMicroPotential requires nS > 0, got {:g}.",
                  nS);
    }
    if (!(rho_SR > 0.0))
    {
        OGS_FATAL(
            "computeVanDerWaalsMicroPotential requires rho_SR > 0, got {:g}.",
            rho_SR);
    }
    if (!(hamaker_constant > 0.0))
    {
        OGS_FATAL(
            "computeVanDerWaalsMicroPotential requires hamaker_constant > 0, "
            "got {:g}.",
            hamaker_constant);
    }
    if (!(specific_surface > 0.0))
    {
        OGS_FATAL(
            "computeVanDerWaalsMicroPotential requires specific_surface > 0, "
            "got {:g}.",
            specific_surface);
    }
    if (!(potential_augmentation_prefactor >= 0.0))
    {
        OGS_FATAL(
            "computeVanDerWaalsMicroPotential requires "
            "potential_augmentation_prefactor >= 0, got {:g}.",
            potential_augmentation_prefactor);
    }
    if (potential_augmentation_prefactor > 0.0 &&
        !(potential_augmentation_exponent > 0.0))
    {
        OGS_FATAL(
            "computeVanDerWaalsMicroPotential requires "
            "potential_augmentation_exponent > 0 when "
            "potential_augmentation_prefactor > 0, got {:g}.",
            potential_augmentation_exponent);
    }

    if (!(n_l_floor >= 0.0))
    {
        OGS_FATAL(
            "computeVanDerWaalsMicroPotential requires n_l_floor >= 0, got "
            "{:g}.",
            n_l_floor);
    }

    constexpr double pi = 3.141592653589793238462643383279502884;

    // Disjoining-pressure floor: evaluate every VALUE formula at n_l_eff =
    // max(n_l, n_l_floor). `clamped` flags that n_l sits below the floor, where
    // Pi is held at Pi(n_l_floor) and is therefore FLAT in n_l (n_l-derivatives
    // -> 0). With the default n_l_floor = 0 the OGS_FATAL above already enforced
    // n_l > 0, so clamped == false and n_l_eff == n_l exactly (byte-identical).
    bool const clamped = (n_l_floor > 0.0) && (n_l < n_l_floor);
    double const n_l_eff = clamped ? n_l_floor : n_l;  // [-]

    VanDerWaalsMicroPotentialData out;

    out.omega_l = n_l_eff * rho_lR / (nS * rho_SR);

    // omega_l value uses n_l_eff; its n_l-derivative is 0 when clamped (the
    // clamped omega_l is constant in n_l). Unclamped: rho_lR/(nS*rho_SR) as before.
    out.domega_l_dnl = clamped ? 0.0 : rho_lR / (nS * rho_SR);
    out.domega_l_drho_lR = n_l_eff / (nS * rho_SR);
    out.domega_l_dnS = -out.omega_l / nS;
    out.domega_l_drho_SR = -out.omega_l / rho_SR;

    double const prefactor = hamaker_constant * specific_surface *
                                 specific_surface * specific_surface /
                             (6.0 * pi);
    // Units (n_l, nS are dimensionless volume fractions, so [1]):
    //   A · Sa^3 · nS^3 · rho_SR^3 / (n_l^3 · rho_lR)
    //   [J] · [m^2/kg]^3 · [1] · [kg/m^3]^3 / ([1] · [kg/m^3])
    //   = [J · m^6/kg^3 · kg^3/m^9] / [kg/m^3]
    //   = [J/m^3] / [kg/m^3]
    //   = J/kg  ✓
    // The leading /rho_lR converts J/m^3 (Pa) to J/kg — dimensionally
    // required by the exchange equation rho_l_hat = alpha * (mu_LR - mu_lR).
    // Value uses n_l_eff (the floored water content). When clamped this is the
    // capped Pi(n_l_floor); when unclamped (incl. floor = 0) it is the exact
    // unfloored mu_lR_vdW.
    out.mu_lR = potential_sign_factor * prefactor * (nS * nS * nS) *
                (rho_SR * rho_SR * rho_SR) /
                (n_l_eff * n_l_eff * n_l_eff * rho_lR);  // J/kg

    // PARTIAL derivatives (nS held fixed). The TOTAL d mu_lR/d n_l under
    // current_porosity_split (nS = 1 - n_l LIVE) additionally picks up the
    // dmu_lR_dnS * dnS/dnl chain, applied after the augmentation contributions
    // below via the dnS_dnl multiplier (F2, 2026-06-06; tangent-only).
    //
    // n_l-DERIVATIVES: 0 when clamped (capped mu_lR is flat in n_l), else the
    // exact cubic-core forms; divide by n_l_eff so the unclamped branch matches
    // the value's n_l_eff (== n_l unclamped). The rho_lR / nS / rho_SR partials
    // are NOT n_l-derivatives and stay exact (the capped value still depends on
    // them through n_l_eff).
    out.dmu_lR_dnl = clamped ? 0.0 : -3.0 * out.mu_lR / n_l_eff;  // J/kg
    out.dmu_lR_drho_lR = -out.mu_lR / rho_lR;  // non-zero after /rho_lR fix
    out.dmu_lR_dnS = 3.0 * out.mu_lR / nS;
    out.dmu_lR_drho_SR = 3.0 * out.mu_lR / rho_SR;
    // PARTIAL second derivative of the cubic core (mu_core ~ C/n_l^3, nS fixed):
    //   d^2 mu_core/d n_l^2 = 12*C/n_l^5 = 12*mu_core/n_l^2.   [J/kg]
    // 0 when clamped (flat in n_l).
    out.d2mu_lR_dnl2 =
        clamped ? 0.0 : 12.0 * out.mu_lR / (n_l_eff * n_l_eff);  // J/kg

    // Lumped exponential force augmentation:
    // h = n_l / (nS * rho_SR * Sa)  [mean water film thickness, m]
    // mu_lR_aug = sign * K * exp(-h / lambda)
    if (potential_augmentation_prefactor > 0.0)
    {
        // Film thickness h ~ n_l uses n_l_eff, so the augmentation value mu_aug
        // is capped at the floor exactly like the vdW core.
        double const xi = n_l_eff / (potential_augmentation_exponent * nS *
                                     rho_SR * specific_surface);
        double const mu_aug =
            potential_sign_factor * potential_augmentation_prefactor *
            std::exp(-xi);

        out.mu_lR += mu_aug;  // J/kg
        // n_l-derivative: 0 when clamped (mu_aug flat in n_l), else the exact
        // -mu_aug*xi/n_l_eff form.
        out.dmu_lR_dnl += clamped ? 0.0 : -mu_aug * xi / n_l_eff;  // J/kg
        // xi = n_l_eff / (lambda * nS * rho_SR * Sa) — independent of rho_lR,
        // so mu_aug has no rho_lR dependence; dmu_lR_drho_lR from vdW term only
        out.dmu_lR_dnS += mu_aug * xi / nS;
        out.dmu_lR_drho_SR += mu_aug * xi / rho_SR;
        // PARTIAL second derivative of the augmentation (nS fixed): xi linear in
        // n_l so xi/n_l = const, d mu_aug/d n_l = -mu_aug*(xi/n_l), hence
        //   d^2 mu_aug/d n_l^2 = mu_aug*(xi/n_l)^2 = mu_aug*xi^2/n_l^2.   [J/kg]
        // 0 when clamped (flat in n_l).
        out.d2mu_lR_dnl2 +=
            clamped ? 0.0 : mu_aug * xi * xi / (n_l_eff * n_l_eff);  // J/kg
        // Exact K-partials (mu_aug LINEAR in K; live-K(rho_d) Jacobian chain,
        // K_OF_RHO_D_LIVE.md). Clamped: value flat in n_l -> mixed partial 0.
        out.dmu_lR_dK = potential_sign_factor * std::exp(-xi);  // [-]
        out.ddmu_lR_dnl_dK =
            clamped ? 0.0 : -out.dmu_lR_dK * xi / n_l_eff;  // [1/n_l]
    }

    // ── F2 (2026-06-06, tangent-only): live-nS chain for dmu_lR/dnl ──────────
    // Under current_porosity_split nS = 1 - n_l is a function of n_l, so the
    // TOTAL derivative is dmu_lR/dnl + dmu_lR/dnS * dnS/dnl. The caller passes
    // dnS_dnl = -1 in that mode (and 0 in reference mode, where nS is constant
    // -> NO change, exact). out.dmu_lR_dnS here is the COMPLETE partial (vdW core
    // 3*mu/nS plus any augmentation term), so the chain is consistent across
    // both contributions. Tangent-only: the converged forward solve uses a
    // finite-difference Jacobian and is unaffected.
    //
    // Disjoining floor: when clamped, the disjoining law is held at Pi(n_l_floor)
    // and is FLAT in n_l (task spec: the n_l-derivatives of the clamped value are
    // 0), so the TOTAL n_l-tangent is 0 -> skip the live-nS chain as well.
    // Unclamped (incl. the default floor = 0) -> unchanged.
    if (!clamped)
    {
        out.dmu_lR_dnl += out.dmu_lR_dnS * dnS_dnl;
    }

    return out;
}

// RETIRED 2026-06-08 (Vinay's Option-B): superseded by
// computeIntegrableMechanicalMicroPotential; no callers. The whole
// computeMaxwellConjugateMicroPotential unit below (struct + function) is kept
// for the historical record only — the macro-exchange residual + Jacobian now
// route the mechanical Maxwell partner through the INTEGRABLE, ungated form
// (single mu_lR everywhere, equipresence). Do NOT add new callers.
// ── DSM Maxwell-conjugate term (B1) ─────────────────────────────────────────
// Restores the mean-effective-stress dependence of mu_lR — the Maxwell partner
// of the swelling eigenstress sigma_sw = -phi_m * Pi — so that (sigma, mu_lR)
// derive from ONE free energy Psi. Design + derivation:
//   ProcessLib/RichardsMechanics/DSM/MAXWELL_CONJUGATE_IMPLEMENTATION.md
//
//   mu_lR_mech = (1/rho_lR) * S1 * eps_v        for p_conf >= Pi   (sharp gate)
//              = 0                              otherwise
//   S1 = d sigma_sw,m / d n_l  (frozen phi, B1) [Pa] — supplied by the caller
//        from the eigenstress site (S1 = -n_S*(Pi + n_l*dPi_dnl)).
//
// Decisions (Vinay 2026-06-02): load EXPELS micro water; OGS effective stress,
// tension-positive; gate on confining pressure p_conf = -tr(sigma')/3 >= Pi;
// gate SHARP (the snap-drain at p_conf = Pi is the expected physical
// repercussion, not a numerical artifact); freeze phi when forming S1 (B1;
// B1.5 keep-phi-live deferred — see design doc ss.6.1).
struct MaxwellConjugateMicroPotentialData
{
    double mu_lR_mech = 0.0;           // J/kg
    double dmu_lR_mech_deps_v = 0.0;   // J/kg            per unit volumetric strain
    double dmu_lR_mech_dnl = 0.0;      // J/kg            per unit n_l
    double dmu_lR_mech_drho_lR = 0.0;  // (J/kg)/(kg/m^3)
    bool gate_open = false;
};

// RETIRED 2026-06-08 (Vinay's Option-B): superseded by
// computeIntegrableMechanicalMicroPotential; no callers. Definition retained
// (historical); do NOT re-wire.
// N3 (review 2026-06-14): re-confirmed fully dead — grep shows zero live
// callers in ProcessLib/ and Tests/ (only docs and this helper's own FATAL
// strings reference the name). Kept on disk per CLAUDE.md §6.3 (historical
// record, never delete); this banner is the deprecation marker.
inline MaxwellConjugateMicroPotentialData computeMaxwellConjugateMicroPotential(
    double const S1, double const dS1_dnl, double const eps_v,
    double const p_conf, double const Pi, double const rho_lR,
    double const n_S)
{
    if (!(rho_lR > 0.0))
    {
        OGS_FATAL(
            "computeMaxwellConjugateMicroPotential requires rho_lR > 0, got "
            "{:g}.",
            rho_lR);
    }
    // n_S = 1 - phi_M = REV macro-solid (aggregate) fraction. mu_lR is a
    // SPECIFIC potential [J/kg]; its thermodynamic conjugate is the micro liquid
    // mass per REV volume m_l = rho_lR * phi_m = rho_lR * (1 - phi_M) * n_l, with
    // d m_l / d n_l = rho_lR * n_S. The eigenstress half is already REV-referenced
    // (sigma_sw = -(1 - phi_M)*n_l*Pi), so the conjugate half must divide by the
    // REV micro-liquid mass density rho_lR * n_S, NOT rho_lR alone. Omitting n_S
    // under-references the term by (1 - phi_M): the two cross-partials
    // d sigma_sw,m/d m_l and d mu_lR_mech/d eps_v then differ by (1 - phi_M),
    // the integrability condition d sigma/d m_l = d mu_lR/d eps fails, and the
    // (sigma, mu_lR) response is non-conservative (loop integral != 0). See
    // DSM/MAXWELL_CONJUGATE_REV_REFERENCING.md.
    if (!(n_S > 0.0))
    {
        OGS_FATAL(
            "computeMaxwellConjugateMicroPotential requires n_S = 1 - phi_M > 0, "
            "got {:g}.",
            n_S);
    }
    MaxwellConjugateMicroPotentialData out;
    // Sharp gate (B1): the term switches on only once the confining pressure
    // reaches the disjoining pressure. Below it the films out-push the wall, so
    // stress cannot expel water and the term is EXACTLY zero (the model stays
    // identical to the pre-Maxwell code there). This branch is a constitutive
    // choice (sharp Pi), not a numerical guard.
    out.gate_open = (p_conf >= Pi);
    if (!out.gate_open)
    {
        return out;
    }
    // REV-referenced measure: divide by rho_lR * n_S (per-REV micro-liquid mass
    // density), not rho_lR. This scales the conjugate UP by 1/(1 - phi_M), the
    // reciprocal of the (1 - phi_M) the eigenstress carries -> the pair derives
    // from one Psi (energy-conserving). phi_M is frozen at the GP (B1), so
    // d n_S / d rho_lR = 0 and the rho_lR-derivative keeps its form.
    double const rho_rev = rho_lR * n_S;                 // kg/m^3, per-REV measure
    out.mu_lR_mech = S1 * eps_v / rho_rev;               // J/kg
    out.dmu_lR_mech_deps_v = S1 / rho_rev;               // J/kg per unit strain
    out.dmu_lR_mech_dnl = dS1_dnl * eps_v / rho_rev;     // J/kg per unit n_l
    out.dmu_lR_mech_drho_lR = -out.mu_lR_mech / rho_lR;  // (J/kg)/(kg/m^3)
    return out;
}

// ── Film-pressure micro potential (maxwell beamer sec.5, consolidated form) ──
// SUPERSEDES the strain-view computeMaxwellConjugateMicroPotential above when
// the film_pressure_coupling flag is ON. The two halves (disjoining + the
// mechanical Maxwell partner) fuse into ONE potential of the *film pressure*
// p_film = Pi - b*p_conf:
//
//   mu_lR = -(Pi - b*p_conf) / rho_lR = -p_film / rho_lR,
//
// with Pi = -rho_lR*mu_lR_vdw > 0 (disjoining) and p_conf = -tr(sigma_eff)/3 > 0
// in compression (confining; OGS sigma_eff tension-positive). The vdW helper
// already returns mu_lR_vdw = -Pi/rho_lR, so this routine returns ONLY the FILM
// DELTA to ADD on top of mu_lR_vdw:
//
//   delta = mu_lR_film - mu_lR_vdw = +g * b * p_conf / rho_lR,
//
// where g is a C1 activation (smooth gate) of x = (p_conf - Pi_gate),
// Pi_gate = phi_m*Pi = n_S*n_l*Pi (option-1 REV-consistent threshold, decision
// #4 2026-06-02), of width w = film_pressure_gate_width. As p_conf rises past
// Pi_gate, mu_lR rises -> rho_l_hat = alpha*(mu_LR - mu_lR) turns negative ->
// the micro DRAINS (the smooth replacement of the sharp Macaulay snap).
//
//   w == 0  -> g = Heaviside(x)  (sharp; reproduces the old gate SHAPE exactly).
//   w  > 0  -> g = smoothstep cubic t^2*(3-2t), t = x/w  (C1: g' = 0 at both
//              ends), so the drain emerges continuously.
//
// Derivatives (the film delta D = g(x)*b*p_conf/rho_lR; p_conf and rho_lR are
// independent of n_l in this partial, x depends on n_l only through Pi_gate):
//   dD/dp_conf = g'(x)*b*p_conf/rho_lR + g*b/rho_lR
//   dD/dn_l    = g'(x)*(-dPi_gate/dn_l)*b*p_conf/rho_lR,
//                dPi_gate/dn_l = n_S*(Pi + n_l*dPi_dnl)
//   dD/drho_lR = -D/rho_lR
//
// Equipresence note: the (1 - phi_M) contact-area factor that the eigenstress
// carries CANCELS the per-REV-mass referencing of mu_lR (sec.5 slide "Contact
// area and density"), so the SPECIFIC potential mu_lR = -p_film/rho_lR carries
// NO explicit (1-phi_M) and divides by the intrinsic rho_lR — exactly mirroring
// the disjoining term mu_lR_vdw = -Pi/rho_lR it augments. n_S enters ONLY
// through the gate threshold Pi_gate = phi_m*Pi. (Contrast the strain-view
// helper above, whose conjugate divides by rho_lR*n_S because its driver S1 is
// the REV-scaled eigenstress slope.)
struct FilmPressureMicroPotentialData
{
    double mu_lR_film_delta = 0.0;        // J/kg            (add to mu_lR_vdw)
    double dmu_lR_film_dnl = 0.0;         // J/kg            per unit n_l
    double dmu_lR_film_drho_lR = 0.0;     // (J/kg)/(kg/m^3)
    double dmu_lR_film_dp_conf = 0.0;     // (J/kg)/Pa
    double gate_value = 0.0;              // g in [0, 1]
    bool gate_open = false;              // g > 0
};

inline FilmPressureMicroPotentialData computeFilmPressureMicroPotential(
    double const Pi, double const p_conf, double const rho_lR,
    double const n_S, double const n_l, double const dmu_lR_vdw_dnl,
    double const gate_width_w, double const biot_b)
{
    (void)dmu_lR_vdw_dnl;  // reserved (the film delta's n_l-dependence is via the
                           // gate threshold only; vdW's own dmu_lR_dnl is folded
                           // separately at the call site).
    if (!(rho_lR > 0.0))
    {
        OGS_FATAL(
            "computeFilmPressureMicroPotential requires rho_lR > 0, got {:g}.",
            rho_lR);
    }
    if (!(n_S > 0.0))
    {
        OGS_FATAL(
            "computeFilmPressureMicroPotential requires n_S = 1 - phi_M > 0, "
            "got {:g}.",
            n_S);
    }
    FilmPressureMicroPotentialData out;

    // dPi/dn_l from the vdW chain: Pi = -rho_lR*mu_lR_vdw, density-agnostic form
    // (matches the assembly-site convention dPi_dnl = (Pi/mu_lR)*dmu_lR_dnl).
    double const mu_lR_vdw =
        (rho_lR > 0.0) ? -Pi / rho_lR : 0.0;  // mu_lR_vdw = -Pi/rho_lR
    double const dPi_dnl =
        (std::abs(mu_lR_vdw) > 1e-300) ? (Pi / mu_lR_vdw) * dmu_lR_vdw_dnl : 0.0;

    // REV-consistent gate threshold Pi_gate = phi_m*Pi = n_S*n_l*Pi.
    double const Pi_gate = n_S * n_l * Pi;
    double const dPi_gate_dnl = n_S * (Pi + n_l * dPi_dnl);
    double const x = p_conf - Pi_gate;  // gate argument [Pa]

    // C1 activation g(x) and dg/dx.
    double g = 0.0;
    double dg_dx = 0.0;
    if (!(gate_width_w > 0.0))
    {
        // Sharp fallback (w == 0): Heaviside step. dg/dx = 0 a.e. (the delta at
        // x = 0 is not represented — the constitutive choice of a sharp gate).
        g = (x >= 0.0) ? 1.0 : 0.0;
        dg_dx = 0.0;
    }
    else if (x <= 0.0)
    {
        g = 0.0;
        dg_dx = 0.0;
    }
    else if (x >= gate_width_w)
    {
        g = 1.0;
        dg_dx = 0.0;
    }
    else
    {
        double const t = x / gate_width_w;          // in (0, 1)
        g = t * t * (3.0 - 2.0 * t);                // smoothstep, C1
        dg_dx = 6.0 * t * (1.0 - t) / gate_width_w;  // g'(x)
    }

    out.gate_value = g;
    out.gate_open = (g > 0.0);

    // Film delta D = g * b * p_conf / rho_lR  (= +b*p_conf/rho_lR when gate open
    // and saturated; turns mu_lR -Pi/rho_lR into -(Pi - b*p_conf)/rho_lR).
    double const D = g * biot_b * p_conf / rho_lR;
    out.mu_lR_film_delta = D;
    // dx/dp_conf = +1, dx/dn_l = -dPi_gate/dn_l.
    out.dmu_lR_film_dp_conf =
        dg_dx * biot_b * p_conf / rho_lR + g * biot_b / rho_lR;
    out.dmu_lR_film_dnl =
        dg_dx * (-dPi_gate_dnl) * biot_b * p_conf / rho_lR;
    out.dmu_lR_film_drho_lR = -D / rho_lR;
    return out;
}

// ── INTEGRABLE Maxwell mechanical micro potential (the "web", spec item 2) ───
// REPLACES the non-integrable film bolt-on (computeFilmPressureMicroPotential's
// +g*b*p_conf/rho_lR delta) when film_pressure_coupling is ON. It is the RIGHT
// half of the Maxwell pair: the strain (eps_v) dependence of mu_lR that derives
// from the SAME free energy Psi as the swelling eigenstress
//   sigma_sw,m = -phi_m*p_film = -n_S*n_l*(Pi - b*p_conf)   (the WIP closure,
// kept UNCHANGED). The integrable partner is
//
//   mu_lR_mech = -[ (Pi + n_l*Pi')*eps_v + 0.5*b*K_drained*eps_v^2 ] / rho_lR
//                                                                    [J/kg]
//   mu_lR = mu_lR_vdw + mu_lR_mech   (ADDITIVE; never overwrites mu_lR_vdw).
//
// with Pi  = -rho_lR*mu_lR_vdw > 0 (disjoining), Pi' = dPi/dn_l, Pi'' = d2Pi/dn_l2
// (from the vdW + augmentation potential), b = biot_coefficient, K_drained the
// drained bulk modulus (from the elastic stiffness), eps_v the volumetric strain.
//
// Maxwell identity (the GP test asserts it FD-vs-analytic, with the drained
// skeleton p_conf = -K_drained*eps_v):
//   d sigma_sw,m/d n_l = n_S*rho_lR*d mu_lR/d eps_v.
//   LHS = -n_S*[ p_film + n_l*Pi' ] = -n_S*[ Pi + b*K_drained*eps_v + n_l*Pi' ]
//         (p_film = Pi - b*p_conf = Pi + b*K_drained*eps_v on the drained line),
//   d mu_lR/d eps_v = d mu_lR_mech/d eps_v
//                   = -[ (Pi + n_l*Pi') + b*K_drained*eps_v ] / rho_lR,
//   so n_S*rho_lR*(d mu_lR/d eps_v) = -n_S*[ (Pi + n_l*Pi') + b*K_drained*eps_v ]
//   = LHS.  EXACT (one Psi; loop integral zero), PARAMETER-FREE.
//
// Analytic derivatives (units in the comments below). NOTE: the linear-in-eps_v
// piece -(Pi + n_l*Pi')/rho_lR is the leading conjugate; the quadratic
// -0.5*b*K_drained*eps_v^2/rho_lR carries the p_conf(eps_v) = -K_drained*eps_v
// chain so that d mu_lR_mech/d eps_v reproduces the p_film slope exactly.
struct IntegrableMechanicalMicroPotentialData
{
    double mu_lR_mech = 0.0;           // J/kg
    double dmu_lR_mech_deps_v = 0.0;   // J/kg            per unit volumetric strain
    double dmu_lR_mech_dnl = 0.0;      // J/kg            per unit n_l
    double dmu_lR_mech_drho_lR = 0.0;  // (J/kg)/(kg/m^3)
};

inline IntegrableMechanicalMicroPotentialData
computeIntegrableMechanicalMicroPotential(
    double const Pi, double const dPi_dnl, double const d2Pi_dnl2,
    double const n_l, double const eps_v, double const biot_b,
    double const K_drained, double const rho_lR)
{
    if (!(rho_lR > 0.0))
    {
        OGS_FATAL(
            "computeIntegrableMechanicalMicroPotential requires rho_lR > 0, got "
            "{:g}.",
            rho_lR);
    }
    IntegrableMechanicalMicroPotentialData out;
    // A := Pi + n_l*Pi'  [Pa]; B := b*K_drained  [Pa].
    double const A = Pi + n_l * dPi_dnl;
    double const Bmod = biot_b * K_drained;
    // mu_lR_mech = -[ A*eps_v + 0.5*B*eps_v^2 ] / rho_lR   [J/kg]
    out.mu_lR_mech = -(A * eps_v + 0.5 * Bmod * eps_v * eps_v) / rho_lR;  // J/kg
    // d(mu_lR_mech)/d(eps_v) = -[ A + B*eps_v ] / rho_lR   [J/kg per unit strain]
    out.dmu_lR_mech_deps_v = -(A + Bmod * eps_v) / rho_lR;
    // d(mu_lR_mech)/d(n_l) = -[ (2*Pi' + n_l*Pi'')*eps_v ] / rho_lR
    //   (dA/dn_l = 2*Pi' + n_l*Pi''; B carries no n_l dependence). [J/kg per n_l]
    out.dmu_lR_mech_dnl =
        -((2.0 * dPi_dnl + n_l * d2Pi_dnl2) * eps_v) / rho_lR;
    // d(mu_lR_mech)/d(rho_lR) = -mu_lR_mech/rho_lR   [(J/kg)/(kg/m^3)]
    out.dmu_lR_mech_drho_lR = -out.mu_lR_mech / rho_lR;
    return out;
}

// ── Strained-film disjoining state — h(w_m, eps_v) ──────────────────────────
// Design: ProcessLib/RichardsMechanics/DSM/STRAINED_FILM_IMPLEMENTATION.md.
// Both variants reduce to evaluating the EXISTING bare law at an effective
// micro water content w_eff (the law depends on n_l only through the film
// thickness h = n_l/(nS*rho_SR*Sa), so straining h is straining the
// evaluation point):
//   Kinematic   (A): w_eff = n_l*(1 + kappa*eps_v), kappa = active_nS
//                    (Aggregate, the integrable completion of the existing
//                    eigenstress scale) or 1 (Unity). kappa is FROZEN at the
//                    GP (B1) — no d(kappa)/d(eps_v) chain.
//   Equilibrium (B): on the loaded branch (p_conf > Pi(n_l) > 0) w_eff solves
//                    Pi(w_eff) = p_conf (film force balance; emergent branch
//                    point, no bolted-on gate); else w_eff = n_l.
// Mass bookkeeping stays at n_l everywhere; only the disjoining evaluation
// point is strained.
struct StrainedFilmStateData
{
    double w_eff = 0.0;          // effective content fed to the law [-]
    double dw_eff_dnl = 1.0;     // d w_eff / d n_l [-]
    double dw_eff_deps_v = 0.0;  // d w_eff / d eps_v [-]
    bool loaded_branch = false;  // equilibrium mode: Pi(w_eff) = p_conf branch
};

// Invert the bare disjoining law Pi(w) = -rho_pi * mu_lR_bare(w) for w at a
// target pressure (equilibrium-spacing variant B). Pi is strictly decreasing
// in w (vdW core ~ w^-3 plus exponential augmentation), so the root in
// (w_floor, w_upper] is unique when it exists. Newton on f(w) = Pi(w) -
// p_target with the analytic dPi/dw from the law, seeded by the cubic-core
// inverse and guarded by bisection on the bracket. All law arguments mirror
// computeVanDerWaalsMicroPotential.
inline double invertDisjoiningPressure(
    double const p_target, double const w_upper, double const rho_lR,
    double const nS, double const rho_SR, double const hamaker_constant,
    double const specific_surface, double const potential_sign_factor,
    double const potential_augmentation_prefactor,
    double const potential_augmentation_exponent, double const n_l_floor,
    double const rho_pi)
{
    auto const Pi_and_slope = [&](double const w)
    {
        auto const v = computeVanDerWaalsMicroPotential(
            w, rho_lR, nS, rho_SR, hamaker_constant, specific_surface,
            potential_sign_factor, potential_augmentation_prefactor,
            potential_augmentation_exponent, 0.0 /*dnS_dnl*/, n_l_floor);
        return std::pair{-rho_pi * v.mu_lR, -rho_pi * v.dmu_lR_dnl};
    };

    // Bracket: Pi(w_upper) < p_target (caller-checked loaded branch). Lower
    // end: the law floor (Pi capped there) or a structural tiny fraction of
    // w_upper. 1e-8 is a bracket-width guard, not a physical value.
    double w_lo = std::max(n_l_floor, 1e-8 * w_upper);
    double w_hi = w_upper;
    auto const [Pi_lo, dPi_lo] = Pi_and_slope(w_lo);
    if (!(Pi_lo > p_target))
    {
        // Law cannot reach the target (e.g. floored Pi cap below p_target):
        // the film is squeezed to its cap; return the lower end.
        return w_lo;
    }

    // Seed: cubic-core inverse Pi ~ w^-3 => w0 = w_upper*(Pi(w_upper)/p)^(1/3).
    auto const [Pi_up, dPi_up] = Pi_and_slope(w_upper);
    double w = w_upper *
               std::cbrt(std::max(1e-300, Pi_up) / std::max(1e-300, p_target));
    w = std::clamp(w, w_lo, w_hi);

    // Iteration cap 50 (structural bound, not physics); relative tolerance on
    // the pressure residual scaled by the target per the tolerance-from-
    // problem-scale rule.
    for (int i = 0; i < 50; ++i)
    {
        auto const [Pi_w, dPi_dw] = Pi_and_slope(w);
        double const f = Pi_w - p_target;
        if (std::abs(f) <= 1e-12 * std::abs(p_target))
        {
            break;
        }
        // Maintain the bracket (Pi decreasing in w).
        if (f > 0.0)
        {
            w_lo = w;  // Pi too high -> root at larger w
        }
        else
        {
            w_hi = w;
        }
        double const step = (std::isfinite(dPi_dw) && std::abs(dPi_dw) > 0.0)
                                ? -f / dPi_dw
                                : 0.0;
        double w_next = w + step;
        if (!(w_next > w_lo && w_next < w_hi))
        {
            w_next = 0.5 * (w_lo + w_hi);  // bisection fallback
        }
        if (std::abs(w_next - w) <= 1e-15 * std::max(1.0, std::abs(w)))
        {
            w = w_next;
            break;
        }
        w = w_next;
    }
    return w;
}

inline StrainedFilmStateData computeStrainedFilmState(
    FilmStrainCouplingMode const mode, FilmStrainKappaMode const kappa_mode,
    double const n_l, double const active_nS, double const eps_v,
    double const p_conf, double const rho_lR, double const rho_SR,
    double const hamaker_constant, double const specific_surface,
    double const potential_sign_factor,
    double const potential_augmentation_prefactor,
    double const potential_augmentation_exponent, double const n_l_floor,
    double const rho_pi)
{
    StrainedFilmStateData out;
    out.w_eff = n_l;

    switch (mode)
    {
        case FilmStrainCouplingMode::Off:
            return out;
        case FilmStrainCouplingMode::Kinematic:
        {
            double const kappa =
                kappa_mode == FilmStrainKappaMode::Aggregate ? active_nS : 1.0;
            // Positivity guard on the spacing factor (1e-6 is a numeric
            // floor against w_eff <= 0 at extreme compression, not physics).
            double const f =
                std::max(1e-6, 1.0 + kappa * (std::isfinite(eps_v) ? eps_v
                                                                   : 0.0));
            out.w_eff = n_l * f;
            out.dw_eff_dnl = f;
            // kappa frozen at the GP (B1): d w_eff/d eps_v = n_l*kappa on the
            // unclamped branch, 0 when the positivity guard clamps.
            out.dw_eff_deps_v =
                (f > 1e-6) ? n_l * kappa : 0.0;
            return out;
        }
        case FilmStrainCouplingMode::Equilibrium:
        {
            if (!(std::isfinite(p_conf) && p_conf > 0.0))
            {
                return out;  // no load supplied -> unloaded branch
            }
            auto const bare = computeVanDerWaalsMicroPotential(
                n_l, rho_lR, active_nS, rho_SR, hamaker_constant,
                specific_surface, potential_sign_factor,
                potential_augmentation_prefactor,
                potential_augmentation_exponent, 0.0 /*dnS_dnl*/, n_l_floor);
            double const Pi_unloaded = -rho_pi * bare.mu_lR;
            if (!(Pi_unloaded > 0.0 && p_conf > Pi_unloaded))
            {
                return out;  // below the branch point: film carries the load
            }
            out.loaded_branch = true;
            out.w_eff = invertDisjoiningPressure(
                p_conf, n_l, rho_lR, active_nS, rho_SR, hamaker_constant,
                specific_surface, potential_sign_factor,
                potential_augmentation_prefactor,
                potential_augmentation_exponent, n_l_floor, rho_pi);
            // On the loaded branch the film state is pinned by the load:
            // w_eff solves Pi(w_eff) = p_conf, independent of n_l and eps_v
            // (implicit-function chains enter via p_conf only, which the
            // outer Newton iteration carries).
            out.dw_eff_dnl = 0.0;
            out.dw_eff_deps_v = 0.0;
            return out;
        }
    }
    return out;
}

// ── EXACT one-Psi strained-film energy pair (film_energy_route = exact) ─────
// Design: ProcessLib/RichardsMechanics/DSM/PI_OF_NL_EV_IMPLEMENTATION.md §2.1;
// closes STRAINED_FILM_IMPLEMENTATION.md §9a. Kinematic h-law only:
// w(e) = n_l*(1 + kappa*e). One energy
//
//   Psi_film(n_l, eps_v) = -(1-phi_M)*n_l*[ I_vdw + I_aug + S ],
//   I_T = int_0^{eps_v} Pi_T(w(e)) de  (closed form per term),
//   S   = 0.5*b*K_drained*eps_v^2      (transmitted-load work, route R3;
//                                       caller may drop it via include_S)
//
// gives BOTH halves by differentiation; Maxwell holds identically.
//
// WEIGHT vs FILM GEOMETRY (§6.7.5 split, Vinay's ruling 2026-09-02). The
// (1-phi_M) written above is the REV weight and is now ACTUALLY (1-phi_M):
// it is supplied by the caller as `eigenstress_rev_weight`. Before this split
// the routine had ONE `nS` argument doing two jobs, and the caller passed the
// film-geometry fraction (active_nS) into both, so the eigenstress and the
// energy silently carried the wrong weight. The two are numerically distinct:
//   1 - phi_M = (1-phi)/(1-n_l)   runs 0.5755 -> 1.0 over an MS33 run,
//   active_nS (Reference mode)    is frozen at n_s_ref = 0.5755395683453237.
// They coincide at t=0 only because the MS33 decks set n_s_ref = 1 - phi0.
// Because W now depends on n_l, dsigma_sw_dnl gains the (W'/W)*sigma_sw_m
// chain term; the caller supplies W'/W as `eigenstress_rev_weight_dln_dnl`.
// Both parameters default to the pre-split fused behaviour (W = nS_film,
// W' = 0), so every call site that does not opt in is bit-for-bit unchanged.
//
// In mu-space
// (mu_T = -Pi_T/rho_lR; sign conventions inherited from the bare law):
//
//   M_v = mu_v(n_l)*G3,  G3 = [1-(1+x)^-2]/(2*kappa),  x = kappa*eps_v
//   M_a = mu_a(n_l)*Gx,  Gx = -expm1(-xi0*x)/(xi0*kappa),
//                        xi0 = n_l_eff/(lambda*nS*rho_SR*Sa)
//   mu_mech       = -2*M_v + mu_a*[x*E/kappa - xi0*Gx] - 0.5*b*K_d*eps_v^2/rho_lR
//                   with E = exp(-xi0*x)                                  [J/kg]
//   dmu_mech/deps = -2*mu_v*(1+x)^-3 + mu_a*E*(1-xi0-xi0*x) - b*K_d*eps_v/rho_lR
//   dmu_mech/dnl  = 6*mu_v*G3/n_l
//                   - (xi0*mu_a/n_l)*[ (x*E/kappa)*(2+x) - xi0*Gx ]
//   dmu_mech/drho = -(vdW part + S part)/rho_lR   (mu_a carries no rho_lR)
//
// kappa->0 limits: G3, Gx -> eps_v and the pair reduces EXACTLY to
// computeIntegrableMechanicalMicroPotential (unit-tested) — unlike the
// OPERATIONAL route. Eigenstress half (for tests; the FEM eigenstress site is
// unchanged — it already evaluates Pi(w_eff) with the actual p_conf):
//   sigma_sw_m = -(1-phi_M)*n_l*[ Pi(w_eff) + b*K_d*eps_v ]   (drained line)
struct StrainedFilmEnergyPairData
{
    double Psi_film = 0.0;     // J/m^3 REV (drained-line S-form)
    double mu_mech = 0.0;      // J/kg
    double dmu_mech_dnl = 0.0;       // J/kg per unit n_l
    double dmu_mech_deps_v = 0.0;    // J/kg per unit strain
    double dmu_mech_drho_lR = 0.0;   // (J/kg)/(kg/m^3)
    double sigma_sw_m = 0.0;         // Pa (drained-line eigenstress half)
    double dsigma_sw_dnl = 0.0;      // Pa per unit n_l
    // Pre-cutoff bare-law reference values, exposed so the fold point can
    // recover the macro-floor cutoff factor g = mu_post/mu_pre and its chain.
    double mu_bare_pre = 0.0;        // J/kg
    double dmu_bare_dnl_pre = 0.0;   // J/kg per unit n_l
};

inline StrainedFilmEnergyPairData computeStrainedFilmEnergyPair(
    double const n_l, double const eps_v, double const kappa,
    double const biot_b, double const K_drained, bool const include_S,
    double const rho_lR, double const nS_film, double const rho_SR,
    double const hamaker_constant, double const specific_surface,
    double const potential_sign_factor,
    double const potential_augmentation_prefactor,
    double const potential_augmentation_exponent, double const dnS_dnl,
    double const n_l_floor,
    // ── §6.7.5 split (Vinay's ruling 2026-09-02) ────────────────────────────
    // W = the REV eigenstress/energy WEIGHT = 1 - phi_M [-]. NaN sentinel
    // (default) -> W = nS_film, i.e. the pre-split fused behaviour,
    // bit-for-bit. See the block comment above for why the two differ.
    double const eigenstress_rev_weight =
        std::numeric_limits<double>::quiet_NaN(),
    // d ln W / d n_l [1/n_l], the chain term the n_l-dependent weight
    // requires; 0 (default) = frozen weight, the legacy tangent.
    double const eigenstress_rev_weight_dln_dnl = 0.0)
{
    StrainedFilmEnergyPairData out;

    // ── Resolve the two DISTINCT solid fractions this routine consumes ─────
    // nS_film  : FILM GEOMETRY. Sets the film thickness h = n_l/(nS*rho_SR*Sa)
    //            and the decay ratio xi0, and enters the vdW core as nS^3. It
    //            is the AGGREGATE solid fraction that references the
    //            gravimetric water content (impl.h:444-468, the 2026-05-26
    //            active_nS fix; CLAUDE.md §2).
    // W        : REV EIGENSTRESS/ENERGY WEIGHT. The fraction of the REV the
    //            micro film occupies, phi_m = W*n_l. Documented as (1 - phi_M)
    //            in this header's own derivation above and in
    //            impl.h:2599 / PI_OF_NL_EV_IMPLEMENTATION.md:140, but silently
    //            inherited nS_film's value before this split.
    // They are NOT the same number: 1 - phi_M = (1-phi)/(1-n_l) runs 0.5755 ->
    // 1.0 over an MS33 run (MEASURED, exact-route VTUs 2026-09-02) while the
    // frozen reference nS is constant at 0.5755395683453237.
    bool const have_W = std::isfinite(eigenstress_rev_weight);
    double const W = have_W ? eigenstress_rev_weight : nS_film;  // [-]
    double const dlnW_dnl =
        have_W ? eigenstress_rev_weight_dln_dnl : 0.0;  // [1/n_l]

    // Per-term bare values at the TRUE n_l (floor handled inside the law).
    // vdW-only call (augmentation off) + full call; aug term by subtraction —
    // keeps the split in sync with any future change of the bare law.
    auto const full = computeVanDerWaalsMicroPotential(
        n_l, rho_lR, nS_film, rho_SR, hamaker_constant, specific_surface,
        potential_sign_factor, potential_augmentation_prefactor,
        potential_augmentation_exponent, dnS_dnl, n_l_floor);
    auto const vdw = computeVanDerWaalsMicroPotential(
        n_l, rho_lR, nS_film, rho_SR, hamaker_constant, specific_surface,
        potential_sign_factor, 0.0, 0.0, dnS_dnl, n_l_floor);
    double const mu_v = vdw.mu_lR;                // J/kg (cubic core)
    double const mu_a = full.mu_lR - vdw.mu_lR;   // J/kg (augmentation)
    out.mu_bare_pre = full.mu_lR;                 // J/kg
    out.dmu_bare_dnl_pre = full.dmu_lR_dnl;       // J/kg per n_l

    // Strained geometry x = kappa*eps_v with the SAME positivity guard and
    // freeze convention as computeStrainedFilmState (1e-6 numeric floor, not
    // physics; strain derivatives 0 when clamped).
    double const eps = std::isfinite(eps_v) ? eps_v : 0.0;
    double const x_raw = kappa * eps;
    bool const clamped_f = !(1.0 + x_raw > 1e-6);
    double const f = std::max(1e-6, 1.0 + x_raw);  // 1 + x, guarded
    double const x = f - 1.0;
    // x/kappa == eps algebraically; computing it as (f-1)/kappa amplifies the
    // rounding of 1+kappa*eps by 1/kappa (catastrophic for kappa -> 0 when
    // multiplied by a large mu_a). Use eps exactly on the unclamped branch.
    double const x_over_kappa = clamped_f ? x / kappa : eps;  // [-]

    // Floor-consistent xi0 (mirrors the bare law's n_l_eff clamping).
    bool const floored = (n_l_floor > 0.0) && (n_l < n_l_floor);
    double const n_l_eff = floored ? n_l_floor : n_l;
    double const xi0 =
        (potential_augmentation_prefactor > 0.0)
            ? n_l_eff / (potential_augmentation_exponent * nS_film * rho_SR *
                         specific_surface)
            : 0.0;  // unused when no augmentation

    // Stable strain integrals (series switch 1e-5: relative series error
    // < 1e-10 at the switch — scoped numeric default, PI_OF_NL_EV §4.6).
    double const y = xi0 * x;
    double G3;  // [-]; G3 -> eps_v as kappa -> 0
    if (std::abs(x) > 1e-5)
    {
        G3 = (1.0 - 1.0 / (f * f)) / (2.0 * kappa);
    }
    else
    {
        G3 = eps * (1.0 - 1.5 * x + 2.0 * x * x);  // series, O(x^3)
    }
    double Gx = 0.0;  // [-]; Gx -> eps_v as kappa -> 0
    double E = 1.0;   // exp(-xi0*x) [-]
    if (potential_augmentation_prefactor > 0.0)
    {
        E = std::exp(-y);
        if (std::abs(y) > 1e-5)
        {
            Gx = -std::expm1(-y) / (xi0 * kappa);
        }
        else
        {
            Gx = eps * (1.0 - 0.5 * y + y * y / 6.0);  // series, O(y^3)
        }
    }

    // Transmitted-load work S (route R3; NaN K_drained -> 0, mirrors the
    // shipped partner's sentinel convention).
    double const K_d =
        (include_S && std::isfinite(K_drained)) ? K_drained : 0.0;

    // ── mu half (J/kg; see derivation in the comment block above) ──────────
    double const mu_mech_vdw = -2.0 * mu_v * G3;                       // J/kg
    double const mu_mech_aug =
        (potential_augmentation_prefactor > 0.0)
            ? mu_a * (x_over_kappa * E - xi0 * Gx)
            : 0.0;                                                     // J/kg
    double const mu_mech_S = -0.5 * biot_b * K_d * eps * eps / rho_lR;  // J/kg
    out.mu_mech = mu_mech_vdw + mu_mech_aug + mu_mech_S;               // J/kg

    // d(mu_mech)/d(eps_v): strain derivatives 0 when the positivity guard
    // clamps (frozen geometry there, same convention as the kinematic state).
    if (!clamped_f)
    {
        out.dmu_mech_deps_v =
            -2.0 * mu_v / (f * f * f)                       // vdW    [J/kg per -]
            + ((potential_augmentation_prefactor > 0.0)
                   ? mu_a * E * (1.0 - xi0 - xi0 * x)
                   : 0.0)                                   // aug    [J/kg per -]
            - biot_b * K_d * eps / rho_lR;                  // S      [J/kg per -]
    }
    // d(mu_mech)/d(n_l): flat when the disjoining floor clamps (bare
    // derivatives are 0 there, same convention as the bare law).
    if (!floored)
    {
        out.dmu_mech_dnl =
            6.0 * mu_v * G3 / n_l_eff  // vdW [J/kg per n_l]
            - ((potential_augmentation_prefactor > 0.0)
                   ? (xi0 * mu_a / n_l_eff) *
                         ((x_over_kappa * E) * (2.0 + x) - xi0 * Gx)
                   : 0.0);  // aug [J/kg per n_l]
    }
    // d(mu_mech)/d(rho_lR): mu_v ~ 1/rho_lR and the S term ~ 1/rho_lR; the
    // augmentation mu_a carries no rho_lR.   [(J/kg)/(kg/m^3)]
    out.dmu_mech_drho_lR = -(mu_mech_vdw + mu_mech_S) / rho_lR;

    // ── eigenstress half (drained line; for the Maxwell/loop tests) ────────
    // Pi(w_eff) = -rho_lR*[mu_v*(1+x)^-3 + mu_a*E]   [Pa]
    double const Pi_weff =
        -rho_lR * (mu_v / (f * f * f) +
                   ((potential_augmentation_prefactor > 0.0) ? mu_a * E
                                                             : 0.0));  // Pa
    // Units: [-]*[-]*(Pa + [-]*Pa*[-]) = Pa  ✓  (§4.2)
    // Magnitude at the III exact-route endpoint (MEASURED 2026-09-02 from the
    // run VTU ts_982, nodal range over the mesh): W = 1 - phi_M = 1.000000
    // everywhere (phi_M has reached its 0 clamp) and n_l in [0.50766,
    // 0.510416], so phi_m = W*n_l in [0.508, 0.510] of the REV carries the
    // film — against phi_m in [0.292, 0.294] under the pre-split frozen
    // nS_film = 0.5755395683453237. The weight ratio at that state is
    // 1/0.5755395683 = 1.7375 = rho_s/rho_d = 2780/1600.  (§4.3)
    out.sigma_sw_m = -W * n_l * (Pi_weff + biot_b * K_d * eps);        // Pa
    // d(sigma_sw)/d(n_l) — exact chain through the n_l prefactor, through
    // w_eff = n_l*(1+x), AND (new, 2026-09-02) through the WEIGHT itself:
    //   d/dn_l[-W(n_l)*n_l*P] = -W*d/dn_l[n_l*P] + (W'/W)*sigma_sw_m
    // where d/dn_l[n_l*Pi_T(n_l*f)] is
    //   vdW: (1-3)*Pi_v(w_eff) = -2*Pi_v(w_eff)
    //   aug: Pi_a(w_eff)*(1 - xi0*f)
    // The (W'/W)*sigma_sw_m term is what makes the pair Maxwell-exact once W
    // depends on n_l; with the frozen weight (dlnW_dnl = 0, the default) it
    // vanishes and the legacy tangent is recovered bit-for-bit. It is the
    // SAME expression on both the floored and unfloored branches because
    // sigma_sw_m already carries the full -W*n_l*P product.
    if (!floored)
    {
        double const Pi_v_weff = -rho_lR * mu_v / (f * f * f);  // Pa
        double const Pi_a_weff =
            (potential_augmentation_prefactor > 0.0) ? -rho_lR * mu_a * E
                                                     : 0.0;  // Pa
        // Units: [-]*Pa + [1/n_l]*Pa = Pa per n_l  ✓  (§4.2)
        out.dsigma_sw_dnl =
            -W * (-2.0 * Pi_v_weff + Pi_a_weff * (1.0 - xi0 * f) +
                  biot_b * K_d * eps) +
            dlnW_dnl * out.sigma_sw_m;  // Pa per n_l
    }
    else
    {
        // Units: [-]*Pa + [1/n_l]*Pa = Pa per n_l  ✓  (§4.2)
        out.dsigma_sw_dnl = -W * (Pi_weff + biot_b * K_d * eps) +
                            dlnW_dnl * out.sigma_sw_m;  // Pa per n_l
    }

    // ── energy (drained-line S-form; for the loop test) ────────────────────
    // Psi_film = -(1-phi_M)*n_l*[I_vdw + I_aug + S], I_T = -rho_lR*M_T.
    // The SAME REV weight as sigma_sw_m — sigma_sw_m is dPsi_film/deps_v, so
    // if the two used different weights the pair would no longer come from
    // one Psi and the exact route would lose the property it exists for.
    // Units: [-]*[-]*((kg/m^3)*(J/kg) + [-]*Pa*[-]) = J/m^3 REV  ✓  (§4.2)
    // Magnitude at the III endpoint (MEASURED, ts_982): W = 1.000000 against
    // the pre-split 0.5755395683 — the stored film energy per REV rises by
    // the same 1.7375 factor as the eigenstress, as it must for one Psi.  (§4.3)
    out.Psi_film =
        -W * n_l *
        (-rho_lR * (mu_v * G3 + ((potential_augmentation_prefactor > 0.0)
                                     ? mu_a * Gx
                                     : 0.0)) +
         0.5 * biot_b * K_d * eps * eps);  // J/m^3 REV

    // ── AUDIT DIAGNOSTIC (2026-09-03, NOT a proposed fix) ──────────────────
    // Restore the W' channel the mu half is missing once W = W(n_l):
    //   Psi = -W*n_l*A  =>  dPsi/dn_l = W*rho*mu_mech - W'*n_l*A
    //                                 = W*rho*(mu_mech + (W'/W)*Psi/(W*rho))
    // so the n_l-conjugate of the SAME Psi is mu_mech + dlnW*Psi/(W*rho).
    // Its eps_v derivative is dlnW*sigma_sw_m/(W*rho). dmu_mech_dnl is left
    // alone (local-solve tangent only; costs iterations, not the converged
    // residual). Inert when dlnW_dnl == 0 (every default call).
    if (dlnW_dnl != 0.0)
    {
        double const inv = 1.0 / (W * rho_lR);            // m^3/kg per [-]
        out.mu_mech += dlnW_dnl * out.Psi_film * inv;     // J/kg
        out.dmu_mech_deps_v += dlnW_dnl * out.sigma_sw_m * inv;  // J/kg per -
    }

    return out;
}

struct PotentialDrivenMassExchangeData
{
    double rho_l_hat = 0.0;
    double rho_L_hat = 0.0;

    double drho_l_hat_dmu_LR = 0.0;
    double drho_l_hat_dmu_lR = 0.0;
    double drho_l_hat_dalpha_M = 0.0;
};

// DSM dsm_micromacro sign convention:
// rho_l_hat = alpha_M * (mu_LR - mu_lR)
// rho_L_hat = -rho_l_hat
inline PotentialDrivenMassExchangeData computePotentialDrivenMassExchange(
    double const alpha_M, double const mu_LR, double const mu_lR)
{
    PotentialDrivenMassExchangeData out;
    out.rho_l_hat = alpha_M * (mu_LR - mu_lR);
    out.rho_L_hat = -out.rho_l_hat;
    out.drho_l_hat_dmu_LR = alpha_M;
    out.drho_l_hat_dmu_lR = -alpha_M;
    out.drho_l_hat_dalpha_M = mu_LR - mu_lR;
    return out;
}
}  // namespace ProcessLib::RichardsMechanics
