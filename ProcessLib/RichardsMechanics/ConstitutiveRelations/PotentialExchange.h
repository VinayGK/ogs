// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <cmath>

#include "BaseLib/Error.h"

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
    double const p_LR, double const rho_LR, double const pressure_tolerance = 1.0)
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
    double dmu_lR_drho_lR = 0.0;
    double dmu_lR_dnS = 0.0;
    double dmu_lR_drho_SR = 0.0;
};

// DSM dsm_micromacro microscale vdW potential helper:
// omega_l = n_l * rho_lR / (nS * rho_SR)
// mu_lR = (A * rho_lR^3 / (6*pi)) * Sa^3
//          / (omega_l^3 + omega_min_vdw^3)
inline VanDerWaalsMicroPotentialData computeVanDerWaalsMicroPotential(
    double const n_l, double const rho_lR, double const nS, double const rho_SR,
    double const hamaker_constant, double const specific_surface,
    double const potential_sign_factor = 1.0)
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

    constexpr double pi = 3.141592653589793238462643383279502884;

    VanDerWaalsMicroPotentialData out;

    out.omega_l = n_l * rho_lR / (nS * rho_SR);

    out.domega_l_dnl = rho_lR / (nS * rho_SR);
    out.domega_l_drho_lR = n_l / (nS * rho_SR);
    out.domega_l_dnS = -out.omega_l / nS;
    out.domega_l_drho_SR = -out.omega_l / rho_SR;

    constexpr double h_min_vdw = 5e-11;
    double const omega_min_vdw = h_min_vdw * specific_surface;
    double const omega_min_vdw3 =
        omega_min_vdw * omega_min_vdw * omega_min_vdw;
    double const omega2 = out.omega_l * out.omega_l;
    double const omega3 = omega2 * out.omega_l;
    double const denominator = omega3 + omega_min_vdw3;
    double const prefactor = hamaker_constant * specific_surface *
                                 specific_surface * specific_surface *
                                 rho_lR * rho_lR * rho_lR /
                             (6.0 * pi);

    out.mu_lR = potential_sign_factor * prefactor / denominator;

    double const ddenominator_dnl = 3.0 * omega2 * out.domega_l_dnl;
    double const ddenominator_drho_lR =
        3.0 * omega2 * out.domega_l_drho_lR;
    double const ddenominator_dnS = 3.0 * omega2 * out.domega_l_dnS;
    double const ddenominator_drho_SR =
        3.0 * omega2 * out.domega_l_drho_SR;

    out.dmu_lR_dnl = -out.mu_lR * ddenominator_dnl / denominator;
    out.dmu_lR_drho_lR =
        out.mu_lR * (3.0 / rho_lR - ddenominator_drho_lR / denominator);
    out.dmu_lR_dnS = -out.mu_lR * ddenominator_dnS / denominator;
    out.dmu_lR_drho_SR = -out.mu_lR * ddenominator_drho_SR / denominator;

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
    if (!(alpha_M >= 0.0))
    {
        OGS_FATAL(
            "computePotentialDrivenMassExchange requires alpha_M >= 0, got "
            "{:g}.",
            alpha_M);
    }
    if (!std::isfinite(mu_LR) || !std::isfinite(mu_lR))
    {
        OGS_FATAL(
            "computePotentialDrivenMassExchange requires finite mu_LR and "
            "mu_lR, got mu_LR={:g}, mu_lR={:g}.",
            mu_LR, mu_lR);
    }

    PotentialDrivenMassExchangeData out;
    out.rho_l_hat = alpha_M * (mu_LR - mu_lR);
    out.rho_L_hat = -out.rho_l_hat;
    out.drho_l_hat_dmu_LR = alpha_M;
    out.drho_l_hat_dmu_lR = -alpha_M;
    out.drho_l_hat_dalpha_M = mu_LR - mu_lR;
    return out;
}
}  // namespace ProcessLib::RichardsMechanics
