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
    double dmu_lR_drho_lR = 0.0;   // exactly zero in the reduced algebraic form.
    double dmu_lR_dnS = 0.0;
    double dmu_lR_drho_SR = 0.0;
};

// DSM dsm_micromacro microscale vdW potential helper:
// omega_l = n_l * rho_lR / (nS * rho_SR)
// mu_lR = (A * rho_lR^3 / (6*pi)) * (Sa / omega_l)^3
//       = (A * Sa^3 / (6*pi)) * (nS^3 * rho_SR^3) / n_l^3
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

    double const prefactor = hamaker_constant * specific_surface *
                                 specific_surface * specific_surface /
                             (6.0 * pi);
    out.mu_lR = potential_sign_factor * prefactor * (nS * nS * nS) *
                (rho_SR * rho_SR * rho_SR) / (n_l * n_l * n_l);

    out.dmu_lR_dnl = -3.0 * out.mu_lR / n_l;
    out.dmu_lR_drho_lR = 0.0;
    out.dmu_lR_dnS = 3.0 * out.mu_lR / nS;
    out.dmu_lR_drho_SR = 3.0 * out.mu_lR / rho_SR;

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
