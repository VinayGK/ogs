// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <limits>

#include "ProcessLib/RichardsMechanics/RichardsMechanicsFEM-impl.h"

using namespace ProcessLib::RichardsMechanics;

namespace
{
struct ReferenceVKSinglePointData
{
    double n_l = 0.0;
    VanDerWaalsMicroPotentialData micro_potential;
    PotentialDrivenMassExchangeData exchange;
    double p_L_m = 0.0;
    double S_L_m = 0.0;
    double phi_M = 0.0;
    double phi_m = 0.0;
};

double comparisonTolerance(double const a, double const b,
                           double const rel = 1e-10,
                           double const abs = 1e-14)
{
    return abs + rel * std::max(std::abs(a), std::abs(b));
}

ReferenceVKSinglePointData solveReferenceVKSinglePoint(
    double const p_L, double const n_l_prev, double const dt,
    double const rho_LR, double const alpha_bar, double const mu,
    double const phi, VKPotentialExchangeParameters const& vkp)
{
    constexpr double n_l_floor = 1e-16;

    auto const macro_potential =
        computeYoungLaplaceMacroPotential(p_L, rho_LR, vkp.pressure_tolerance);
    double const alpha_M_effective = alpha_bar * rho_LR / mu;

    auto const eval_exchange = [&](double const n_l)
    {
        auto const micro_potential = computeVanDerWaalsMicroPotential(
            n_l, rho_LR, vkp.micro_solid_volume_fraction_reference,
            vkp.micro_solid_density_reference, vkp.hamaker_constant,
            vkp.specific_surface);
        auto const exchange = computePotentialDrivenMassExchange(
            alpha_M_effective, macro_potential.mu_LR, micro_potential.mu_lR);
        return std::pair{micro_potential, exchange};
    };

    auto const residual = [&](double const n_l)
    {
        auto const [micro_potential, exchange] = eval_exchange(n_l);
        (void)micro_potential;
        return n_l - n_l_prev - dt * exchange.rho_l_hat / rho_LR;
    };

    double n_l = std::max(n_l_floor, n_l_prev);
    constexpr int max_iterations = 40;
    constexpr double residual_tolerance = 1e-14;
    constexpr double increment_tolerance = 1e-14;
    bool converged = false;

    for (int iter = 0; iter < max_iterations; ++iter)
    {
        double const r = residual(n_l);
        if (std::abs(r) <=
            residual_tolerance * std::max(1.0, std::abs(n_l_prev)))
        {
            converged = true;
            break;
        }

        double const h = 1e-8 * std::max(1.0, std::abs(n_l));
        double const n_l_plus = n_l + h;
        double const n_l_minus = std::max(n_l_floor, n_l - h);
        double const denom = n_l_plus - n_l_minus;
        EXPECT_GT(denom, 0.0);
        if (!(denom > 0.0))
        {
            return {};
        }

        double const jacobian =
            (residual(n_l_plus) - residual(n_l_minus)) / denom;
        EXPECT_TRUE(std::isfinite(jacobian));
        EXPECT_GT(std::abs(jacobian), 1e-20);
        if (!(std::isfinite(jacobian) && std::abs(jacobian) > 1e-20))
        {
            return {};
        }

        double step = -r / jacobian;
        double n_l_candidate = std::max(n_l_floor, n_l + step);

        // Basic backtracking to keep the independently coded reference solve
        // robust while remaining distinct from the production helper.
        double candidate_residual = residual(n_l_candidate);
        int backtracking_steps = 0;
        while (std::abs(candidate_residual) > std::abs(r) &&
               backtracking_steps < 12)
        {
            step *= 0.5;
            n_l_candidate = std::max(n_l_floor, n_l + step);
            candidate_residual = residual(n_l_candidate);
            ++backtracking_steps;
        }

        if (std::abs(n_l_candidate - n_l) <=
            increment_tolerance * std::max(1.0, std::abs(n_l)))
        {
            n_l = n_l_candidate;
            converged = true;
            break;
        }

        n_l = n_l_candidate;
    }

    EXPECT_TRUE(converged);
    if (!converged)
    {
        return {};
    }

    auto const [micro_potential, exchange] = eval_exchange(n_l);
    double const n_l_ref = std::max(
        1e-16, vkp.initial_micro_water_content.value_or(
                   vkp.micro_solid_volume_fraction_reference));
    double const phi_safe = std::max(0.0, phi);
    double const phi_m = std::clamp(n_l, 0.0, phi_safe);

    return {
        .n_l = n_l,
        .micro_potential = micro_potential,
        .exchange = exchange,
        .p_L_m = -rho_LR * micro_potential.mu_lR,
        .S_L_m = n_l / n_l_ref,
        .phi_M = phi_safe - phi_m,
        .phi_m = phi_m,
    };
}

double referenceDnLDpL(double const p_L, double const n_l_prev, double const dt,
                       double const rho_LR, double const alpha_bar,
                       double const mu, double const phi,
                       VKPotentialExchangeParameters const& vkp)
{
    double const h = 1e-8 * std::max(1.0, std::abs(p_L));
    auto const plus = solveReferenceVKSinglePoint(
        p_L + h, n_l_prev, dt, rho_LR, alpha_bar, mu, phi, vkp);
    auto const minus = solveReferenceVKSinglePoint(
        p_L - h, n_l_prev, dt, rho_LR, alpha_bar, mu, phi, vkp);
    return (plus.n_l - minus.n_l) / (2.0 * h);
}

double referenceDrhoLHatDpL(double const p_L, double const n_l_prev,
                            double const dt, double const rho_LR,
                            double const alpha_bar, double const mu,
                            double const phi,
                            VKPotentialExchangeParameters const& vkp)
{
    double const h = 1e-8 * std::max(1.0, std::abs(p_L));
    auto const plus = solveReferenceVKSinglePoint(
        p_L + h, n_l_prev, dt, rho_LR, alpha_bar, mu, phi, vkp);
    auto const minus = solveReferenceVKSinglePoint(
        p_L - h, n_l_prev, dt, rho_LR, alpha_bar, mu, phi, vkp);
    return ((-plus.exchange.rho_l_hat) - (-minus.exchange.rho_l_hat)) /
           (2.0 * h);
}
}  // namespace

TEST(RichardsMechanics, VKSingleIntegrationPointReferencePath)
{
    VKPotentialExchangeParameters vkp;
    vkp.enabled = true;
    vkp.pressure_tolerance = 0.0;
    vkp.hamaker_constant = 1e-30;
    vkp.specific_surface = 1.0;
    vkp.micro_solid_density_reference = 2650.0;
    vkp.micro_solid_volume_fraction_reference = 0.6;
    vkp.initial_micro_water_content = 0.1;
    vkp.local_jacobian_perturbation = 1e-8;

    double const p_L = -1.0e7;
    double const n_l_prev = 0.1;
    double const dt = 100.0;
    double const rho_LR = 1000.0;
    double const drho_LR_dpL = 0.0;
    double const alpha_bar = 1.0e-13;
    double const mu = 1.0e-3;
    double const phi = 0.4;
    double const phi_prev = 0.4;

    auto const macro_potential =
        computeYoungLaplaceMacroPotential(p_L, rho_LR, vkp.pressure_tolerance);
    auto const ogs_update = solveVKImplicitMicroWaterContent(
        n_l_prev, dt, rho_LR, alpha_bar, mu, macro_potential, vkp);
    ASSERT_TRUE(ogs_update.converged);

    auto const reference = solveReferenceVKSinglePoint(
        p_L, n_l_prev, dt, rho_LR, alpha_bar, mu, phi, vkp);

    EXPECT_NEAR(ogs_update.n_l, reference.n_l,
                comparisonTolerance(ogs_update.n_l, reference.n_l));
    EXPECT_NEAR(ogs_update.micro_potential.mu_lR, reference.micro_potential.mu_lR,
                comparisonTolerance(ogs_update.micro_potential.mu_lR,
                                    reference.micro_potential.mu_lR,
                                    1e-10, 1e-18));
    EXPECT_NEAR(ogs_update.exchange.rho_l_hat, reference.exchange.rho_l_hat,
                comparisonTolerance(ogs_update.exchange.rho_l_hat,
                                    reference.exchange.rho_l_hat,
                                    1e-10, 1e-18));
    EXPECT_NEAR(ogs_update.exchange.rho_L_hat, reference.exchange.rho_L_hat,
                comparisonTolerance(ogs_update.exchange.rho_L_hat,
                                    reference.exchange.rho_L_hat,
                                    1e-10, 1e-18));

    auto const compatibility_output =
        computeVKCompatibilityMicroHydraulicOutput(ogs_update.n_l, rho_LR, vkp);
    EXPECT_NEAR(compatibility_output.p_L_m, reference.p_L_m,
                comparisonTolerance(compatibility_output.p_L_m,
                                    reference.p_L_m, 1e-10, 1e-12));
    EXPECT_NEAR(compatibility_output.S_L_m, reference.S_L_m,
                comparisonTolerance(compatibility_output.S_L_m,
                                    reference.S_L_m));

    auto const transport_porosity_update =
        computeVKTransportPorosityUpdate(phi, phi_prev, ogs_update.n_l, n_l_prev);
    EXPECT_NEAR(transport_porosity_update.phi_M, reference.phi_M,
                comparisonTolerance(transport_porosity_update.phi_M,
                                    reference.phi_M));
    EXPECT_NEAR(transport_porosity_update.phi_m, reference.phi_m,
                comparisonTolerance(transport_porosity_update.phi_m,
                                    reference.phi_m));
    EXPECT_NEAR(transport_porosity_update.phi_M_prev, phi_prev - n_l_prev,
                comparisonTolerance(transport_porosity_update.phi_M_prev,
                                    phi_prev - n_l_prev));
    EXPECT_NEAR(transport_porosity_update.phi_m_prev, n_l_prev,
                comparisonTolerance(transport_porosity_update.phi_m_prev,
                                    n_l_prev));

    double const analytic_dn_l_dpL = computeVKImplicitNlDpL(
        dt, rho_LR, drho_LR_dpL, alpha_bar, mu, macro_potential,
        ogs_update.micro_potential, ogs_update.exchange);
    double const reference_dn_l_dpL = referenceDnLDpL(
        p_L, n_l_prev, dt, rho_LR, alpha_bar, mu, phi, vkp);

    EXPECT_NEAR(analytic_dn_l_dpL, reference_dn_l_dpL,
                comparisonTolerance(analytic_dn_l_dpL, reference_dn_l_dpL,
                                    5e-5, 1e-18));

    auto const fd_diagnostic = computeVKLocalJacobianDiagnosticData(
        n_l_prev, p_L, dt, rho_LR, drho_LR_dpL, alpha_bar, mu,
        vkp.pressure_tolerance, vkp);
    EXPECT_NEAR(fd_diagnostic.fd_dn_l_dpL, reference_dn_l_dpL,
                comparisonTolerance(fd_diagnostic.fd_dn_l_dpL,
                                    reference_dn_l_dpL, 5e-5, 1e-18));

    double const reference_drho_L_hat_dpL = referenceDrhoLHatDpL(
        p_L, n_l_prev, dt, rho_LR, alpha_bar, mu, phi, vkp);
    EXPECT_NEAR(fd_diagnostic.fd_drho_L_hat_dpL, reference_drho_L_hat_dpL,
                comparisonTolerance(fd_diagnostic.fd_drho_L_hat_dpL,
                                    reference_drho_L_hat_dpL,
                                    5e-5, 1e-18));
}
