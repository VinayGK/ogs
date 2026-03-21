// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
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

double referenceMicroSolidVolumeFraction(
    double const n_l, double const phi, double const phi_M_prev,
    double const phi_m_prev, double const volumetric_strain,
    double const volumetric_strain_prev,
    VKPotentialExchangeParameters const& vkp)
{
    if (vkp.micro_solid_volume_fraction_mode ==
        VKMicroSolidVolumeFractionMode::Reference)
    {
        return std::max(1e-16, vkp.micro_solid_volume_fraction_reference);
    }

    auto const split = computeVKTransportPorosityUpdate(
        phi, phi_M_prev, phi_m_prev, n_l, volumetric_strain,
        volumetric_strain_prev, vkp.macro_porosity_update_mode);
    return std::max(1e-16, 1.0 - split.phi_M - split.phi_m);
}

VKReducedMicroLiquidDensityData solveReferenceReducedMicroLiquidDensity(
    double const n_l, double const rho_LR, double const nS,
    VKPotentialExchangeParameters const& vkp)
{
    auto const solve_rho = [&](double const n_eval)
    {
        double const n_l_safe = std::max(1e-16, n_eval);
        double const nS_safe = std::max(1e-16, nS);
        double const rho_SR =
            std::max(1e-16, vkp.micro_solid_density_reference);
        double const rho_l0 =
            std::max(1e-16, vkp.micro_liquid_density_reference);
        double const a_rho = std::max(1e-16, vkp.micro_liquid_density_a);
        double const b_rho = std::max(1e-16, vkp.micro_liquid_density_b);
        double const denominator = nS_safe * rho_SR;

        auto const rhs = [&](double const rho_lR)
        {
            double const omega_l =
                std::max(1e-16, n_l_safe * rho_lR / denominator);
            return std::pair{
                omega_l,
                rho_l0 * std::exp(-a_rho * std::pow(omega_l, b_rho)) +
                    rho_LR};
        };

        double rho_lR =
            rho_LR +
            rho_l0 *
                std::exp(-a_rho *
                         std::pow(std::max(1e-16, n_l_safe * rho_LR /
                                                      denominator),
                                  b_rho));
        constexpr int max_iterations = 40;
        for (int iter = 0; iter < max_iterations; ++iter)
        {
            auto const [omega_l, target] = rhs(rho_lR);
            double const residual = rho_lR - target;
            if (std::abs(residual) <=
                1e-14 * std::max(1.0, std::abs(rho_lR)))
            {
                return std::pair{rho_lR, omega_l};
            }

            double const h = 1e-8 * std::max(1.0, std::abs(rho_lR));
            double const rho_plus = rho_lR + h;
            double const rho_minus = std::max(1e-16, rho_lR - h);
            auto const [omega_plus, target_plus] = rhs(rho_plus);
            auto const [omega_minus, target_minus] = rhs(rho_minus);
            (void)omega_l;
            (void)omega_plus;
            (void)omega_minus;
            double const g_plus = rho_plus - target_plus;
            double const g_minus = rho_minus - target_minus;
            double const jacobian = (g_plus - g_minus) / (rho_plus - rho_minus);
            EXPECT_TRUE(std::isfinite(jacobian));
            EXPECT_GT(std::abs(jacobian), 1e-20);
            if (!(std::isfinite(jacobian) && std::abs(jacobian) > 1e-20))
            {
                break;
            }

            double const candidate =
                std::max(1e-16, rho_lR - residual / jacobian);
            if (std::abs(candidate - rho_lR) <=
                1e-14 * std::max(1.0, std::abs(rho_lR)))
            {
                auto const [omega_candidate, _] = rhs(candidate);
                (void)_;
                return std::pair{candidate, omega_candidate};
            }
            rho_lR = candidate;
        }
        auto const [omega_l, _] = rhs(rho_lR);
        (void)_;
        return std::pair{rho_lR, omega_l};
    };

    double const n_l_safe = std::max(1e-16, n_l);
    auto const [rho_lR, omega_l] = solve_rho(n_l_safe);
    double const h = 1e-8 * std::max(1.0, std::abs(n_l_safe));
    double const n_plus = n_l_safe + h;
    double const n_minus = std::max(1e-16, n_l_safe - h);
    double const rho_plus = solve_rho(n_plus).first;
    double const rho_minus = solve_rho(n_minus).first;
    double const drho_lR_dnl = (rho_plus - rho_minus) / (n_plus - n_minus);

    return {.rho_lR = rho_lR,
            .omega_l = omega_l,
            .drho_lR_dnl = drho_lR_dnl,
            .drho_l_dn_l = rho_lR + n_l_safe * drho_lR_dnl};
}

ReferenceVKSinglePointData solveReferenceVKSinglePoint(
    double const p_L, double const n_l_prev, double const dt,
    double const rho_LR, double const alpha_bar, double const mu,
    double const phi, VKPotentialExchangeParameters const& vkp,
    double const volumetric_strain = 0.0,
    double const volumetric_strain_prev = 0.0)
{
    constexpr double n_l_floor = 1e-16;
    double const phi_ceiling =
        vkp.local_nonlinear_solve_mode !=
                VKLocalNonlinearSolveMode::ScalarExchange &&
            std::isfinite(phi)
            ? std::max(n_l_floor, phi)
            : std::numeric_limits<double>::infinity();
    double const volumetric_strain_rate =
        dt > 0.0 ? (volumetric_strain - volumetric_strain_prev) / dt : 0.0;
    bool const use_mass_storage =
        vkp.local_nonlinear_solve_mode ==
        VKLocalNonlinearSolveMode::ScalarNotebookMassStorage;
    double const nS_prev = vkp.micro_solid_volume_fraction_mode ==
                                   VKMicroSolidVolumeFractionMode::Reference
                               ? vkp.micro_solid_volume_fraction_reference
                               : std::max(1e-16, 1.0 - 0.0 - n_l_prev);
    auto const prev_micro_liquid_density =
        use_mass_storage
            ? std::optional<VKReducedMicroLiquidDensityData>{
                  solveReferenceReducedMicroLiquidDensity(
                      n_l_prev, rho_LR, nS_prev, vkp)}
            : std::nullopt;
    double const rho_l_prev =
        prev_micro_liquid_density
            ? n_l_prev * prev_micro_liquid_density->rho_lR
            : 0.0;

    auto const macro_potential =
        computeYoungLaplaceMacroPotential(p_L, rho_LR, vkp.pressure_tolerance);
    double const alpha_M_effective = alpha_bar * rho_LR / mu;

    auto const eval_exchange = [&](double const n_l)
    {
        double const active_nS = referenceMicroSolidVolumeFraction(
            n_l, phi, 0.0, n_l_prev, volumetric_strain, volumetric_strain_prev,
            vkp);
        double const rho_lR_for_potential =
            use_mass_storage
                ? solveReferenceReducedMicroLiquidDensity(
                      n_l, rho_LR, active_nS, vkp)
                      .rho_lR
                : rho_LR;
        auto const micro_potential = computeVanDerWaalsMicroPotential(
            n_l, rho_lR_for_potential, active_nS,
            vkp.micro_solid_density_reference, vkp.hamaker_constant,
            vkp.specific_surface,
            microPotentialSignFactor(vkp.micro_potential_convention));
        auto const exchange = computePotentialDrivenMassExchange(
            alpha_M_effective, macro_potential.mu_LR, micro_potential.mu_lR);
        return std::pair{micro_potential, exchange};
    };

    auto const residual = [&](double const n_l)
    {
        auto const [micro_potential, exchange] = eval_exchange(n_l);
        (void)micro_potential;
        if (use_mass_storage)
        {
            double const active_nS = referenceMicroSolidVolumeFraction(
                n_l, phi, 0.0, n_l_prev, volumetric_strain,
                volumetric_strain_prev, vkp);
            auto const micro_liquid_density =
                solveReferenceReducedMicroLiquidDensity(
                    n_l, rho_LR, active_nS, vkp);
            double residual = n_l * micro_liquid_density.rho_lR - rho_l_prev -
                              dt * exchange.rho_l_hat;
            residual -= dt * n_l * micro_liquid_density.rho_lR *
                        volumetric_strain_rate;
            return residual;
        }

        double residual = n_l - n_l_prev - dt * exchange.rho_l_hat / rho_LR;
        if (vkp.local_nonlinear_solve_mode !=
            VKLocalNonlinearSolveMode::ScalarExchange)
        {
            residual -= dt * n_l * volumetric_strain_rate;
        }
        return residual;
    };

    double n_l = std::clamp(n_l_prev, n_l_floor, phi_ceiling);
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
        double n_l_candidate = std::clamp(n_l + step, n_l_floor, phi_ceiling);

        // Basic backtracking to keep the independently coded reference solve
        // robust while remaining distinct from the production helper.
        double candidate_residual = residual(n_l_candidate);
        int backtracking_steps = 0;
        while (std::abs(candidate_residual) > std::abs(r) &&
               backtracking_steps < 12)
        {
            step *= 0.5;
            n_l_candidate = std::clamp(n_l + step, n_l_floor, phi_ceiling);
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
                       VKPotentialExchangeParameters const& vkp,
                       double const volumetric_strain = 0.0,
                       double const volumetric_strain_prev = 0.0)
{
    double const h = 1e-8 * std::max(1.0, std::abs(p_L));
    auto const plus = solveReferenceVKSinglePoint(
        p_L + h, n_l_prev, dt, rho_LR, alpha_bar, mu, phi, vkp,
        volumetric_strain, volumetric_strain_prev);
    auto const minus = solveReferenceVKSinglePoint(
        p_L - h, n_l_prev, dt, rho_LR, alpha_bar, mu, phi, vkp,
        volumetric_strain, volumetric_strain_prev);
    return (plus.n_l - minus.n_l) / (2.0 * h);
}

double referenceDrhoLHatDpL(double const p_L, double const n_l_prev,
                            double const dt, double const rho_LR,
                            double const alpha_bar, double const mu,
                            double const phi,
                            VKPotentialExchangeParameters const& vkp,
                            double const volumetric_strain = 0.0,
                            double const volumetric_strain_prev = 0.0)
{
    double const h = 1e-8 * std::max(1.0, std::abs(p_L));
    auto const plus = solveReferenceVKSinglePoint(
        p_L + h, n_l_prev, dt, rho_LR, alpha_bar, mu, phi, vkp,
        volumetric_strain, volumetric_strain_prev);
    auto const minus = solveReferenceVKSinglePoint(
        p_L - h, n_l_prev, dt, rho_LR, alpha_bar, mu, phi, vkp,
        volumetric_strain, volumetric_strain_prev);
    return ((-plus.exchange.rho_l_hat) - (-minus.exchange.rho_l_hat)) /
           (2.0 * h);
}

enum class CoupledExchangeReferenceMode
{
    legacy_placeholder,
    full_potential_vdw
};

struct RepresentativeCoupledExchangeState
{
    char const* name = "";
    CoupledExchangeReferenceMode mode =
        CoupledExchangeReferenceMode::legacy_placeholder;
    double p_L = 0.0;
    double p_L_m = 0.0;
    double pressure_tolerance = 0.0;
    double n_l_prev = 0.0;
    double dt = 0.0;
    double rho_LR = 0.0;
    double drho_LR_dpL = 0.0;
    double alpha_bar = 0.0;
    double mu = 0.0;
    double phi = 0.0;
    double volumetric_strain = 0.0;
    double volumetric_strain_prev = 0.0;
};

double linearizedDensityAtPressure(double const p_L_eval, double const p_L_ref,
                                   double const rho_LR_ref,
                                   double const drho_LR_dpL)
{
    return std::max(1e-16, rho_LR_ref + drho_LR_dpL * (p_L_eval - p_L_ref));
}

double referenceCoupledRhoLHat(
    RepresentativeCoupledExchangeState const& state, double const p_L_eval,
    VKPotentialExchangeParameters const& vkp)
{
    double const rho_LR_eval = linearizedDensityAtPressure(
        p_L_eval, state.p_L, state.rho_LR, state.drho_LR_dpL);

    if (state.mode == CoupledExchangeReferenceMode::legacy_placeholder)
    {
        double const alpha_M_effective = state.alpha_bar * rho_LR_eval / state.mu;
        double const mu_LR_active = p_L_eval / rho_LR_eval;
        double const mu_lR_active = state.p_L_m / rho_LR_eval;
        auto const exchange = computePotentialDrivenMassExchange(
            alpha_M_effective, mu_LR_active, mu_lR_active);
        return exchange.rho_L_hat;
    }

    auto vkp_eval = vkp;
    vkp_eval.pressure_tolerance = state.pressure_tolerance;
    auto const reference = solveReferenceVKSinglePoint(
        p_L_eval, state.n_l_prev, state.dt, rho_LR_eval, state.alpha_bar,
        state.mu, state.phi, vkp_eval, state.volumetric_strain,
        state.volumetric_strain_prev);
    return reference.exchange.rho_L_hat;
}

double referenceCoupledDrhoLHatDpL(
    RepresentativeCoupledExchangeState const& state,
    VKPotentialExchangeParameters const& vkp)
{
    double const h = 1e-8 * std::max(1.0, std::abs(state.p_L));
    double const plus = referenceCoupledRhoLHat(state, state.p_L + h, vkp);
    double const minus = referenceCoupledRhoLHat(state, state.p_L - h, vkp);
    return (plus - minus) / (2.0 * h);
}

struct ProductionCoupledExchangeData
{
    double rho_L_hat = 0.0;
    double drho_L_hat_dpL = 0.0;
    bool converged = true;
};

ProductionCoupledExchangeData productionCoupledExchangeData(
    RepresentativeCoupledExchangeState const& state,
    VKPotentialExchangeParameters const& vkp)
{
    double const beta_LR = state.drho_LR_dpL / state.rho_LR;

    if (state.mode == CoupledExchangeReferenceMode::legacy_placeholder)
    {
        auto const data = computeVKPhase2CPlaceholderExchange(
            state.alpha_bar, state.mu, state.p_L, state.p_L_m, state.rho_LR,
            beta_LR, state.pressure_tolerance, false, false, 0.0, 0.0,
            false, 0.0, VKPotentialExchangeRoleMapping::CurrentOgs, false,
            vkp.fd_jacobian_perturbation);
        return {
            .rho_L_hat = data.exchange.rho_L_hat,
            .drho_L_hat_dpL = data.drho_L_hat_dpL_direct,
            .converged = true,
        };
    }

    auto const macro_potential = computeYoungLaplaceMacroPotential(
        state.p_L, state.rho_LR, state.pressure_tolerance);
    auto const n_l_update = solveVKImplicitMicroWaterContent(
        state.n_l_prev, state.dt, state.rho_LR, state.alpha_bar, state.mu,
        macro_potential,
        {.phi = state.phi,
         .volumetric_strain = state.volumetric_strain,
         .volumetric_strain_prev = state.volumetric_strain_prev},
        vkp);
    double const dn_l_dpL = computeVKImplicitNlDpL(
        state.n_l_prev, state.p_L, state.dt, state.rho_LR,
        state.drho_LR_dpL, state.alpha_bar, state.mu,
        macro_potential, n_l_update.micro_potential, n_l_update.exchange,
        {.phi = state.phi,
         .volumetric_strain = state.volumetric_strain,
         .volumetric_strain_prev = state.volumetric_strain_prev},
        vkp);
    double const dmu_lR_vdw_dpL =
        n_l_update.micro_potential.dmu_lR_dnl * dn_l_dpL +
        n_l_update.micro_potential.dmu_lR_drho_lR * state.drho_LR_dpL;

    auto const data = computeVKPhase2CPlaceholderExchange(
        state.alpha_bar, state.mu, state.p_L, state.p_L_m, state.rho_LR,
        beta_LR, state.pressure_tolerance, true, true,
        n_l_update.micro_potential.mu_lR,
        n_l_update.micro_potential.dmu_lR_drho_lR, true, dmu_lR_vdw_dpL,
        vkp.potential_role_mapping,
        false, vkp.fd_jacobian_perturbation);

    return {
        .rho_L_hat = data.exchange.rho_L_hat,
        .drho_L_hat_dpL = data.drho_L_hat_dpL_direct,
        .converged = n_l_update.converged,
    };
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
        n_l_prev, dt, rho_LR, alpha_bar, mu, macro_potential,
        {.phi = phi, .volumetric_strain = 0.0, .volumetric_strain_prev = 0.0},
        vkp);
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
        computeVKTransportPorosityUpdate(
            phi, phi_prev - n_l_prev, n_l_prev, ogs_update.n_l,
            /*volumetric_strain=*/0.0, /*volumetric_strain_prev=*/0.0,
            VKMacroPorosityUpdateMode::AlgebraicSplit);
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
        n_l_prev, p_L, dt, rho_LR, drho_LR_dpL, alpha_bar, mu, macro_potential,
        ogs_update.micro_potential, ogs_update.exchange,
        {.phi = phi, .volumetric_strain = 0.0, .volumetric_strain_prev = 0.0},
        vkp);
    double const reference_dn_l_dpL = referenceDnLDpL(
        p_L, n_l_prev, dt, rho_LR, alpha_bar, mu, phi, vkp);

    EXPECT_NEAR(analytic_dn_l_dpL, reference_dn_l_dpL,
                comparisonTolerance(analytic_dn_l_dpL, reference_dn_l_dpL,
                                    5e-5, 1e-18));

    auto const fd_diagnostic = computeVKLocalJacobianDiagnosticData(
        n_l_prev, p_L, dt, rho_LR, drho_LR_dpL, alpha_bar, mu,
        vkp.pressure_tolerance,
        {.phi = phi, .volumetric_strain = 0.0, .volumetric_strain_prev = 0.0},
        vkp);
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

TEST(RichardsMechanics, VKBranchSensitivityNearMacroPotentialTransition)
{
    VKPotentialExchangeParameters vkp;
    vkp.enabled = true;
    vkp.pressure_tolerance = 100.0;
    vkp.hamaker_constant = 1e-30;
    vkp.specific_surface = 1.0;
    vkp.micro_solid_density_reference = 2650.0;
    vkp.micro_solid_volume_fraction_reference = 0.6;
    vkp.initial_micro_water_content = 0.1;
    vkp.local_jacobian_perturbation = 1e-8;

    double const n_l_prev = 0.1;
    double const dt = 100.0;
    double const rho_LR = 1000.0;
    double const alpha_bar = 1.0e-13;
    double const mu = 1.0e-3;
    double const phi = 0.4;

    std::array<double, 5> const pressures = {
        -150.0,
        -100.0,
        -99.999,
        -50.0,
        0.0,
    };
    std::array<bool, 5> const saturated_expectation = {
        false,
        false,
        true,
        true,
        true,
    };

    struct CaseResult
    {
        double p_L = 0.0;
        bool saturated_branch = false;
        double mu_LR = 0.0;
        double n_l = 0.0;
        double rho_l_hat = 0.0;
        double p_L_m = 0.0;
        double S_L_m = 0.0;
    };

    std::array<CaseResult, 5> results;

    for (std::size_t i = 0; i < pressures.size(); ++i)
    {
        double const p_L = pressures[i];
        auto const macro_potential = computeYoungLaplaceMacroPotential(
            p_L, rho_LR, vkp.pressure_tolerance);

        EXPECT_EQ(macro_potential.saturated_branch, saturated_expectation[i]);
        if (saturated_expectation[i])
        {
            EXPECT_DOUBLE_EQ(macro_potential.mu_LR, 0.0);
            EXPECT_DOUBLE_EQ(macro_potential.dmu_LR_dpLR, 0.0);
        }
        else
        {
            EXPECT_LT(macro_potential.mu_LR, 0.0);
            EXPECT_DOUBLE_EQ(macro_potential.mu_LR, p_L / rho_LR);
        }

        auto const ogs_update = solveVKImplicitMicroWaterContent(
            n_l_prev, dt, rho_LR, alpha_bar, mu, macro_potential,
            {.phi = phi, .volumetric_strain = 0.0, .volumetric_strain_prev = 0.0},
            vkp);
        ASSERT_TRUE(ogs_update.converged);

        auto const reference = solveReferenceVKSinglePoint(
            p_L, n_l_prev, dt, rho_LR, alpha_bar, mu, phi, vkp);
        auto const compatibility_output =
            computeVKCompatibilityMicroHydraulicOutput(ogs_update.n_l, rho_LR,
                                                       vkp);

        EXPECT_NEAR(ogs_update.n_l, reference.n_l,
                    comparisonTolerance(ogs_update.n_l, reference.n_l));
        EXPECT_NEAR(ogs_update.exchange.rho_l_hat, reference.exchange.rho_l_hat,
                    comparisonTolerance(ogs_update.exchange.rho_l_hat,
                                        reference.exchange.rho_l_hat,
                                        1e-10, 1e-18));
        EXPECT_NEAR(compatibility_output.p_L_m, reference.p_L_m,
                    comparisonTolerance(compatibility_output.p_L_m,
                                        reference.p_L_m, 1e-10, 1e-12));
        EXPECT_NEAR(compatibility_output.S_L_m, reference.S_L_m,
                    comparisonTolerance(compatibility_output.S_L_m,
                                        reference.S_L_m));

        // Current kept VK branch: mu_LR <= 0 on the macro side, mu_lR > 0 on
        // the vdW microscale side. The exchange law therefore stays
        // sign-locked to non-increasing micro water content.
        EXPECT_LE(ogs_update.exchange.rho_l_hat, 0.0);
        EXPECT_LE(ogs_update.n_l,
                  n_l_prev + comparisonTolerance(ogs_update.n_l, n_l_prev,
                                                 0.0, 1e-18));

        results[i] = {
            .p_L = p_L,
            .saturated_branch = macro_potential.saturated_branch,
            .mu_LR = macro_potential.mu_LR,
            .n_l = ogs_update.n_l,
            .rho_l_hat = ogs_update.exchange.rho_l_hat,
            .p_L_m = compatibility_output.p_L_m,
            .S_L_m = compatibility_output.S_L_m,
        };
    }

    // Once the macro state is on the saturated helper branch, the active macro
    // potential is identically zero, so the local VK update becomes invariant
    // with respect to further increases in p_L as long as rho_LR stays fixed.
    for (std::size_t i = 3; i < results.size(); ++i)
    {
        EXPECT_NEAR(results[2].n_l, results[i].n_l,
                    comparisonTolerance(results[2].n_l, results[i].n_l));
        EXPECT_NEAR(results[2].rho_l_hat, results[i].rho_l_hat,
                    comparisonTolerance(results[2].rho_l_hat,
                                        results[i].rho_l_hat, 1e-10, 1e-18));
        EXPECT_NEAR(results[2].p_L_m, results[i].p_L_m,
                    comparisonTolerance(results[2].p_L_m, results[i].p_L_m,
                                        1e-10, 1e-12));
        EXPECT_NEAR(results[2].S_L_m, results[i].S_L_m,
                    comparisonTolerance(results[2].S_L_m, results[i].S_L_m));
    }

    // The kept branch remains monotone with respect to p_L: less negative
    // pressures produce less drying, but never net wetting.
    for (std::size_t i = 1; i < results.size(); ++i)
    {
        double const n_l_tolerance =
            comparisonTolerance(results[i - 1].n_l, results[i].n_l);
        EXPECT_LE(results[i - 1].n_l, results[i].n_l + n_l_tolerance);

        double const rho_l_hat_tolerance = comparisonTolerance(
            results[i - 1].rho_l_hat, results[i].rho_l_hat, 1e-10, 1e-18);
        EXPECT_LE(results[i - 1].rho_l_hat,
                  results[i].rho_l_hat + rho_l_hat_tolerance);
    }
}

TEST(RichardsMechanics, VKNegativeAttractiveMicroPotentialAdmitsWetting)
{
    VKPotentialExchangeParameters vkp;
    vkp.enabled = true;
    vkp.pressure_tolerance = 1.0;
    vkp.hamaker_constant = 6.0e-20;
    vkp.specific_surface = 1000.0;
    vkp.micro_solid_density_reference = 2650.0;
    vkp.micro_solid_volume_fraction_reference = 0.6;
    vkp.micro_potential_convention =
        VKMicroPotentialConvention::NegativeAttractive;
    vkp.initial_micro_water_content = 0.03;
    vkp.local_jacobian_perturbation = 1e-8;

    double const p_L = 0.0;
    double const n_l_prev = 0.03;
    double const dt = 1.0;
    double const rho_LR = 1000.0;
    double const alpha_bar = 1.0e-9;
    double const mu = 1.0e-3;
    double const phi = 0.4;

    auto const macro_potential =
        computeYoungLaplaceMacroPotential(p_L, rho_LR, vkp.pressure_tolerance);
    ASSERT_TRUE(macro_potential.saturated_branch);
    ASSERT_DOUBLE_EQ(macro_potential.mu_LR, 0.0);

    auto const ogs_update = solveVKImplicitMicroWaterContent(
        n_l_prev, dt, rho_LR, alpha_bar, mu, macro_potential,
        {.phi = phi, .volumetric_strain = 0.0, .volumetric_strain_prev = 0.0},
        vkp);
    ASSERT_TRUE(ogs_update.converged);

    auto const reference = solveReferenceVKSinglePoint(
        p_L, n_l_prev, dt, rho_LR, alpha_bar, mu, phi, vkp);
    auto const compatibility_output =
        computeVKCompatibilityMicroHydraulicOutput(ogs_update.n_l, rho_LR, vkp);

    EXPECT_LT(ogs_update.micro_potential.mu_lR, 0.0);
    EXPECT_GT(ogs_update.exchange.rho_l_hat, 0.0);
    EXPECT_LT(ogs_update.exchange.rho_L_hat, 0.0);
    EXPECT_GT(ogs_update.n_l, n_l_prev);
    EXPECT_GT(compatibility_output.p_L_m, 0.0);

    EXPECT_NEAR(ogs_update.n_l, reference.n_l,
                comparisonTolerance(ogs_update.n_l, reference.n_l));
    EXPECT_NEAR(ogs_update.exchange.rho_l_hat, reference.exchange.rho_l_hat,
                comparisonTolerance(ogs_update.exchange.rho_l_hat,
                                    reference.exchange.rho_l_hat,
                                    1e-10, 1e-18));
    EXPECT_NEAR(compatibility_output.p_L_m, reference.p_L_m,
                comparisonTolerance(compatibility_output.p_L_m,
                                    reference.p_L_m, 1e-10, 1e-12));
    EXPECT_NEAR(compatibility_output.S_L_m, reference.S_L_m,
                comparisonTolerance(compatibility_output.S_L_m,
                                    reference.S_L_m));
}

TEST(RichardsMechanics, VKVdWRelaxationStressIncrement)
{
    VKPotentialExchangeParameters vkp;
    vkp.enabled = true;
    vkp.micro_potential_convention =
        VKMicroPotentialConvention::NegativeAttractive;
    vkp.vdw_relaxation_stress_gain = 10.0;

    auto const& identity2 = MathLib::KelvinVector::Invariants<
        MathLib::KelvinVector::kelvin_vector_dimensions(2)>::identity2;

    auto const compressive_increment =
        computeVdWRelaxationStressIncrement<2>(5.0, 3.0, vkp);
    EXPECT_NEAR((compressive_increment + 20.0 * identity2).norm(), 0.0, 1e-14);

    auto const no_relaxation_increment =
        computeVdWRelaxationStressIncrement<2>(3.0, 5.0, vkp);
    EXPECT_NEAR(no_relaxation_increment.norm(), 0.0, 1e-14);

    vkp.vdw_relaxation_stress_gain = 0.0;
    auto const zero_gain_increment =
        computeVdWRelaxationStressIncrement<2>(5.0, 3.0, vkp);
    EXPECT_NEAR(zero_gain_increment.norm(), 0.0, 1e-14);

    vkp.vdw_relaxation_stress_gain = 10.0;
    vkp.micro_potential_convention = VKMicroPotentialConvention::PositiveReduced;
    auto const unsupported_convention_increment =
        computeVdWRelaxationStressIncrement<2>(5.0, 3.0, vkp);
    EXPECT_NEAR(unsupported_convention_increment.norm(), 0.0, 1e-14);
}

TEST(RichardsMechanics, VKMicroWaterContentStressIncrement)
{
    VKPotentialExchangeParameters vkp;
    vkp.enabled = true;
    vkp.micro_water_content_stress_gain = 10.0;

    auto const& identity2 = MathLib::KelvinVector::Invariants<
        MathLib::KelvinVector::kelvin_vector_dimensions(2)>::identity2;

    auto const compressive_increment =
        computeMicroWaterContentStressIncrement<2>(0.2, 0.3, vkp);
    EXPECT_NEAR((compressive_increment + 1.0 * identity2).norm(), 0.0, 1e-14);

    auto const no_growth_increment =
        computeMicroWaterContentStressIncrement<2>(0.3, 0.2, vkp);
    EXPECT_NEAR(no_growth_increment.norm(), 0.0, 1e-14);

    vkp.micro_water_content_stress_gain = 0.0;
    auto const zero_gain_increment =
        computeMicroWaterContentStressIncrement<2>(0.2, 0.3, vkp);
    EXPECT_NEAR(zero_gain_increment.norm(), 0.0, 1e-14);
}

TEST(RichardsMechanics, VKNotebookMicroPorositySwellingStressIncrement)
{
    using KM = MathLib::KelvinVector::KelvinMatrixType<2>;

    VKPotentialExchangeParameters vkp;
    vkp.enabled = true;
    vkp.micro_water_content_swelling_slope = 0.1;

    auto const& identity2 = MathLib::KelvinVector::Invariants<
        MathLib::KelvinVector::kelvin_vector_dimensions(2)>::identity2;
    KM C_el = KM::Identity();

    auto const loading_increment =
        computeNotebookMicroPorositySwellingStressIncrement<2>(
            0.2, 0.3, C_el, vkp);
    auto const expected_loading = -(0.1 * (0.3 - 0.2) / 3.0) * identity2;
    EXPECT_NEAR((loading_increment - expected_loading).norm(), 0.0, 1e-14);

    auto const unloading_increment =
        computeNotebookMicroPorositySwellingStressIncrement<2>(
            0.3, 0.2, C_el, vkp);
    auto const expected_unloading = -(0.1 * (0.2 - 0.3) / 3.0) * identity2;
    EXPECT_NEAR((unloading_increment - expected_unloading).norm(), 0.0, 1e-14);

    vkp.micro_water_content_swelling_slope = 0.0;
    auto const disabled_increment =
        computeNotebookMicroPorositySwellingStressIncrement<2>(
            0.2, 0.3, C_el, vkp);
    EXPECT_NEAR(disabled_increment.norm(), 0.0, 1e-14);
}

TEST(RichardsMechanics, VKNotebookAlignedSwellingIgnoresExploratoryGains)
{
    using KM = MathLib::KelvinVector::KelvinMatrixType<2>;
    using KV = MathLib::KelvinVector::KelvinVectorType<2>;

    VKPotentialExchangeParameters vkp;
    vkp.enabled = true;
    vkp.potential_role_mapping = VKPotentialExchangeRoleMapping::NotebookRoles;
    vkp.local_nonlinear_solve_mode =
        VKLocalNonlinearSolveMode::ScalarNotebookMassStorage;
    vkp.micro_water_content_swelling_slope = 0.1;
    vkp.vdw_relaxation_stress_gain = 100.0;
    vkp.micro_water_content_stress_gain = 100.0;

    KM C_el = KM::Identity();
    double const phi_m_prev = 0.2;
    double const phi_m = 0.3;
    double const n_l_prev = 0.2;
    double const n_l = 0.3;
    double const p_L_m_prev = 5.0;
    double const p_L_m = 3.0;

    KV const expected =
        computeNotebookMicroPorositySwellingStressIncrement<2>(
            phi_m_prev, phi_m, C_el, vkp);
    KV const actual =
        computeVKSwellingStressIncrement<2>(
            phi_m_prev, phi_m, n_l_prev, n_l, p_L_m_prev, p_L_m, C_el, vkp);

    EXPECT_NEAR((actual - expected).norm(), 0.0, 1e-14);
}

TEST(RichardsMechanics, VKTransportPorositySplitRecomposesTotalPorosity)
{
    auto const split = computeVKTransportPorosityUpdate(
        0.4, 0.27, 0.08, 0.1, 0.0, 0.0,
        VKMacroPorosityUpdateMode::AlgebraicSplit);

    EXPECT_NEAR(split.phi_m, 0.1, 1e-14);
    EXPECT_NEAR(split.phi_M, 0.3, 1e-14);
    EXPECT_NEAR(split.phi_m_prev, 0.08, 1e-14);
    EXPECT_NEAR(split.phi_M_prev, 0.27, 1e-14);
    EXPECT_NEAR(split.phi_M + split.phi_m, 0.4, 1e-14);
    EXPECT_NEAR(split.phi_M_prev + split.phi_m_prev, 0.35, 1e-14);

    auto const clamped = computeVKTransportPorosityUpdate(
        0.25, 0.0, 0.2, 0.4, 0.0, 0.0,
        VKMacroPorosityUpdateMode::AlgebraicSplit);
    EXPECT_NEAR(clamped.phi_m, 0.25, 1e-14);
    EXPECT_NEAR(clamped.phi_M, 0.0, 1e-14);
    EXPECT_NEAR(clamped.phi_M + clamped.phi_m, 0.25, 1e-14);
}

TEST(RichardsMechanics, VKNotebookAdditiveMacroPorosityRateUpdate)
{
    double const phi_M_prev = 0.30;
    double const phi_m_prev = 0.10;
    double const phi_m = 0.11;
    double const volumetric_strain_prev = 0.0;
    double const volumetric_strain = 1.0e-3;

    auto const split = computeVKTransportPorosityUpdate(
        0.4, phi_M_prev, phi_m_prev, phi_m, volumetric_strain,
        volumetric_strain_prev,
        VKMacroPorosityUpdateMode::NotebookAdditiveRate);

    double const delta_eps_v = volumetric_strain - volumetric_strain_prev;
    double const expected_phi_M =
        (phi_M_prev + (1.0 - phi_m) * delta_eps_v - (phi_m - phi_m_prev)) /
        (1.0 + delta_eps_v);

    EXPECT_NEAR(split.phi_m, phi_m, 1e-14);
    EXPECT_NEAR(split.phi_m_prev, phi_m_prev, 1e-14);
    EXPECT_NEAR(split.phi_M_prev, phi_M_prev, 1e-14);
    EXPECT_NEAR(split.phi_M, expected_phi_M, 1e-14);
    EXPECT_NEAR(split.phi_M + split.phi_m, expected_phi_M + phi_m, 1e-14);
    EXPECT_GT(split.phi_M + split.phi_m, phi_M_prev + phi_m_prev);
}

TEST(RichardsMechanics, VKCurrentPorositySplitMicroSolidFractionMode)
{
    VKPotentialExchangeParameters vkp;
    vkp.hamaker_constant = 6.0e-20;
    vkp.specific_surface = 4000.0;
    vkp.micro_solid_density_reference = 2650.0;
    vkp.micro_solid_volume_fraction_reference = 0.6;
    vkp.micro_solid_volume_fraction_mode =
        VKMicroSolidVolumeFractionMode::CurrentPorositySplit;
    vkp.macro_porosity_update_mode = VKMacroPorosityUpdateMode::AlgebraicSplit;
    vkp.initial_micro_water_content = 0.1;

    double const n_l = 0.1;
    double const rho_LR = 1000.0;
    VKLocalSolveContext const local_context{
        .phi = 0.25,
        .phi_M_prev = 0.15,
        .phi_m_prev = 0.1,
        .volumetric_strain = 0.0,
        .volumetric_strain_prev = 0.0};

    auto const active_nS =
        computeVKActiveMicroSolidVolumeFraction(n_l, local_context, vkp);
    EXPECT_NEAR(active_nS, 0.75, 1e-12);

    auto const active_output = computeVKCompatibilityMicroHydraulicOutput(
        n_l, rho_LR, local_context, vkp);
    auto const reference_output =
        computeVKCompatibilityMicroHydraulicOutput(n_l, rho_LR, vkp);

    double const expected_ratio =
        std::pow(active_nS / vkp.micro_solid_volume_fraction_reference, 3.0);
    EXPECT_NEAR(active_output.micro_potential.mu_lR /
                    reference_output.micro_potential.mu_lR,
                expected_ratio, 1e-12);
}

TEST(RichardsMechanics, VKReducedMicroLiquidDensityEOSReferencePath)
{
    VKPotentialExchangeParameters vkp;
    vkp.micro_solid_density_reference = 2650.0;
    vkp.micro_solid_volume_fraction_reference = 0.8;
    vkp.micro_liquid_density_reference = 1300.0;
    vkp.micro_liquid_density_a = 1.3;
    vkp.micro_liquid_density_b = 1.0;
    vkp.local_nonlinear_solve_mode =
        VKLocalNonlinearSolveMode::ScalarNotebookMassStorage;

    double const n_l = 0.1;
    double const rho_LR = 1000.0;
    double const nS = 0.8;

    auto const production =
        computeVKReducedMicroLiquidDensity(n_l, rho_LR, nS, vkp);
    auto const reference =
        solveReferenceReducedMicroLiquidDensity(n_l, rho_LR, nS, vkp);

    EXPECT_NEAR(production.rho_lR, reference.rho_lR,
                comparisonTolerance(production.rho_lR, reference.rho_lR,
                                    1e-9, 1e-14));
    EXPECT_NEAR(production.omega_l, reference.omega_l,
                comparisonTolerance(production.omega_l, reference.omega_l,
                                    1e-9, 1e-14));
    EXPECT_NEAR(production.drho_l_dn_l, reference.drho_l_dn_l,
                comparisonTolerance(production.drho_l_dn_l,
                                    reference.drho_l_dn_l, 1e-7, 1e-12));
}

TEST(RichardsMechanics, VKScalarNotebookStorageLocalSolveReferencePath)
{
    VKPotentialExchangeParameters vkp;
    vkp.enabled = true;
    vkp.pressure_tolerance = 0.0;
    vkp.hamaker_constant = 6.0e-20;
    vkp.specific_surface = 1000.0;
    vkp.micro_solid_density_reference = 2650.0;
    vkp.micro_solid_volume_fraction_reference = 0.6;
    vkp.micro_potential_convention =
        VKMicroPotentialConvention::NegativeAttractive;
    vkp.local_nonlinear_solve_mode =
        VKLocalNonlinearSolveMode::ScalarNotebookStorage;
    vkp.initial_micro_water_content = 0.03;

    double const p_L = 0.0;
    double const n_l_prev = 0.03;
    double const dt = 100.0;
    double const rho_LR = 1000.0;
    double const alpha_bar = 1.0e-9;
    double const mu = 1.0e-3;
    double const phi = 0.031;
    double const volumetric_strain_prev = 0.0;
    double const volumetric_strain = 1.0e-3;

    auto const macro_potential =
        computeYoungLaplaceMacroPotential(p_L, rho_LR, vkp.pressure_tolerance);
    auto const ogs_update = solveVKImplicitMicroWaterContent(
        n_l_prev, dt, rho_LR, alpha_bar, mu, macro_potential,
        {.phi = phi,
         .volumetric_strain = volumetric_strain,
         .volumetric_strain_prev = volumetric_strain_prev},
        vkp);
    ASSERT_TRUE(ogs_update.converged);

    auto const reference = solveReferenceVKSinglePoint(
        p_L, n_l_prev, dt, rho_LR, alpha_bar, mu, phi, vkp,
        volumetric_strain, volumetric_strain_prev);
    EXPECT_NEAR(ogs_update.n_l, reference.n_l,
                comparisonTolerance(ogs_update.n_l, reference.n_l));
    EXPECT_LE(ogs_update.n_l, phi + comparisonTolerance(ogs_update.n_l, phi));

    auto vkp_scalar = vkp;
    vkp_scalar.local_nonlinear_solve_mode =
        VKLocalNonlinearSolveMode::ScalarExchange;
    auto const scalar_update = solveVKImplicitMicroWaterContent(
        n_l_prev, dt, rho_LR, alpha_bar, mu, macro_potential,
        {.phi = phi,
         .volumetric_strain = volumetric_strain,
         .volumetric_strain_prev = volumetric_strain_prev},
        vkp_scalar);
    ASSERT_TRUE(scalar_update.converged);
    EXPECT_GT(scalar_update.n_l, phi);
    EXPECT_LE(ogs_update.n_l,
              scalar_update.n_l +
                  comparisonTolerance(ogs_update.n_l, scalar_update.n_l));
}

TEST(RichardsMechanics, VKScalarNotebookMassStorageLocalSolveReferencePath)
{
    VKPotentialExchangeParameters vkp;
    vkp.enabled = true;
    vkp.pressure_tolerance = 0.0;
    vkp.hamaker_constant = 6.0e-20;
    vkp.specific_surface = 1000.0;
    vkp.micro_solid_density_reference = 2650.0;
    vkp.micro_solid_volume_fraction_reference = 0.6;
    vkp.micro_liquid_density_reference = 1300.0;
    vkp.micro_liquid_density_a = 1.3;
    vkp.micro_liquid_density_b = 1.0;
    vkp.micro_potential_convention =
        VKMicroPotentialConvention::NegativeAttractive;
    vkp.potential_role_mapping =
        VKPotentialExchangeRoleMapping::NotebookRoles;
    vkp.local_nonlinear_solve_mode =
        VKLocalNonlinearSolveMode::ScalarNotebookMassStorage;
    vkp.initial_micro_water_content = 0.03;

    double const p_L = 0.0;
    double const n_l_prev = 0.03;
    double const dt = 100.0;
    double const rho_LR = 1000.0;
    double const alpha_bar = 1.0e-9;
    double const mu = 1.0e-3;
    double const phi = 0.031;
    double const volumetric_strain_prev = 0.0;
    double const volumetric_strain = 1.0e-3;

    auto const macro_potential =
        computeYoungLaplaceMacroPotential(p_L, rho_LR, vkp.pressure_tolerance);
    auto const ogs_update = solveVKImplicitMicroWaterContent(
        n_l_prev, dt, rho_LR, alpha_bar, mu, macro_potential,
        {.phi = phi,
         .volumetric_strain = volumetric_strain,
         .volumetric_strain_prev = volumetric_strain_prev},
        vkp);
    ASSERT_TRUE(ogs_update.converged);

    auto const reference = solveReferenceVKSinglePoint(
        p_L, n_l_prev, dt, rho_LR, alpha_bar, mu, phi, vkp,
        volumetric_strain, volumetric_strain_prev);
    EXPECT_NEAR(ogs_update.n_l, reference.n_l,
                comparisonTolerance(ogs_update.n_l, reference.n_l,
                                    1e-8, 1e-14));

    double const analytic_dn_l_dpL = computeVKImplicitNlDpL(
        n_l_prev, p_L, dt, rho_LR, 0.0, alpha_bar, mu, macro_potential,
        ogs_update.micro_potential, ogs_update.exchange,
        {.phi = phi,
         .volumetric_strain = volumetric_strain,
         .volumetric_strain_prev = volumetric_strain_prev},
        vkp);
    double const reference_dn_l_dpL = referenceDnLDpL(
        p_L, n_l_prev, dt, rho_LR, alpha_bar, mu, phi, vkp,
        volumetric_strain, volumetric_strain_prev);
    EXPECT_NEAR(analytic_dn_l_dpL, reference_dn_l_dpL,
                comparisonTolerance(analytic_dn_l_dpL, reference_dn_l_dpL,
                                    1e-6, 1e-12));
}

TEST(RichardsMechanics, VKNotebookMassStorageCoupledSolveResidual)
{
    VKPotentialExchangeParameters vkp;
    vkp.enabled = true;
    vkp.pressure_tolerance = 0.0;
    vkp.hamaker_constant = 6.0e-20;
    vkp.specific_surface = 1000.0;
    vkp.micro_solid_density_reference = 2650.0;
    vkp.micro_solid_volume_fraction_reference = 0.6;
    vkp.micro_liquid_density_reference = 1300.0;
    vkp.micro_liquid_density_a = 1.3;
    vkp.micro_liquid_density_b = 1.0;
    vkp.micro_potential_convention =
        VKMicroPotentialConvention::NegativeAttractive;
    vkp.potential_role_mapping = VKPotentialExchangeRoleMapping::NotebookRoles;
    vkp.local_nonlinear_solve_mode =
        VKLocalNonlinearSolveMode::ScalarNotebookMassStorage;
    vkp.initial_micro_water_content = 0.03;

    double const p_L = 0.0;
    double const n_l_prev = 0.03;
    double const dt = 100.0;
    double const rho_LR = 1000.0;
    double const alpha_bar = 1.0e-9;
    double const mu = 1.0e-3;
    double const phi = 0.031;
    double const volumetric_strain_prev = 0.0;
    double const volumetric_strain = 1.0e-3;

    auto const macro_potential =
        computeYoungLaplaceMacroPotential(p_L, rho_LR, vkp.pressure_tolerance);
    VKLocalSolveContext const local_context{
        .phi = phi,
        .phi_M_prev = phi - n_l_prev,
        .phi_m_prev = n_l_prev,
        .volumetric_strain = volumetric_strain,
        .volumetric_strain_prev = volumetric_strain_prev};
    auto const prev_micro_liquid_density =
        computeVKPreviousMicroLiquidDensity(n_l_prev, rho_LR, local_context,
                                            vkp);
    double const rho_l_prev = n_l_prev * prev_micro_liquid_density.rho_lR;

    auto const coupled_update = solveVKNotebookMassStorageCoupledState(
        n_l_prev, rho_l_prev, prev_micro_liquid_density.rho_lR, dt, rho_LR,
        alpha_bar, mu, macro_potential, local_context, vkp);
    ASSERT_TRUE(coupled_update.converged);

    double const active_nS =
        computeVKActiveMicroSolidVolumeFraction(coupled_update.n_l,
                                                 local_context, vkp);
    auto const micro_potential = computeVanDerWaalsMicroPotential(
        coupled_update.n_l, coupled_update.rho_lR, active_nS,
        vkp.micro_solid_density_reference, vkp.hamaker_constant,
        vkp.specific_surface, vkMicroPotentialSignFactor(vkp));
    double const mu_LR_active = micro_potential.mu_lR;
    double const mu_lR_active = macro_potential.mu_LR;
    double const alpha_M_effective = alpha_bar * rho_LR / mu;
    auto const exchange = computePotentialDrivenMassExchange(
        alpha_M_effective, mu_LR_active, mu_lR_active);

    double const volumetric_strain_rate =
        (volumetric_strain - volumetric_strain_prev) / dt;
    double const rho_l = coupled_update.n_l * coupled_update.rho_lR;
    double const mass_residual = rho_l - rho_l_prev -
                                 dt * exchange.rho_l_hat -
                                 dt * rho_l * volumetric_strain_rate;
    auto const micro_density = computeVKReducedMicroLiquidDensity(
        coupled_update.n_l, rho_LR, active_nS, vkp);
    double const density_residual =
        coupled_update.rho_lR - micro_density.rho_lR;

    EXPECT_NEAR(mass_residual, 0.0, 1e-10);
    EXPECT_NEAR(density_residual, 0.0, 1e-10);
}

TEST(RichardsMechanics, VKNotebookMassStorageCoupledSolveResiduals)
{
    VKPotentialExchangeParameters vkp;
    vkp.enabled = true;
    vkp.pressure_tolerance = 0.0;
    vkp.hamaker_constant = 6.0e-20;
    vkp.specific_surface = 1000.0;
    vkp.micro_solid_density_reference = 2650.0;
    vkp.micro_solid_volume_fraction_reference = 0.6;
    vkp.micro_liquid_density_reference = 1300.0;
    vkp.micro_liquid_density_a = 1.3;
    vkp.micro_liquid_density_b = 1.0;
    vkp.micro_potential_convention =
        VKMicroPotentialConvention::NegativeAttractive;
    vkp.local_nonlinear_solve_mode =
        VKLocalNonlinearSolveMode::ScalarNotebookMassStorage;
    vkp.macro_porosity_update_mode =
        VKMacroPorosityUpdateMode::NotebookAdditiveRate;
    vkp.micro_solid_volume_fraction_mode =
        VKMicroSolidVolumeFractionMode::CurrentPorositySplit;
    vkp.initial_micro_water_content = 0.05;

    double const p_L = 0.0;
    double const n_l_prev = 0.05;
    double const dt = 100.0;
    double const rho_LR = 1000.0;
    double const alpha_bar = 1.0e-9;
    double const mu = 1.0e-3;
    double const phi = 0.26;
    double const phi_M_prev = 0.18;
    double const phi_m_prev = 0.08;
    double const volumetric_strain_prev = 0.0;
    double const volumetric_strain = 1.0e-3;

    VKLocalSolveContext const local_context{
        .phi = phi,
        .phi_M_prev = phi_M_prev,
        .phi_m_prev = phi_m_prev,
        .volumetric_strain = volumetric_strain,
        .volumetric_strain_prev = volumetric_strain_prev,
    };

    auto const macro_potential =
        computeYoungLaplaceMacroPotential(p_L, rho_LR, vkp.pressure_tolerance);
    auto const prev_micro_liquid_density =
        computeVKPreviousMicroLiquidDensity(n_l_prev, rho_LR, local_context,
                                            vkp);
    double const rho_l_prev = n_l_prev * prev_micro_liquid_density.rho_lR;

    auto const coupled_update = solveVKNotebookMassStorageCoupledState(
        n_l_prev, rho_l_prev, prev_micro_liquid_density.rho_lR, dt, rho_LR,
        alpha_bar, mu, macro_potential, local_context, vkp);
    ASSERT_TRUE(coupled_update.converged);
    EXPECT_GT(coupled_update.n_l, 0.0);
    EXPECT_LE(coupled_update.n_l, phi + comparisonTolerance(coupled_update.n_l,
                                                            phi, 1e-10, 1e-12));
    EXPECT_GT(coupled_update.rho_lR, 0.0);

    double const active_nS = computeVKActiveMicroSolidVolumeFraction(
        coupled_update.n_l, local_context, vkp);
    auto const micro_potential = computeVanDerWaalsMicroPotential(
        coupled_update.n_l, coupled_update.rho_lR, active_nS,
        vkp.micro_solid_density_reference, vkp.hamaker_constant,
        vkp.specific_surface, vkMicroPotentialSignFactor(vkp));
    auto const exchange = computePotentialDrivenMassExchange(
        alpha_bar * rho_LR / mu,
        vkp.potential_role_mapping ==
                VKPotentialExchangeRoleMapping::NotebookRoles
            ? micro_potential.mu_lR
            : macro_potential.mu_LR,
        vkp.potential_role_mapping ==
                VKPotentialExchangeRoleMapping::NotebookRoles
            ? macro_potential.mu_LR
            : micro_potential.mu_lR);

    double const rho_l = coupled_update.n_l * coupled_update.rho_lR;
    double const volumetric_strain_rate =
        (volumetric_strain - volumetric_strain_prev) / dt;
    double const mass_residual =
        rho_l - rho_l_prev - dt * exchange.rho_l_hat -
        dt * rho_l * volumetric_strain_rate;
    auto const density = computeVKReducedMicroLiquidDensity(
        coupled_update.n_l, rho_LR, active_nS, vkp);
    double const density_residual = coupled_update.rho_lR - density.rho_lR;

    EXPECT_NEAR(mass_residual, 0.0,
                comparisonTolerance(mass_residual, 0.0, 1e-8, 1e-12));
    EXPECT_NEAR(density_residual, 0.0,
                comparisonTolerance(density_residual, 0.0, 1e-8, 1e-12));
}

TEST(RichardsMechanics, VKCoupledExchangeTangentRepresentativeStates)
{
    VKPotentialExchangeParameters vkp;
    vkp.enabled = true;
    vkp.pressure_tolerance = 0.0;
    vkp.hamaker_constant = 1e-30;
    vkp.specific_surface = 1.0;
    vkp.micro_solid_density_reference = 2650.0;
    vkp.micro_solid_volume_fraction_reference = 0.6;
    vkp.initial_micro_water_content = 0.1;
    vkp.fd_jacobian_perturbation = 1e-8;

    std::array<RepresentativeCoupledExchangeState, 3> const states = {{
        {
            .name = "legacy_placeholder_unsaturated",
            .mode = CoupledExchangeReferenceMode::legacy_placeholder,
            .p_L = -1.0e7,
            .p_L_m = -2.0e7,
            .pressure_tolerance = 0.0,
            .n_l_prev = 0.1,
            .dt = 100.0,
            .rho_LR = 1000.0,
            .drho_LR_dpL = 1.0e-7,
            .alpha_bar = 1.0e-13,
            .mu = 1.0e-3,
            .phi = 0.4,
        },
        {
            .name = "vk_full_potential_unsaturated",
            .mode = CoupledExchangeReferenceMode::full_potential_vdw,
            .p_L = -1.0e7,
            .p_L_m = 0.0,
            .pressure_tolerance = 0.0,
            .n_l_prev = 0.1,
            .dt = 100.0,
            .rho_LR = 1000.0,
            .drho_LR_dpL = 1.0e-7,
            .alpha_bar = 1.0e-13,
            .mu = 1.0e-3,
            .phi = 0.4,
        },
        {
            .name = "vk_full_potential_saturated_helper_branch",
            .mode = CoupledExchangeReferenceMode::full_potential_vdw,
            .p_L = -50.0,
            .p_L_m = 0.0,
            .pressure_tolerance = 100.0,
            .n_l_prev = 0.1,
            .dt = 100.0,
            .rho_LR = 1000.0,
            .drho_LR_dpL = 1.0e-7,
            .alpha_bar = 1.0e-13,
            .mu = 1.0e-3,
            .phi = 0.4,
        },
    }};

    for (auto const& state : states)
    {
        auto const reference_rho_L_hat =
            referenceCoupledRhoLHat(state, state.p_L, vkp);
        auto const reference_drho_L_hat_dpL =
            referenceCoupledDrhoLHatDpL(state, vkp);
        auto const production = productionCoupledExchangeData(state, vkp);

        ASSERT_TRUE(production.converged) << state.name;
        EXPECT_NEAR(production.rho_L_hat, reference_rho_L_hat,
                    comparisonTolerance(production.rho_L_hat,
                                        reference_rho_L_hat, 1e-10, 1e-18))
            << state.name;
        EXPECT_NEAR(production.drho_L_hat_dpL, reference_drho_L_hat_dpL,
                    comparisonTolerance(production.drho_L_hat_dpL,
                                        reference_drho_L_hat_dpL, 5e-5, 1e-16))
            << state.name;
    }
}
