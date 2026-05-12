// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <algorithm>
#include <cmath>
#include <Eigen/LU>
#include <cassert>
#include <limits>
#include <mutex>

#include "BaseLib/Logging.h"
#include "ComputeMicroPorosity.h"
#include "ConstitutiveRelations/ConstitutiveModels.h"
#include "ConstitutiveRelations/PotentialExchange.h"
#include "IntegrationPointData.h"
#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Utils/FormEigenTensor.h"
#include "MaterialLib/SolidModels/SelectSolidConstitutiveRelation.h"
#include "MathLib/EigenBlockMatrixView.h"
#include "MathLib/KelvinVector.h"
#include "NumLib/Fem/Interpolation.h"
#include "ProcessLib/Utils/SetOrGetIntegrationPointData.h"
#include "ProcessLib/Utils/TransposeInPlace.h"
#include "RichardsMechanicsFEM.h"

namespace ProcessLib
{
namespace RichardsMechanics
{
inline bool isPotentialExchangeEnabled(
    PotentialExchangeParameters const* const potential_exchange_parameters)
{
    return potential_exchange_parameters &&
           potential_exchange_parameters->enabled &&
           potential_exchange_parameters->mode ==
               PotentialExchangeMode::FullPotential;
}

inline bool isPotentialExchangeEnabled(
    std::optional<PotentialExchangeParameters> const&
        potential_exchange_parameters)
{
    return isPotentialExchangeEnabled(
        potential_exchange_parameters ? &*potential_exchange_parameters
                                         : nullptr);
}

inline double getPotentialPressureTolerance(
    PotentialExchangeParameters const* const potential_exchange_parameters)
{
    if (!isPotentialExchangeEnabled(potential_exchange_parameters))
    {
        return 0.0;
    }

    return potential_exchange_parameters->pressure_tolerance;
}

inline double getPotentialPressureTolerance(
    std::optional<PotentialExchangeParameters> const&
        potential_exchange_parameters)
{
    return getPotentialPressureTolerance(
        potential_exchange_parameters ? &*potential_exchange_parameters
                                         : nullptr);
}

struct PotentialExchangeUpdateData
{
    YoungLaplaceMacroPotentialData macro_potential;
    PotentialDrivenMassExchangeData exchange;

    double alpha_M_effective = 0.0;
    double mu_LR_active = 0.0;
    double mu_lR_exchange_input = 0.0;
    bool use_macro_potential_for_active_exchange = false;
    bool use_vdw_micro_potential_for_active_exchange = false;
    bool use_fd_jacobian_for_direct_macro_derivative = false;
    double fd_jacobian_perturbation = 0.0;

    // Direct macro derivative (with density dependence through rho_LR), while
    // keeping the microscale state lagged.
    double drho_L_hat_dpL_direct = 0.0;
};

inline PotentialExchangeUpdateData computePotentialExchangeUpdate(
    double const alpha_bar, double const mu, double const p_L_ip,
    double const p_L_m, double const rho_LR, double const beta_LR,
    double const pressure_tolerance = 0.0,
    bool const use_macro_potential_for_active_exchange = false,
    bool const use_vdw_micro_potential_for_active_exchange = false,
    double const mu_lR_vdw = 0.0,
    double const dmu_lR_vdw_drho_lR = 0.0,
    bool const use_custom_dmu_lR_vdw_dpL = false,
    double const dmu_lR_vdw_dpL = 0.0,
    bool const use_fd_jacobian_for_direct_macro_derivative = false,
    double const fd_jacobian_perturbation = 1e-8)
{
    if (!(mu > 0.0))
    {
        OGS_FATAL(
            "computePotentialExchangeUpdate requires mu > 0, got {:g}.",
            mu);
    }

    PotentialExchangeUpdateData out;

    // Keep the exchange coefficient scaling in mass-density units.
    out.alpha_M_effective = alpha_bar * rho_LR / mu;

    out.macro_potential =
        computeYoungLaplaceMacroPotential(p_L_ip, rho_LR, pressure_tolerance);
    out.use_macro_potential_for_active_exchange =
        use_macro_potential_for_active_exchange;
    out.use_vdw_micro_potential_for_active_exchange =
        use_vdw_micro_potential_for_active_exchange;
    out.use_fd_jacobian_for_direct_macro_derivative =
        use_fd_jacobian_for_direct_macro_derivative;
    out.fd_jacobian_perturbation = fd_jacobian_perturbation;

    out.mu_lR_exchange_input = use_vdw_micro_potential_for_active_exchange
                                   ? mu_lR_vdw
                                   : p_L_m / rho_LR;
    out.mu_LR_active = use_macro_potential_for_active_exchange
                           ? out.macro_potential.mu_LR
                           : p_L_ip / rho_LR;

    out.exchange = computePotentialDrivenMassExchange(
        out.alpha_M_effective, out.mu_LR_active, out.mu_lR_exchange_input);

    // rho_LR depends on liquid pressure in RM through beta_LR = (1/rho) drho/dp.
    double const drho_LR_dpL = rho_LR * beta_LR;

    if (use_fd_jacobian_for_direct_macro_derivative)
    {
        auto const compute_rho_L_hat = [&](double const p_L_ip_eval,
                                           double const rho_LR_eval)
        {
            auto const macro_potential_eval = computeYoungLaplaceMacroPotential(
                p_L_ip_eval, rho_LR_eval, pressure_tolerance);
            double const alpha_M_effective_eval =
                alpha_bar * rho_LR_eval / mu;
            double const mu_LR_active_eval = use_macro_potential_for_active_exchange
                                                 ? macro_potential_eval.mu_LR
                                                 : p_L_ip_eval / rho_LR_eval;
            double const mu_lR_active_eval = use_vdw_micro_potential_for_active_exchange
                                                 ? mu_lR_vdw
                                                 : p_L_m / rho_LR_eval;
            auto const exchange_eval = computePotentialDrivenMassExchange(
                alpha_M_effective_eval, mu_LR_active_eval, mu_lR_active_eval);
            return -exchange_eval.rho_l_hat;
        };

        double const h =
            fd_jacobian_perturbation * std::max(1.0, std::abs(p_L_ip));
        if (!(h > 0.0) || !std::isfinite(h))
        {
            OGS_FATAL(
                "computePotentialExchangeUpdate requires finite h > 0 for FD Jacobian, got {:g} (from fd_jacobian_perturbation={:g}, p_L_ip={:g}).",
                h, fd_jacobian_perturbation, p_L_ip);
        }

        constexpr double rho_floor = 1e-16;
        double const rho_plus = std::max(rho_floor, rho_LR + drho_LR_dpL * h);
        double const rho_minus = rho_LR - drho_LR_dpL * h;
        double const rho_L_hat_plus = compute_rho_L_hat(p_L_ip + h, rho_plus);
        if (rho_minus > rho_floor)
        {
            double const rho_L_hat_minus =
                compute_rho_L_hat(p_L_ip - h, rho_minus);
            out.drho_L_hat_dpL_direct =
                (rho_L_hat_plus - rho_L_hat_minus) / (2.0 * h);
        }
        else
        {
            double const rho_L_hat = -out.exchange.rho_l_hat;
            out.drho_L_hat_dpL_direct = (rho_L_hat_plus - rho_L_hat) / h;
        }

        return out;
    }

    // alpha_M_effective = alpha_bar * rho_LR / mu (mu dependence is lagged).
    double const dalpha_M_effective_dpL = alpha_bar / mu * drho_LR_dpL;

    double const dmu_LR_dpL = use_macro_potential_for_active_exchange
                                  ? out.macro_potential.dmu_LR_dpLR +
                                        out.macro_potential.dmu_LR_drho_LR *
                                            drho_LR_dpL
                                  : 1.0 / rho_LR -
                                        p_L_ip / (rho_LR * rho_LR) * drho_LR_dpL;

    double const dmu_lR_exchange_input_dpL =
        use_vdw_micro_potential_for_active_exchange
            ? (use_custom_dmu_lR_vdw_dpL
                   ? dmu_lR_vdw_dpL
                   : dmu_lR_vdw_drho_lR * drho_LR_dpL)
            : -p_L_m / (rho_LR * rho_LR) * drho_LR_dpL;

    double const drho_l_hat_dpL_direct =
        out.exchange.drho_l_hat_dalpha_M * dalpha_M_effective_dpL +
        out.exchange.drho_l_hat_dmu_LR * dmu_LR_dpL +
        out.exchange.drho_l_hat_dmu_lR * dmu_lR_exchange_input_dpL;

    out.drho_L_hat_dpL_direct = -drho_l_hat_dpL_direct;
    return out;
}

struct ImplicitMicroWaterContentUpdateData
{
    double n_l = 0.0;
    VanDerWaalsMicroPotentialData micro_potential;
    PotentialDrivenMassExchangeData exchange;
    bool converged = true;
};

struct CompatibilityMicroHydraulicOutputData
{
    double p_L_m = 0.0;
    double S_L_m = 0.0;
    double n_l_ref = 0.0;
    VanDerWaalsMicroPotentialData micro_potential;
};

inline double microPotentialSignFactorFromParameters(
    PotentialExchangeParameters const& potential_exchange_params)
{
    return microPotentialSignFactor(potential_exchange_params.micro_potential_convention);
}

inline CompatibilityMicroHydraulicOutputData
computeCompatibilityMicroHydraulicOutput(
    double const n_l, double const rho_LR,
    PotentialExchangeParameters const& potential_exchange_params)
{
    double const n_l_safe = std::max(1e-16, n_l);
    double const n_l_ref = std::max(
        1e-16, potential_exchange_params.initial_micro_water_content.value_or(
                   potential_exchange_params.micro_solid_volume_fraction_reference));

    auto const micro_potential = computeVanDerWaalsMicroPotential(
        n_l_safe, rho_LR, potential_exchange_params.micro_solid_volume_fraction_reference,
        potential_exchange_params.micro_solid_density_reference, potential_exchange_params.hamaker_constant,
        potential_exchange_params.specific_surface,
        microPotentialSignFactorFromParameters(potential_exchange_params),
            potential_exchange_params.vdw_augmentation_prefactor,
            potential_exchange_params.vdw_augmentation_decay_length);

    return {
        .p_L_m = -rho_LR * micro_potential.mu_lR,
        .S_L_m = n_l_safe / n_l_ref,
        .n_l_ref = n_l_ref,
        .micro_potential = micro_potential,
    };
}

struct TransportPorosityUpdateData
{
    double phi_M = 0.0;
    double phi_M_prev = 0.0;
    double phi_m = 0.0;
    double phi_m_prev = 0.0;
};

struct PotentialExchangeLocalSolveContext
{
    double phi = std::numeric_limits<double>::infinity();
    double phi_M_prev = 0.0;
    double phi_m_prev = 0.0;
    double volumetric_strain = 0.0;
    double volumetric_strain_prev = 0.0;
};

inline double boundedMicroWaterContentCeiling(
    PotentialExchangeLocalSolveContext const& local_context,
    double const n_l_floor)
{
    constexpr double porosity_upper = 1.0 - 1e-12;
    auto const compute_total_porosity_bound = [&]()
    {
        double const phi_prev_sum = std::clamp(
            std::max(0.0, local_context.phi_M_prev) +
                std::max(0.0, local_context.phi_m_prev),
            0.0, porosity_upper);
        if (std::isfinite(local_context.phi))
        {
            return std::clamp(std::max(0.0, local_context.phi), 0.0,
                              porosity_upper);
        }

        double const delta_eps_v =
            local_context.volumetric_strain - local_context.volumetric_strain_prev;
        double const denominator = 1.0 + delta_eps_v;
        if (std::isfinite(denominator) && std::abs(denominator) > 1e-12)
        {
            double const phi_from_kinematics =
                (phi_prev_sum + delta_eps_v) / denominator;
            if (std::isfinite(phi_from_kinematics))
            {
                return std::clamp(phi_from_kinematics, 0.0, porosity_upper);
            }
        }

        return phi_prev_sum;
    };

    return std::max(n_l_floor, compute_total_porosity_bound());
}

inline TransportPorosityUpdateData computeTransportPorosityUpdate(
    double const phi, double const phi_M_prev, double const phi_m_prev,
    double const n_l, double const volumetric_strain,
    double const volumetric_strain_prev,
    MacroPorosityUpdateMode const macro_porosity_update_mode)
{
    constexpr double porosity_upper = 1.0 - 1e-12;
    double const phi_prev_sum = std::clamp(
        std::max(0.0, phi_M_prev) + std::max(0.0, phi_m_prev), 0.0,
        porosity_upper);
    double const delta_eps_v = volumetric_strain - volumetric_strain_prev;
    double const denominator = 1.0 + delta_eps_v;
    double phi_safe = phi_prev_sum;
    if (std::isfinite(phi))
    {
        phi_safe = std::clamp(std::max(0.0, phi), 0.0, porosity_upper);
    }
    else if (std::isfinite(denominator) && std::abs(denominator) > 1e-12)
    {
        double const phi_from_kinematics =
            (phi_prev_sum + delta_eps_v) / denominator;
        if (std::isfinite(phi_from_kinematics))
        {
            phi_safe = std::clamp(phi_from_kinematics, 0.0, porosity_upper);
        }
    }
    double const phi_M_prev_safe = std::min(std::max(0.0, phi_M_prev), phi_safe);
    double const phi_m_prev_safe =
        std::min(std::max(0.0, phi_m_prev), std::max(0.0, phi_safe - phi_M_prev_safe));

    if (macro_porosity_update_mode ==
        MacroPorosityUpdateMode::AlgebraicSplit)
    {
        double const phi_m = std::clamp(n_l, 0.0, phi_safe);
        return {
            .phi_M = phi_safe - phi_m,
            .phi_M_prev = phi_M_prev_safe,
            .phi_m = phi_m,
            .phi_m_prev = phi_m_prev_safe,
        };
    }

    double const phi_m = std::clamp(n_l, 0.0, phi_safe);

    if (!(std::isfinite(denominator) && std::abs(denominator) > 1e-12))
    {
        static std::once_flag once;
        std::call_once(once, []
        {
            WARN(
                "[RM Phase6I] additive_macro_porosity_rate_mode porosity update encountered a near-singular denominator and fell back to algebraic_split at least once.");
        });
        return {
            .phi_M = std::max(0.0, phi_safe - phi_m),
            .phi_M_prev = phi_M_prev_safe,
            .phi_m = phi_m,
            .phi_m_prev = phi_m_prev_safe,
        };
    }

    double const phi_M_candidate =
        (phi_M_prev_safe + (1.0 - phi_m) * delta_eps_v -
         (phi_m - phi_m_prev_safe)) /
        denominator;
    double const phi_M = std::clamp(
        phi_M_candidate, 0.0, std::max(0.0, phi_safe - phi_m));

    return {
        .phi_M = phi_M,
        .phi_M_prev = phi_M_prev_safe,
        .phi_m = phi_m,
        .phi_m_prev = phi_m_prev_safe,
    };
}

inline double computeActiveMicroSolidVolumeFraction(
    double const n_l, PotentialExchangeLocalSolveContext const& local_context,
    PotentialExchangeParameters const& potential_exchange_params)
{
    if (potential_exchange_params.micro_solid_volume_fraction_mode ==
        MicroSolidVolumeFractionMode::Reference)
    {
        return std::max(1e-16, potential_exchange_params.micro_solid_volume_fraction_reference);
    }

    auto const split = computeTransportPorosityUpdate(
        local_context.phi, local_context.phi_M_prev, local_context.phi_m_prev,
        n_l, local_context.volumetric_strain,
        local_context.volumetric_strain_prev,
        potential_exchange_params.macro_porosity_update_mode);
    return std::max(1e-16, 1.0 - split.phi_M - split.phi_m);
}

inline double computePreviousMicroSolidVolumeFraction(
    PotentialExchangeLocalSolveContext const& local_context,
    PotentialExchangeParameters const& potential_exchange_params)
{
    if (potential_exchange_params.micro_solid_volume_fraction_mode ==
        MicroSolidVolumeFractionMode::Reference)
    {
        return std::max(1e-16, potential_exchange_params.micro_solid_volume_fraction_reference);
    }

    double const total_prev_porosity = std::clamp(
        std::max(0.0, local_context.phi_M_prev) +
            std::max(0.0, local_context.phi_m_prev),
        0.0, 1.0 - 1e-12);
    return std::max(1e-16, 1.0 - total_prev_porosity);
}

struct ReducedMicroLiquidDensityData
{
    double rho_lR = 0.0;
    double omega_l = 0.0;
    double drho_lR_dnl = 0.0;
    double drho_l_dn_l = 0.0;
};

inline ReducedMicroLiquidDensityData computeReducedMicroLiquidDensity(
    double const n_l, double const rho_LR, double const nS,
    PotentialExchangeParameters const& potential_exchange_params)
{
    double const n_l_safe = std::max(1e-16, n_l);
    double const nS_safe = std::max(1e-16, nS);
    double const rho_SR = std::max(1e-16, potential_exchange_params.micro_solid_density_reference);
    double const rho_l0 = std::max(1e-16, potential_exchange_params.micro_liquid_density_reference);
    double const a_rho = std::max(1e-16, potential_exchange_params.micro_liquid_density_a);
    double const b_rho = std::max(1e-16, potential_exchange_params.micro_liquid_density_b);
    double const denominator = nS_safe * rho_SR;

    auto const eval_rhs = [&](double const rho_lR)
    {
        double const omega_l =
            std::max(1e-16, n_l_safe * rho_lR / denominator);
        double const exp_term =
            std::exp(-a_rho * std::pow(omega_l, b_rho));
        return std::pair{omega_l, rho_l0 * exp_term + rho_LR};
    };

    double rho_lR = rho_LR +
                    rho_l0 *
                        std::exp(-a_rho *
                                 std::pow(std::max(1e-16, n_l_safe * rho_LR /
                                                              denominator),
                                          b_rho));
    constexpr int max_iterations = 30;
    constexpr double tolerance = 1e-14;
    bool converged = false;

    for (int iter = 0; iter < max_iterations; ++iter)
    {
        auto const [omega_l, rhs] = eval_rhs(rho_lR);
        double const residual = rho_lR - rhs;
        if (std::abs(residual) <=
            tolerance * std::max(1.0, std::abs(rho_lR)))
        {
            converged = true;
            break;
        }

        double const common =
            (rhs - rho_LR) * a_rho * b_rho *
            std::pow(omega_l, b_rho - 1.0);
        double const jacobian =
            1.0 + common * (n_l_safe / denominator);
        if (!(std::isfinite(jacobian) && std::abs(jacobian) > 1e-20))
        {
            break;
        }

        double const rho_candidate =
            std::max(1e-16, rho_lR - residual / jacobian);
        if (std::abs(rho_candidate - rho_lR) <=
            tolerance * std::max(1.0, std::abs(rho_lR)))
        {
            rho_lR = rho_candidate;
            converged = true;
            break;
        }
        rho_lR = rho_candidate;
    }

    if (!converged)
    {
        static std::once_flag once;
        std::call_once(once, []
        {
            WARN(
                "[RM Phase6K] reduced microscale liquid-density EOS did not converge at least once; using the last Newton iterate.");
        });
    }

    auto const [omega_l, rhs] = eval_rhs(rho_lR);
    (void)rhs;
    double const common =
        (rho_lR - rho_LR) * a_rho * b_rho *
        std::pow(omega_l, b_rho - 1.0);
    double const dg_drho =
        1.0 + common * (n_l_safe / denominator);
    double const dg_dn = common * (rho_lR / denominator);
    double const drho_lR_dnl =
        (std::isfinite(dg_drho) && std::abs(dg_drho) > 1e-20)
            ? -dg_dn / dg_drho
            : 0.0;

    return {
        .rho_lR = rho_lR,
        .omega_l = omega_l,
        .drho_lR_dnl = drho_lR_dnl,
        .drho_l_dn_l = rho_lR + n_l_safe * drho_lR_dnl,
    };
}

inline ReducedMicroLiquidDensityData computeActiveMicroLiquidDensity(
    double const n_l, double const rho_LR,
    PotentialExchangeLocalSolveContext const& local_context,
    PotentialExchangeParameters const& potential_exchange_params)
{
    double const active_nS =
        computeActiveMicroSolidVolumeFraction(n_l, local_context, potential_exchange_params);
    return computeReducedMicroLiquidDensity(n_l, rho_LR, active_nS, potential_exchange_params);
}

inline ReducedMicroLiquidDensityData computePreviousMicroLiquidDensity(
    double const n_l_prev, double const rho_LR,
    PotentialExchangeLocalSolveContext const& local_context,
    PotentialExchangeParameters const& potential_exchange_params)
{
    double const previous_nS =
        computePreviousMicroSolidVolumeFraction(local_context, potential_exchange_params);
    return computeReducedMicroLiquidDensity(n_l_prev, rho_LR, previous_nS,
                                              potential_exchange_params);
}

struct MicroMacroMassStorageCoupledSolveData
{
    double n_l = 0.0;
    double rho_lR = 0.0;
    double phi_m = 0.0;
    double phi_M = 0.0;
    double p_L_m = 0.0;
    double S_L_m = 0.0;
    VanDerWaalsMicroPotentialData micro_potential;
    PotentialDrivenMassExchangeData exchange;
    bool converged = true;
};

inline MicroMacroMassStorageCoupledSolveData
solveReferenceMassStoragePredictorState(
    double const n_l_prev, double const rho_l_prev, double const rho_lR_prev,
    double const dt, double const rho_LR, double const alpha_bar,
    double const mu, YoungLaplaceMacroPotentialData const& macro_potential,
    PotentialExchangeLocalSolveContext const& local_context,
    PotentialExchangeParameters const& potential_exchange_params)
{
    constexpr double n_l_floor = 1e-16;
    constexpr double rho_floor = 1e-16;
    double const dt_safe = std::isfinite(dt) && dt > 0.0 ? dt : 0.0;
    double const alpha_M_effective = alpha_bar * rho_LR / mu;
    double const volumetric_strain_rate =
        dt_safe > 0.0
            ? (local_context.volumetric_strain -
               local_context.volumetric_strain_prev) /
                  dt_safe
            : 0.0;
    double const n_l_ceiling =
        boundedMicroWaterContentCeiling(local_context, n_l_floor);

    auto evaluate = [&](double const n_l)
    {
        double const active_nS = computeActiveMicroSolidVolumeFraction(
            n_l, local_context, potential_exchange_params);
        auto const micro_liquid_density = computeReducedMicroLiquidDensity(
            n_l, rho_LR, active_nS, potential_exchange_params);
        auto const micro_potential = computeVanDerWaalsMicroPotential(
            n_l, micro_liquid_density.rho_lR, active_nS,
            potential_exchange_params.micro_solid_density_reference, potential_exchange_params.hamaker_constant,
            potential_exchange_params.specific_surface,
            microPotentialSignFactorFromParameters(potential_exchange_params),
            potential_exchange_params.vdw_augmentation_prefactor,
            potential_exchange_params.vdw_augmentation_decay_length);
        double const mu_LR_active = macro_potential.mu_LR;
        double const mu_lR_active = micro_potential.mu_lR;
        auto const exchange = computePotentialDrivenMassExchange(
            alpha_M_effective, mu_LR_active, mu_lR_active);
        double const rho_l = n_l * micro_liquid_density.rho_lR;
        double const residual = rho_l - rho_l_prev -
                                dt_safe * exchange.rho_l_hat -
                                dt_safe * rho_l * volumetric_strain_rate;
        return std::tuple{residual, micro_potential, exchange,
                          micro_liquid_density};
    };

    MicroMacroMassStorageCoupledSolveData out;
    if (dt_safe <= 0.0)
    {
        out.n_l = std::clamp(n_l_prev, n_l_floor, n_l_ceiling);
        out.rho_lR = std::max(rho_floor, rho_lR_prev);
        auto const [residual, micro_potential, exchange, micro_density] =
            evaluate(out.n_l);
        (void)residual;
        (void)micro_density;
        out.micro_potential = micro_potential;
        out.exchange = exchange;
        return out;
    }

    double n_l = std::clamp(n_l_prev, n_l_floor, n_l_ceiling);
    constexpr int max_iterations = 40;
    constexpr double residual_tolerance = 1e-14;
    constexpr double increment_tolerance = 1e-14;

    for (int iter = 0; iter < max_iterations; ++iter)
    {
        auto const [residual, micro_potential, exchange, micro_density] =
            evaluate(n_l);
        if (std::abs(residual) <=
            residual_tolerance * std::max(1.0, std::abs(rho_l_prev)))
        {
            out.n_l = n_l;
            out.rho_lR = micro_density.rho_lR;
            out.micro_potential = micro_potential;
            out.exchange = exchange;
            return out;
        }

        double const drho_l_hat_dn_l =
            exchange.drho_l_hat_dmu_lR * micro_potential.dmu_lR_dnl;
        double const jacobian = micro_density.drho_l_dn_l -
                                dt_safe * drho_l_hat_dn_l -
                                dt_safe * micro_density.drho_l_dn_l *
                                    volumetric_strain_rate;
        if (!(std::isfinite(jacobian) && std::abs(jacobian) > 1e-20))
        {
            break;
        }

        double delta_n_l = -residual / jacobian;
        double n_l_candidate =
            std::clamp(n_l + delta_n_l, n_l_floor, n_l_ceiling);
        auto const [candidate_residual_initial, candidate_micro_potential,
                    candidate_exchange, candidate_micro_density] =
            evaluate(n_l_candidate);
        double candidate_residual = candidate_residual_initial;
        int backtracking_steps = 0;
        while (std::abs(candidate_residual) > std::abs(residual) &&
               backtracking_steps < 12)
        {
            delta_n_l *= 0.5;
            n_l_candidate =
                std::clamp(n_l + delta_n_l, n_l_floor, n_l_ceiling);
            auto const [retry_residual, retry_micro_potential,
                        retry_exchange, retry_micro_density] =
                evaluate(n_l_candidate);
            (void)retry_micro_potential;
            (void)retry_exchange;
            (void)retry_micro_density;
            candidate_residual = retry_residual;
            ++backtracking_steps;
        }

        if (std::abs(n_l_candidate - n_l) <=
            increment_tolerance * std::max(1.0, std::abs(n_l)))
        {
            out.n_l = n_l_candidate;
            auto const [final_residual, final_micro_potential, final_exchange,
                        final_micro_density] =
                evaluate(n_l_candidate);
            (void)final_residual;
            out.rho_lR = final_micro_density.rho_lR;
            out.micro_potential = final_micro_potential;
            out.exchange = final_exchange;
            out.converged = true;
            return out;
        }

        n_l = n_l_candidate;
        out.n_l = n_l;
        out.rho_lR = micro_density.rho_lR;
        out.micro_potential = micro_potential;
        out.exchange = exchange;
    }

    auto const [residual, micro_potential, exchange, micro_density] =
        evaluate(n_l);
    (void)residual;
    out.n_l = n_l;
    out.rho_lR = micro_density.rho_lR;
    out.micro_potential = micro_potential;
    out.exchange = exchange;
    out.converged = false;
    return out;
}

inline MicroMacroMassStorageCoupledSolveData
solveReferenceMassStorageCoupledState(
    double const n_l_prev, double const rho_l_prev, double const rho_lR_prev,
    double const dt, double const rho_LR, double const alpha_bar,
    double const mu, YoungLaplaceMacroPotentialData const& macro_potential,
    PotentialExchangeLocalSolveContext const& local_context,
    PotentialExchangeParameters const& potential_exchange_params)
{
    constexpr double n_l_floor = 1e-16;
    constexpr double rho_floor = 1e-16;
    double const dt_safe = std::isfinite(dt) && dt > 0.0 ? dt : 0.0;
    double const alpha_M_effective = alpha_bar * rho_LR / mu;
    double const volumetric_strain_rate =
        dt_safe > 0.0
            ? (local_context.volumetric_strain -
               local_context.volumetric_strain_prev) /
                  dt_safe
            : 0.0;
    double const n_l_ceiling =
        boundedMicroWaterContentCeiling(local_context, n_l_floor);

    auto evaluate = [&](double const n_l, double const rho_lR)
    {
        double const active_nS = computeActiveMicroSolidVolumeFraction(
            n_l, local_context, potential_exchange_params);
        auto const micro_potential = computeVanDerWaalsMicroPotential(
            n_l, rho_lR, active_nS, potential_exchange_params.micro_solid_density_reference,
            potential_exchange_params.hamaker_constant, potential_exchange_params.specific_surface,
            microPotentialSignFactorFromParameters(potential_exchange_params),
            potential_exchange_params.vdw_augmentation_prefactor,
            potential_exchange_params.vdw_augmentation_decay_length);
        double const mu_LR_active = macro_potential.mu_LR;
        double const mu_lR_active = micro_potential.mu_lR;
        auto const exchange = computePotentialDrivenMassExchange(
            alpha_M_effective, mu_LR_active, mu_lR_active);
        double const rho_l = n_l * rho_lR;
        double const mass_residual = rho_l - rho_l_prev -
                                     dt_safe * exchange.rho_l_hat -
                                     dt_safe * rho_l * volumetric_strain_rate;
        auto const density = computeReducedMicroLiquidDensity(
            n_l, rho_LR, active_nS, potential_exchange_params);
        double const density_residual = rho_lR - density.rho_lR;
        return std::tuple{mass_residual, density_residual, micro_potential,
                          exchange};
    };

    auto const predictor = solveReferenceMassStoragePredictorState(
        n_l_prev, rho_l_prev, rho_lR_prev, dt, rho_LR, alpha_bar, mu,
        macro_potential, local_context, potential_exchange_params);
    if (!predictor.converged)
    {
        return predictor;
    }

    MicroMacroMassStorageCoupledSolveData out = predictor;
    if (dt_safe <= 0.0)
    {
        return out;
    }

    double n_l = std::clamp(n_l_prev, n_l_floor, n_l_ceiling);
    double rho_lR = std::max(rho_floor, rho_lR_prev);
    constexpr int max_iterations = 60;
    constexpr double residual_tolerance = 1e-10;
    constexpr double increment_tolerance = 1e-10;

    for (int iter = 0; iter < max_iterations; ++iter)
    {
        auto const [mass_residual, density_residual, micro_potential, exchange] =
            evaluate(n_l, rho_lR);

        double const residual_norm =
            std::abs(mass_residual) / std::max(1.0, std::abs(rho_l_prev)) +
            std::abs(density_residual) / std::max(1.0, std::abs(rho_lR));
        if (residual_norm <= residual_tolerance)
        {
            out.n_l = n_l;
            out.rho_lR = rho_lR;
            out.micro_potential = micro_potential;
            out.exchange = exchange;
            out.converged = true;
            return out;
        }

        double const h_n = 1e-8 * std::max(1.0, std::abs(n_l));
        double const h_rho = 1e-8 * std::max(1.0, std::abs(rho_lR));
        auto const [r1_n_plus, r2_n_plus] =
            [&]() {
                auto const [r1, r2, _, __] = evaluate(n_l + h_n, rho_lR);
                (void)_;
                (void)__;
                return std::pair{r1, r2};
            }();
        auto const [r1_n_minus, r2_n_minus] =
            [&]() {
                auto const [r1, r2, _, __] =
                    evaluate(std::max(n_l_floor, n_l - h_n), rho_lR);
                (void)_;
                (void)__;
                return std::pair{r1, r2};
            }();
        auto const [r1_rho_plus, r2_rho_plus] =
            [&]() {
                auto const [r1, r2, _, __] = evaluate(n_l, rho_lR + h_rho);
                (void)_;
                (void)__;
                return std::pair{r1, r2};
            }();
        auto const [r1_rho_minus, r2_rho_minus] =
            [&]() {
                auto const [r1, r2, _, __] =
                    evaluate(n_l, std::max(rho_floor, rho_lR - h_rho));
                (void)_;
                (void)__;
                return std::pair{r1, r2};
            }();

        double const denom_n = (n_l + h_n) - std::max(n_l_floor, n_l - h_n);
        double const denom_rho =
            (rho_lR + h_rho) - std::max(rho_floor, rho_lR - h_rho);
        if (!(denom_n > 0.0 && denom_rho > 0.0))
        {
            break;
        }

        double const J11 = (r1_n_plus - r1_n_minus) / denom_n;
        double const J21 = (r2_n_plus - r2_n_minus) / denom_n;
        double const J12 = (r1_rho_plus - r1_rho_minus) / denom_rho;
        double const J22 = (r2_rho_plus - r2_rho_minus) / denom_rho;

        double const det = J11 * J22 - J12 * J21;
        if (!(std::isfinite(det) && std::abs(det) > 1e-24))
        {
            break;
        }

        double const delta_n = (-mass_residual * J22 +
                                density_residual * J12) /
                               det;
        double const delta_rho = (J21 * mass_residual -
                                  J11 * density_residual) /
                                 det;

        double step_scale = 1.0;
        bool accepted = false;
        for (int backtrack = 0; backtrack < 12; ++backtrack)
        {
            double const n_candidate =
                std::clamp(n_l + step_scale * delta_n, n_l_floor, n_l_ceiling);
            double const rho_candidate =
                std::max(rho_floor, rho_lR + step_scale * delta_rho);
            auto const [cand_mass_residual, cand_density_residual,
                        cand_micro_potential, cand_exchange] =
                evaluate(n_candidate, rho_candidate);
            double const current_norm =
                std::abs(mass_residual) /
                    std::max(1.0, std::abs(rho_l_prev)) +
                std::abs(density_residual) / std::max(1.0, std::abs(rho_lR));
            double const candidate_norm =
                std::abs(cand_mass_residual) /
                    std::max(1.0, std::abs(rho_l_prev)) +
                std::abs(cand_density_residual) /
                    std::max(1.0, std::abs(rho_candidate));
            if (candidate_norm <= current_norm || step_scale < 1e-3)
            {
                n_l = n_candidate;
                rho_lR = rho_candidate;
                out.n_l = n_l;
                out.rho_lR = rho_lR;
                out.micro_potential = cand_micro_potential;
                out.exchange = cand_exchange;
                accepted = true;
                break;
            }
            step_scale *= 0.5;
        }

        if (!accepted)
        {
            break;
        }

        if (std::abs(step_scale * delta_n) <=
                increment_tolerance * std::max(1.0, std::abs(n_l)) &&
            std::abs(step_scale * delta_rho) <=
                increment_tolerance * std::max(1.0, std::abs(rho_lR)))
        {
            out.converged = true;
            return out;
        }
    }

    auto const [mass_residual, density_residual, micro_potential, exchange] =
        evaluate(n_l, rho_lR);
    (void)mass_residual;
    (void)density_residual;
    out.n_l = n_l;
    out.rho_lR = rho_lR;
    out.micro_potential = micro_potential;
    out.exchange = exchange;
    out.converged = false;
    return out.converged ? out : predictor;
}

template <int DisplacementDim>
inline void applyReferenceMassStorageLocalState(
    StatefulData<DisplacementDim>& state_current,
    StatefulDataPrev<DisplacementDim> const& state_previous,
    MPL::VariableArray& variables, MPL::VariableArray& variables_prev,
    double const rho_LR, PotentialExchangeLocalSolveContext const& local_context,
    PotentialExchangeParameters const& potential_exchange_params,
    MicroMacroMassStorageCoupledSolveData const& coupled_update)
{
    auto const transport_porosity_update = computeTransportPorosityUpdate(
        local_context.phi, local_context.phi_M_prev, local_context.phi_m_prev,
        coupled_update.n_l, local_context.volumetric_strain,
        local_context.volumetric_strain_prev,
        potential_exchange_params.macro_porosity_update_mode);

    auto const compatibility_output =
        computeCompatibilityMicroHydraulicOutput(
            coupled_update.n_l, rho_LR, potential_exchange_params);

    auto& n_l = std::get<MicroWaterContent>(state_current);
    *n_l = coupled_update.n_l;

    auto& rho_lR = std::get<MicroLiquidDensity>(state_current);
    *rho_lR = coupled_update.rho_lR;

    auto& micro_porosity = std::get<MicroPorosity>(state_current);
    *micro_porosity = transport_porosity_update.phi_m;

    auto& transport_porosity =
        std::get<ProcessLib::ThermoRichardsMechanics::TransportPorosityData>(
            state_current)
            .phi;
    transport_porosity = transport_porosity_update.phi_M;
    variables.transport_porosity = transport_porosity_update.phi_M;
    variables_prev.transport_porosity = transport_porosity_update.phi_M_prev;

    auto& porosity =
        std::get<ProcessLib::ThermoRichardsMechanics::PorosityData>(state_current).phi;
    porosity = transport_porosity_update.phi_M + transport_porosity_update.phi_m;
    variables.porosity = porosity;
    variables_prev.porosity =
        transport_porosity_update.phi_M_prev + transport_porosity_update.phi_m_prev;

    auto& p_L_m = std::get<MicroPressure>(state_current);
    auto& S_L_m = std::get<MicroSaturation>(state_current);
    auto& rho_l_hat = std::get<MicroExchangeSource>(state_current);
    *p_L_m = compatibility_output.p_L_m;
    *S_L_m = compatibility_output.S_L_m;
    rho_l_hat = MicroExchangeSource{coupled_update.exchange.rho_l_hat};

    (void)state_previous;
}

inline VanDerWaalsMicroPotentialData computeActiveMicroPotential(
    double const n_l, double const rho_lR,
    PotentialExchangeLocalSolveContext const& local_context,
    PotentialExchangeParameters const& potential_exchange_params)
{
    double const active_nS =
        computeActiveMicroSolidVolumeFraction(n_l, local_context, potential_exchange_params);
    double const rho_lR_effective =
        potential_exchange_params.local_nonlinear_solve_mode ==
                LocalNonlinearSolveMode::ScalarReferenceMassStorage
            ? computeReducedMicroLiquidDensity(n_l, rho_lR, active_nS, potential_exchange_params)
                  .rho_lR
            : rho_lR;
    return computeVanDerWaalsMicroPotential(
        n_l, rho_lR_effective, active_nS, potential_exchange_params.micro_solid_density_reference,
        potential_exchange_params.hamaker_constant, potential_exchange_params.specific_surface,
        microPotentialSignFactorFromParameters(potential_exchange_params),
            potential_exchange_params.vdw_augmentation_prefactor,
            potential_exchange_params.vdw_augmentation_decay_length);
}

inline CompatibilityMicroHydraulicOutputData
computeCompatibilityMicroHydraulicOutput(
    double const n_l, double const rho_LR,
    PotentialExchangeLocalSolveContext const& local_context,
    PotentialExchangeParameters const& potential_exchange_params)
{
    double const n_l_safe = std::max(1e-16, n_l);
    double const n_l_ref = std::max(
        1e-16, potential_exchange_params.initial_micro_water_content.value_or(
                   potential_exchange_params.micro_solid_volume_fraction_reference));

    auto const micro_potential =
        computeActiveMicroPotential(n_l_safe, rho_LR, local_context, potential_exchange_params);

    return {
        .p_L_m = -rho_LR * micro_potential.mu_lR,
        .S_L_m = n_l_safe / n_l_ref,
        .n_l_ref = n_l_ref,
        .micro_potential = micro_potential,
    };
}

inline ImplicitMicroWaterContentUpdateData solveImplicitMicroWaterContent(
    double const n_l_prev, double const dt, double const rho_LR,
    double const alpha_bar, double const mu,
    YoungLaplaceMacroPotentialData const& macro_potential,
    PotentialExchangeLocalSolveContext const& local_context,
    PotentialExchangeParameters const& potential_exchange_params)
{
    constexpr double n_l_floor = 1e-16;
    double const dt_safe = std::isfinite(dt) && dt > 0.0 ? dt : 0.0;
    double const alpha_M_effective = alpha_bar * rho_LR / mu;
    bool const use_microstate_storage_mode =
        potential_exchange_params.local_nonlinear_solve_mode !=
        LocalNonlinearSolveMode::ScalarExchange;
    bool const use_mass_storage =
        potential_exchange_params.local_nonlinear_solve_mode ==
        LocalNonlinearSolveMode::ScalarReferenceMassStorage;
    double const volumetric_strain_rate =
        dt_safe > 0.0
            ? (local_context.volumetric_strain -
               local_context.volumetric_strain_prev) /
                  dt_safe
            : 0.0;
    double const n_l_ceiling =
        use_microstate_storage_mode
            ? boundedMicroWaterContentCeiling(local_context, n_l_floor)
            : std::max(n_l_floor, 1.0);

    auto eval_at = [&](double const n_l)
    {
        auto const micro_potential =
            computeActiveMicroPotential(n_l, rho_LR, local_context, potential_exchange_params);
        double const mu_LR_active = macro_potential.mu_LR;
        double const mu_lR_active = micro_potential.mu_lR;
        auto const exchange = computePotentialDrivenMassExchange(
            alpha_M_effective, mu_LR_active, mu_lR_active);
        auto const micro_liquid_density =
            use_mass_storage
                ? std::optional<ReducedMicroLiquidDensityData>{
                      computeActiveMicroLiquidDensity(
                          n_l, rho_LR, local_context, potential_exchange_params)}
                : std::nullopt;
        return std::tuple{micro_potential, exchange, micro_liquid_density};
    };

    auto const prev_micro_liquid_density =
        use_mass_storage
            ? std::optional<ReducedMicroLiquidDensityData>{
                  computePreviousMicroLiquidDensity(n_l_prev, rho_LR,
                                                      local_context, potential_exchange_params)}
            : std::nullopt;
    double const rho_l_prev =
        prev_micro_liquid_density
            ? n_l_prev * prev_micro_liquid_density->rho_lR
            : 0.0;

    ImplicitMicroWaterContentUpdateData out;
    if (dt_safe <= 0.0)
    {
        out.n_l = std::clamp(n_l_prev, n_l_floor, n_l_ceiling);
        auto const [micro_potential, exchange, micro_liquid_density] =
            eval_at(out.n_l);
        (void)micro_liquid_density;
        out.micro_potential = micro_potential;
        out.exchange = exchange;
        return out;
    }

    if (use_mass_storage)
    {
        auto const coupled_update = solveReferenceMassStorageCoupledState(
            n_l_prev, rho_l_prev,
            prev_micro_liquid_density ? prev_micro_liquid_density->rho_lR
                                      : rho_LR,
            dt_safe, rho_LR, alpha_bar, mu, macro_potential,
            local_context, potential_exchange_params);
        out.n_l = coupled_update.n_l;
        out.micro_potential = coupled_update.micro_potential;
        out.exchange = coupled_update.exchange;
        out.converged = coupled_update.converged;
        return out;
    }

    double n_l = std::clamp(n_l_prev, n_l_floor, n_l_ceiling);
    constexpr int max_iterations = 25;
    constexpr double residual_tolerance = 1e-12;
    constexpr double increment_tolerance = 1e-12;
    bool converged = false;

    for (int iter = 0; iter < max_iterations; ++iter)
    {
        auto const [micro_potential, exchange, micro_liquid_density] =
            eval_at(n_l);
        double residual = 0.0;
        double jacobian = 0.0;
        double const drho_l_hat_dn_l =
            exchange.drho_l_hat_dmu_lR * micro_potential.dmu_lR_dnl;
        if (use_mass_storage)
        {
            double const rho_l =
                n_l * micro_liquid_density->rho_lR;
            residual =
                rho_l - rho_l_prev - dt_safe * exchange.rho_l_hat;
            residual -= dt_safe * rho_l * volumetric_strain_rate;

            jacobian = micro_liquid_density->drho_l_dn_l -
                       dt_safe * drho_l_hat_dn_l;
            jacobian -=
                dt_safe * micro_liquid_density->drho_l_dn_l *
                volumetric_strain_rate;
        }
        else
        {
            residual =
                n_l - n_l_prev - dt_safe * exchange.rho_l_hat / rho_LR;
            if (use_microstate_storage_mode)
            {
                residual -= dt_safe * n_l * volumetric_strain_rate;
            }

            jacobian = 1.0 - dt_safe * drho_l_hat_dn_l / rho_LR;
            if (use_microstate_storage_mode)
            {
                jacobian -= dt_safe * volumetric_strain_rate;
            }
        }

        if (std::abs(residual) <=
            residual_tolerance * std::max(1.0, std::abs(n_l_prev)))
        {
            converged = true;
            out.n_l = n_l;
            out.micro_potential = micro_potential;
            out.exchange = exchange;
            break;
        }

        if (!(std::isfinite(jacobian) && std::abs(jacobian) > 1e-20))
        {
            break;
        }

        double const delta_n_l = -residual / jacobian;
        double const n_l_candidate =
            std::clamp(n_l + delta_n_l, n_l_floor, n_l_ceiling);
        if (std::abs(n_l_candidate - n_l) <=
            increment_tolerance * std::max(1.0, std::abs(n_l)))
        {
            auto const [micro_potential_candidate, exchange_candidate,
                        micro_density_candidate] =
                eval_at(n_l_candidate);
            (void)micro_density_candidate;
            out.n_l = n_l_candidate;
            out.micro_potential = micro_potential_candidate;
            out.exchange = exchange_candidate;
            converged = true;
            break;
        }

        n_l = n_l_candidate;
    }

    if (!converged)
    {
        // Fallback to explicit update if local scalar Newton does not converge.
        auto const [micro_potential_prev, exchange_prev, micro_density_prev] =
            eval_at(std::clamp(n_l_prev, n_l_floor, n_l_ceiling));
        double explicit_increment = 0.0;
        if (use_mass_storage)
        {
            double const rho_l_prev_fallback =
                std::clamp(n_l_prev, n_l_floor, n_l_ceiling) *
                micro_density_prev->rho_lR;
            explicit_increment =
                dt_safe * exchange_prev.rho_l_hat /
                std::max(1e-16, micro_density_prev->rho_lR);
            explicit_increment +=
                dt_safe * rho_l_prev_fallback * volumetric_strain_rate /
                std::max(1e-16, micro_density_prev->rho_lR);
        }
        else
        {
            explicit_increment =
                dt_safe * exchange_prev.rho_l_hat / rho_LR;
            if (use_microstate_storage_mode)
            {
                explicit_increment +=
                    dt_safe * std::clamp(n_l_prev, n_l_floor, n_l_ceiling) *
                    volumetric_strain_rate;
            }
        }
        out.n_l = std::clamp(n_l_prev + explicit_increment, n_l_floor,
                             n_l_ceiling);
        auto const [micro_potential_fallback, exchange_fallback,
                    micro_density_fallback] =
            eval_at(out.n_l);
        (void)micro_density_fallback;
        out.micro_potential = micro_potential_fallback;
        out.exchange = exchange_fallback;
        out.converged = false;

        static std::once_flag once;
        std::call_once(once, []
        {
            WARN(
                "[RM Phase2C] local implicit n_l solve did not converge at least once; falling back to explicit n_l update for robustness.");
        });
        return out;
    }

    out.converged = true;
    return out;
}

inline double computeImplicitNlDpL(
    double const n_l_prev, double const p_L_ip, double const dt,
    double const rho_LR, double const drho_LR_dpL,
    double const alpha_bar, double const mu,
    YoungLaplaceMacroPotentialData const& macro_potential,
    VanDerWaalsMicroPotentialData const& micro_potential,
    PotentialDrivenMassExchangeData const& exchange,
    PotentialExchangeLocalSolveContext const& local_context,
    PotentialExchangeParameters const& potential_exchange_params)
{
    double const dt_safe = std::isfinite(dt) && dt > 0.0 ? dt : 0.0;
    if (dt_safe <= 0.0)
    {
        return 0.0;
    }

    if (potential_exchange_params.local_nonlinear_solve_mode ==
        LocalNonlinearSolveMode::ScalarReferenceMassStorage)
    {
        constexpr double rho_floor = 1e-16;
        double const perturbation =
            std::max(potential_exchange_params.local_jacobian_perturbation, 1e-8) *
            std::max(1.0, std::abs(p_L_ip));

        auto const eval_at = [&](double const p_L_eval,
                                 double const rho_LR_eval)
        {
            auto const macro_potential_eval = computeYoungLaplaceMacroPotential(
                p_L_eval, rho_LR_eval, potential_exchange_params.pressure_tolerance);
            auto const n_l_update_eval = solveImplicitMicroWaterContent(
                n_l_prev, dt, rho_LR_eval, alpha_bar, mu,
                macro_potential_eval, local_context, potential_exchange_params);
            return n_l_update_eval.n_l;
        };

        double const rho_plus =
            std::max(rho_floor, rho_LR + drho_LR_dpL * perturbation);
        double const n_l_plus = eval_at(p_L_ip + perturbation, rho_plus);
        double const rho_minus = rho_LR - drho_LR_dpL * perturbation;
        if (rho_minus > rho_floor)
        {
            double const n_l_minus =
                eval_at(p_L_ip - perturbation, rho_minus);
            return (n_l_plus - n_l_minus) / (2.0 * perturbation);
        }
        double const n_l_center = eval_at(p_L_ip, rho_LR);
        return (n_l_plus - n_l_center) / perturbation;
    }

    double const dalpha_M_effective_dpL = alpha_bar / mu * drho_LR_dpL;
    double const dmu_first_dpL_fixed_n =
        macro_potential.dmu_LR_dpLR +
        macro_potential.dmu_LR_drho_LR * drho_LR_dpL;
    double const dmu_second_dpL_fixed_n = micro_potential.dmu_lR_drho_lR *
                                          drho_LR_dpL;

    double const drho_l_hat_dpL_fixed_n =
        exchange.drho_l_hat_dalpha_M * dalpha_M_effective_dpL +
        exchange.drho_l_hat_dmu_LR * dmu_first_dpL_fixed_n +
        exchange.drho_l_hat_dmu_lR * dmu_second_dpL_fixed_n;
    double const drho_l_hat_dn_l =
        exchange.drho_l_hat_dmu_lR * micro_potential.dmu_lR_dnl;

    double dr_dn_l = 1.0 - dt_safe * drho_l_hat_dn_l / rho_LR;
    if (potential_exchange_params.local_nonlinear_solve_mode ==
        LocalNonlinearSolveMode::ScalarReferenceStorage)
    {
        double const volumetric_strain_rate =
            (local_context.volumetric_strain -
             local_context.volumetric_strain_prev) /
            dt_safe;
        dr_dn_l -= dt_safe * volumetric_strain_rate;
    }
    if (!(std::isfinite(dr_dn_l) && std::abs(dr_dn_l) > 1e-20))
    {
        return 0.0;
    }

    double const dr_dp_l =
        -dt_safe * (drho_l_hat_dpL_fixed_n / rho_LR -
                    exchange.rho_l_hat / (rho_LR * rho_LR) * drho_LR_dpL);
    return -dr_dp_l / dr_dn_l;
}

struct LocalJacobianDiagnosticData
{
    double fd_dn_l_dpL = 0.0;
    double fd_drho_L_hat_dpL = 0.0;
    double perturbation = 0.0;
};

inline LocalJacobianDiagnosticData computeLocalJacobianDiagnosticData(
    double const n_l_prev, double const p_L_ip, double const dt,
    double const rho_LR, double const drho_LR_dpL, double const alpha_bar,
    double const mu, double const pressure_tolerance,
    PotentialExchangeLocalSolveContext const& local_context,
    PotentialExchangeParameters const& potential_exchange_params)
{
    constexpr double rho_floor = 1e-16;
    double const perturbation =
        potential_exchange_params.local_jacobian_perturbation * std::max(1.0, std::abs(p_L_ip));
    if (!(perturbation > 0.0) || !std::isfinite(perturbation))
    {
        OGS_FATAL(
            "DSM local Jacobian diagnostic requires finite h > 0, got {:g} "
            "(from local_jacobian_perturbation={:g}, p_L_ip={:g}).",
            perturbation, potential_exchange_params.local_jacobian_perturbation, p_L_ip);
    }

    auto const eval_at = [&](double const p_L_eval, double const rho_LR_eval)
    {
        auto const macro_potential_eval = computeYoungLaplaceMacroPotential(
            p_L_eval, rho_LR_eval, pressure_tolerance);
        auto const n_l_update_eval = solveImplicitMicroWaterContent(
            n_l_prev, dt, rho_LR_eval, alpha_bar, mu, macro_potential_eval,
            local_context, potential_exchange_params);
        return std::pair{n_l_update_eval.n_l, -n_l_update_eval.exchange.rho_l_hat};
    };

    double const rho_plus =
        std::max(rho_floor, rho_LR + drho_LR_dpL * perturbation);
    auto const [n_l_plus, rho_L_hat_plus] =
        eval_at(p_L_ip + perturbation, rho_plus);

    double const rho_minus = rho_LR - drho_LR_dpL * perturbation;
    if (rho_minus > rho_floor)
    {
        auto const [n_l_minus, rho_L_hat_minus] =
            eval_at(p_L_ip - perturbation, rho_minus);
        return {
            .fd_dn_l_dpL = (n_l_plus - n_l_minus) / (2.0 * perturbation),
            .fd_drho_L_hat_dpL =
                (rho_L_hat_plus - rho_L_hat_minus) / (2.0 * perturbation),
            .perturbation = perturbation,
        };
    }

    auto const [n_l_center, rho_L_hat_center] = eval_at(p_L_ip, rho_LR);
    return {
        .fd_dn_l_dpL = (n_l_plus - n_l_center) / perturbation,
        .fd_drho_L_hat_dpL =
            (rho_L_hat_plus - rho_L_hat_center) / perturbation,
        .perturbation = perturbation,
    };
}

inline void maybeLogLocalJacobianDiagnostic(
    double const p_L_ip, double const n_l_prev, double const n_l,
    double const analytic_dn_l_dpL, double const analytic_drho_L_hat_dpL,
    LocalJacobianDiagnosticData const& fd_data,
    PotentialExchangeParameters const& potential_exchange_params)
{
    static std::once_flag once;
    std::call_once(once, [=]()
    {
        auto const relative_error = [](double const analytic, double const fd)
        {
            double const scale =
                std::max({1.0, std::abs(analytic), std::abs(fd)});
            return std::abs(analytic - fd) / scale;
        };

        double const rel_dn_l =
            relative_error(analytic_dn_l_dpL, fd_data.fd_dn_l_dpL);
        double const rel_drho = relative_error(
            analytic_drho_L_hat_dpL, fd_data.fd_drho_L_hat_dpL);
        bool const mismatch =
            rel_dn_l > potential_exchange_params.local_jacobian_relative_tolerance ||
            rel_drho > potential_exchange_params.local_jacobian_relative_tolerance;

        auto const* const level_prefix = mismatch ? "[RM Phase3D]" : "[RM Phase3D]";
        if (mismatch)
        {
            WARN(
                "{} local Jacobian diagnostic: pL_ip={} Pa, n_l_prev={}, n_l={}, h={}, analytic dn_l/dpL={}, FD dn_l/dpL={}, rel_err_nl={}, analytic drho_L_hat/dpL={}, FD drho_L_hat/dpL={}, rel_err_exchange={}, tolerance={}.",
                level_prefix, p_L_ip, n_l_prev, n_l, fd_data.perturbation,
                analytic_dn_l_dpL, fd_data.fd_dn_l_dpL, rel_dn_l,
                analytic_drho_L_hat_dpL, fd_data.fd_drho_L_hat_dpL, rel_drho,
                potential_exchange_params.local_jacobian_relative_tolerance);
            return;
        }

        INFO(
            "{} local Jacobian diagnostic: pL_ip={} Pa, n_l_prev={}, n_l={}, h={}, analytic dn_l/dpL={}, FD dn_l/dpL={}, rel_err_nl={}, analytic drho_L_hat/dpL={}, FD drho_L_hat/dpL={}, rel_err_exchange={}, tolerance={}.",
            level_prefix, p_L_ip, n_l_prev, n_l, fd_data.perturbation,
            analytic_dn_l_dpL, fd_data.fd_dn_l_dpL, rel_dn_l,
            analytic_drho_L_hat_dpL, fd_data.fd_drho_L_hat_dpL, rel_drho,
            potential_exchange_params.local_jacobian_relative_tolerance);
    });
}

template <int DisplacementDim>
inline void updateMicroscaleHydraulicState(
    StatefulData<DisplacementDim>& state_current,
    StatefulDataPrev<DisplacementDim> const& state_previous, double const p_cap_ip,
    double const rho_LR, double const mu, double const dt,
    MPL::VariableArray& variables, MPL::VariableArray& variables_prev,
    PotentialExchangeLocalSolveContext const& local_context,
    std::optional<MicroPorosityParameters> const& micro_porosity_parameters,
    PotentialExchangeParameters const* const potential_exchange_parameters)
{
    auto& n_l = std::get<MicroWaterContent>(state_current);
    auto const n_l_prev = std::get<PrevState<MicroWaterContent>>(state_previous);
    auto& rho_lR = std::get<MicroLiquidDensity>(state_current);
    auto const rho_lR_prev = std::get<PrevState<MicroLiquidDensity>>(state_previous);

    double const n_l_prev_value = std::max(1e-16, **n_l_prev);
    *n_l = n_l_prev_value;

    if (!isPotentialExchangeEnabled(potential_exchange_parameters) ||
        !micro_porosity_parameters)
    {
        return;
    }

    auto const& potential_exchange_params = *potential_exchange_parameters;
    if (potential_exchange_params.local_nonlinear_solve_mode ==
        LocalNonlinearSolveMode::ScalarReferenceMassStorage)
    {
        auto const macro_potential = computeYoungLaplaceMacroPotential(
            -p_cap_ip, rho_LR, potential_exchange_params.pressure_tolerance);
        double const rho_lR_prev_value = std::max(1e-16, **rho_lR_prev);
        double const rho_l_prev = n_l_prev_value * rho_lR_prev_value;
        auto const coupled_update = solveReferenceMassStorageCoupledState(
            n_l_prev_value, rho_l_prev, rho_lR_prev_value, dt,
            rho_LR, micro_porosity_parameters->mass_exchange_coefficient, mu,
            macro_potential, local_context, potential_exchange_params);
        applyReferenceMassStorageLocalState<DisplacementDim>(
            state_current, state_previous, variables, variables_prev, rho_LR, local_context, potential_exchange_params,
            coupled_update);
        return;
    }

    auto const macro_potential = computeYoungLaplaceMacroPotential(
        -p_cap_ip, rho_LR, potential_exchange_params.pressure_tolerance);
    auto const n_l_update = solveImplicitMicroWaterContent(
        n_l_prev_value, dt, rho_LR,
        micro_porosity_parameters->mass_exchange_coefficient, mu,
        macro_potential, local_context, potential_exchange_params);

    *n_l = n_l_update.n_l;
    // Keep dsm_micromacro-mode rho_lR evolution consistent with the dsm_micromacro bridge:
    // rho_lR is updated from the active reduced micro EOS.
    *rho_lR = computeActiveMicroLiquidDensity(n_l_update.n_l, rho_LR,
                                                local_context, potential_exchange_params)
                  .rho_lR;

    auto& p_L_m = std::get<MicroPressure>(state_current);
    auto& S_L_m = std::get<MicroSaturation>(state_current);
    auto& rho_l_hat = std::get<MicroExchangeSource>(state_current);
    auto const compatibility_output =
        computeCompatibilityMicroHydraulicOutput(
            n_l_update.n_l, rho_LR, local_context, potential_exchange_params);
    *p_L_m = compatibility_output.p_L_m;
    *S_L_m = compatibility_output.S_L_m;
    rho_l_hat = MicroExchangeSource{n_l_update.exchange.rho_l_hat};
}

template <int DisplacementDim>
inline void updatePorositySplitState(
    StatefulData<DisplacementDim>& state_current,
    StatefulDataPrev<DisplacementDim> const& state_previous, double const phi,
    MPL::VariableArray& variables, MPL::VariableArray& variables_prev,
    PotentialExchangeParameters const* const potential_exchange_parameters)
{
    if (!isPotentialExchangeEnabled(potential_exchange_parameters))
    {
        return;
    }

    auto const mode =
        potential_exchange_parameters->local_nonlinear_solve_mode;
    auto& micro_porosity = std::get<MicroPorosity>(state_current);
    auto& transport_porosity =
        std::get<ProcessLib::ThermoRichardsMechanics::TransportPorosityData>(state_current)
            .phi;
    auto const phi_M_prev = std::get<PrevState<
        ProcessLib::ThermoRichardsMechanics::TransportPorosityData>>(state_previous)
                                ->phi;
    auto const n_l = std::max(1e-16, *std::get<MicroWaterContent>(state_current));

    if (mode == LocalNonlinearSolveMode::ScalarReferenceStorage ||
        mode == LocalNonlinearSolveMode::ScalarReferenceMassStorage)
    {
        // Keep dsm_micromacro support split aligned with the bridge law:
        // phi_m := n_l while transport_porosity remains the process state.
        // The dsm_micromacro support split (phi_m, phi_M) can step outside the
        // algebraic porosity bounds and should not collapse transport porosity.
        *micro_porosity = n_l;
        transport_porosity = phi_M_prev;
        variables.transport_porosity = transport_porosity;
        variables_prev.transport_porosity = phi_M_prev;
        return;
    }

    auto const phi_m_prev = **std::get<PrevState<MicroPorosity>>(state_previous);

    auto const transport_porosity_update =
        computeTransportPorosityUpdate(
            phi, phi_M_prev, phi_m_prev, n_l, variables.volumetric_strain,
            variables_prev.volumetric_strain,
            potential_exchange_parameters->macro_porosity_update_mode);

    *micro_porosity = transport_porosity_update.phi_m;
    transport_porosity = transport_porosity_update.phi_M;
    variables.transport_porosity = transport_porosity_update.phi_M;
    variables_prev.transport_porosity = transport_porosity_update.phi_M_prev;
}

template <int DisplacementDim>
inline void updateTotalPorosityState(
    StatefulData<DisplacementDim>& state_current,
    StatefulDataPrev<DisplacementDim> const& state_previous,
    double& phi, MPL::VariableArray& variables,
    MPL::VariableArray& variables_prev,
    PotentialExchangeParameters const* const potential_exchange_parameters)
{
    if (!isPotentialExchangeEnabled(potential_exchange_parameters))
    {
        return;
    }

    if (potential_exchange_parameters->local_nonlinear_solve_mode ==
            LocalNonlinearSolveMode::ScalarReferenceStorage ||
        potential_exchange_parameters->local_nonlinear_solve_mode ==
            LocalNonlinearSolveMode::ScalarReferenceMassStorage)
    {
        // In scalar dsm_micromacro-storage mode, micro porosity is support-state only.
        // Keep the process porosity state on the medium-law carrier.
        return;
    }

    auto const phi_m = *std::get<MicroPorosity>(state_current);
    auto const phi_m_prev = **std::get<PrevState<MicroPorosity>>(state_previous);
    auto const phi_M =
        std::get<ProcessLib::ThermoRichardsMechanics::TransportPorosityData>(state_current)
            .phi;
    auto const phi_M_prev =
        std::get<PrevState<
            ProcessLib::ThermoRichardsMechanics::TransportPorosityData>>(
            state_previous)
            ->phi;

    auto& porosity =
        std::get<ProcessLib::ThermoRichardsMechanics::PorosityData>(state_current).phi;
    phi = phi_M + phi_m;
    porosity = phi;
    variables.porosity = phi;
    variables_prev.porosity = phi_M_prev + phi_m_prev;
}

template <int DisplacementDim>
inline MathLib::KelvinVector::KelvinVectorType<DisplacementDim>
computeVdWRelaxationStressIncrement(
    double const p_L_m_prev, double const p_L_m,
    PotentialExchangeParameters const& potential_exchange_params)
{
    using KV = MathLib::KelvinVector::KelvinVectorType<DisplacementDim>;

    KV delta_sigma_vdw = KV::Zero();
    if (potential_exchange_params.micro_potential_convention !=
            MicroPotentialConvention::NegativeAttractive ||
        !(potential_exchange_params.vdw_relaxation_stress_gain > 0.0))
    {
        return delta_sigma_vdw;
    }

    double const delta_p_relaxation = std::max(0.0, p_L_m_prev - p_L_m);
    if (!(delta_p_relaxation > 0.0))
    {
        return delta_sigma_vdw;
    }

    auto const& identity2 = MathLib::KelvinVector::Invariants<
        MathLib::KelvinVector::kelvin_vector_dimensions(
            DisplacementDim)>::identity2;
    delta_sigma_vdw.noalias() -=
        potential_exchange_params.vdw_relaxation_stress_gain * delta_p_relaxation * identity2;
    return delta_sigma_vdw;
}

template <int DisplacementDim>
inline MathLib::KelvinVector::KelvinVectorType<DisplacementDim>
computeMicroWaterContentStressIncrement(
    double const n_l_prev, double const n_l,
    PotentialExchangeParameters const& potential_exchange_params)
{
    using KV = MathLib::KelvinVector::KelvinVectorType<DisplacementDim>;

    KV delta_sigma_nl = KV::Zero();
    if (!(potential_exchange_params.micro_water_content_stress_gain > 0.0))
    {
        return delta_sigma_nl;
    }

    double const delta_n_l = std::max(0.0, n_l - n_l_prev);
    if (!(delta_n_l > 0.0))
    {
        return delta_sigma_nl;
    }

    auto const& identity2 = MathLib::KelvinVector::Invariants<
        MathLib::KelvinVector::kelvin_vector_dimensions(
            DisplacementDim)>::identity2;
    delta_sigma_nl.noalias() -=
        potential_exchange_params.micro_water_content_stress_gain * delta_n_l * identity2;
    return delta_sigma_nl;
}

template <int DisplacementDim>
inline MathLib::KelvinVector::KelvinVectorType<DisplacementDim>
computeReferenceMicroPorositySwellingStressIncrement(
    double const phi_m_prev, double const phi_m,
    MathLib::KelvinVector::KelvinMatrixType<DisplacementDim> const& C_el,
    PotentialExchangeParameters const& potential_exchange_params)
{
    using KV = MathLib::KelvinVector::KelvinVectorType<DisplacementDim>;

    KV delta_sigma_sw = KV::Zero();
    if (!(potential_exchange_params.micro_water_content_swelling_slope > 0.0))
    {
        return delta_sigma_sw;
    }

    double const delta_phi_m = phi_m - phi_m_prev;
    if (!(std::isfinite(delta_phi_m) &&
          std::abs(delta_phi_m) > std::numeric_limits<double>::epsilon()))
    {
        return delta_sigma_sw;
    }

    auto const& identity2 = MathLib::KelvinVector::Invariants<
        MathLib::KelvinVector::kelvin_vector_dimensions(
            DisplacementDim)>::identity2;
    double const delta_eps_sw =
        potential_exchange_params.micro_water_content_swelling_slope * delta_phi_m;

    delta_sigma_sw.noalias() -= C_el * ((delta_eps_sw / 3.0) * identity2);
    return delta_sigma_sw;
}

template <int DisplacementDim>
inline MathLib::KelvinVector::KelvinVectorType<DisplacementDim>
computeSwellingStressIncrement(
    double const phi_m_prev, double const phi_m, double const n_l_prev,
    double const n_l, double const p_L_m_prev, double const p_L_m,
    MathLib::KelvinVector::KelvinMatrixType<DisplacementDim> const& C_el,
    PotentialExchangeParameters const& potential_exchange_params)
{
    using KV = MathLib::KelvinVector::KelvinVectorType<DisplacementDim>;

    KV delta_sigma_sw = KV::Zero();
    if (potential_exchange_params.micro_water_content_swelling_slope > 0.0)
    {
        delta_sigma_sw +=
            computeReferenceMicroPorositySwellingStressIncrement<
                DisplacementDim>(phi_m_prev, phi_m, C_el, potential_exchange_params);
    }

    if (potential_exchange_params.vdw_relaxation_stress_gain > 0.0)
    {
        delta_sigma_sw +=
            computeVdWRelaxationStressIncrement<DisplacementDim>(
                p_L_m_prev, p_L_m, potential_exchange_params);
    }

    if (potential_exchange_params.micro_water_content_stress_gain > 0.0)
    {
        delta_sigma_sw +=
            computeMicroWaterContentStressIncrement<DisplacementDim>(
                n_l_prev, n_l, potential_exchange_params);
    }

    return delta_sigma_sw;
}

template <int DisplacementDim>
inline void updateSwellingState(
    MaterialPropertyLib::Phase const& solid_phase,
    MathLib::KelvinVector::KelvinMatrixType<DisplacementDim> const& C_el,
    StatefulData<DisplacementDim>& state_current,
    StatefulDataPrev<DisplacementDim> const& state_previous,
    MPL::VariableArray& variables, MPL::VariableArray& variables_prev,
    ParameterLib::SpatialPosition const& x_position, double const t,
    double const dt,
    PotentialExchangeParameters const* const potential_exchange_parameters)
{
    if (!isPotentialExchangeEnabled(potential_exchange_parameters))
    {
        return;
    }

    auto const& potential_exchange_params = *potential_exchange_parameters;
    (void)solid_phase;
    (void)x_position;
    (void)t;
    (void)dt;

    auto const p_L_m_prev = **std::get<PrevState<MicroPressure>>(state_previous);
    auto const p_L_m = *std::get<MicroPressure>(state_current);
    auto const n_l_prev = **std::get<PrevState<MicroWaterContent>>(state_previous);
    auto const n_l = *std::get<MicroWaterContent>(state_current);
    auto const phi_m_prev = **std::get<PrevState<MicroPorosity>>(state_previous);
    auto const phi_m = *std::get<MicroPorosity>(state_current);

    auto& sigma_sw =
        std::get<ProcessLib::ThermoRichardsMechanics::
                     ConstitutiveStress_StrainTemperature::
                         SwellingDataStateful<DisplacementDim>>(state_current);
    auto const& sigma_sw_prev = std::get<
        PrevState<ProcessLib::ThermoRichardsMechanics::
                      ConstitutiveStress_StrainTemperature::
                          SwellingDataStateful<DisplacementDim>>>(state_previous);

    sigma_sw = *sigma_sw_prev;
    sigma_sw.sigma_sw +=
        computeSwellingStressIncrement<DisplacementDim>(
            phi_m_prev, phi_m, n_l_prev, n_l, p_L_m_prev, p_L_m, C_el, potential_exchange_params);

    auto const& identity2 = MathLib::KelvinVector::Invariants<
        MathLib::KelvinVector::kelvin_vector_dimensions(
            DisplacementDim)>::identity2;
    auto const C_el_inverse = C_el.inverse().eval();

    variables.volumetric_mechanical_strain =
        variables.volumetric_strain +
        identity2.transpose() * C_el_inverse * sigma_sw.sigma_sw;
    variables_prev.volumetric_mechanical_strain =
        variables_prev.volumetric_strain +
        identity2.transpose() * C_el_inverse * sigma_sw_prev->sigma_sw;
}

template <int DisplacementDim>
void updateSwellingStressAndVolumetricStrain(
    MaterialPropertyLib::Medium const& medium,
    MaterialPropertyLib::Phase const& solid_phase,
    MathLib::KelvinVector::KelvinMatrixType<DisplacementDim> const& C_el,
    double const rho_LR, double const mu,
    std::optional<MicroPorosityParameters> micro_porosity_parameters,
    PotentialExchangeParameters const* const potential_exchange_parameters,
    double const alpha, double const phi, double const p_cap_ip,
    MPL::VariableArray& variables, MPL::VariableArray& variables_prev,
    ParameterLib::SpatialPosition const& x_position, double const t,
    double const dt,
    ProcessLib::ThermoRichardsMechanics::ConstitutiveStress_StrainTemperature::
        SwellingDataStateful<DisplacementDim>& sigma_sw,
    PrevState<ProcessLib::ThermoRichardsMechanics::
                  ConstitutiveStress_StrainTemperature::SwellingDataStateful<
                      DisplacementDim>> const& sigma_sw_prev,
    PrevState<ProcessLib::ThermoRichardsMechanics::TransportPorosityData> const
        phi_M_prev,
    PrevState<ProcessLib::ThermoRichardsMechanics::PorosityData> const phi_prev,
    ProcessLib::ThermoRichardsMechanics::TransportPorosityData& phi_M,
    PrevState<MicroPressure> const p_L_m_prev,
    PrevState<MicroSaturation> const S_L_m_prev, MicroPressure& p_L_m,
    MicroSaturation& S_L_m)
{
    auto const& identity2 = MathLib::KelvinVector::Invariants<
        MathLib::KelvinVector::kelvin_vector_dimensions(
            DisplacementDim)>::identity2;
    bool const potential_exchange_enabled =
        isPotentialExchangeEnabled(potential_exchange_parameters);

    if (!medium.hasProperty(MPL::PropertyType::saturation_micro))
    {
        if (potential_exchange_enabled)
        {
            sigma_sw = *sigma_sw_prev;
            variables.volumetric_mechanical_strain =
                variables.volumetric_strain +
                identity2.transpose() * C_el.inverse() * sigma_sw.sigma_sw;
            variables_prev.volumetric_mechanical_strain =
                variables_prev.volumetric_strain + identity2.transpose() *
                                                       C_el.inverse() *
                                                       sigma_sw_prev->sigma_sw;
            return;
        }

        // If there is swelling, compute it. Update volumetric strain rate,
        // s.t. it corresponds to the mechanical part only.
        sigma_sw = *sigma_sw_prev;
        if (solid_phase.hasProperty(MPL::PropertyType::swelling_stress_rate))
        {
            auto const sigma_sw_dot =
                MathLib::KelvinVector::tensorToKelvin<DisplacementDim>(
                    MPL::formEigenTensor<3>(
                        solid_phase[MPL::PropertyType::swelling_stress_rate]
                            .value(variables, variables_prev, x_position, t,
                                   dt)));
            sigma_sw.sigma_sw += sigma_sw_dot * dt;

            variables.volumetric_mechanical_strain =
                variables.volumetric_strain +
                identity2.transpose() * C_el.inverse() * sigma_sw.sigma_sw;
            variables_prev.volumetric_mechanical_strain =
                variables_prev.volumetric_strain + identity2.transpose() *
                                                       C_el.inverse() *
                                                       sigma_sw_prev->sigma_sw;
        }
        else
        {
            variables.volumetric_mechanical_strain =
                variables.volumetric_strain;
            variables_prev.volumetric_mechanical_strain =
                variables_prev.volumetric_strain;
        }
    }

    // TODO (naumov) saturation_micro must be always defined together with
    // the micro_porosity_parameters.
    if (medium.hasProperty(MPL::PropertyType::saturation_micro))
    {
        if (potential_exchange_enabled)
        {
            phi_M.phi = phi_M_prev->phi;
            variables_prev.transport_porosity = phi_M_prev->phi;
            variables.transport_porosity = phi_M.phi;

            *p_L_m = **p_L_m_prev;
            *S_L_m = **S_L_m_prev;
            sigma_sw = *sigma_sw_prev;

            variables.volumetric_mechanical_strain =
                variables.volumetric_strain +
                identity2.transpose() * C_el.inverse() * sigma_sw.sigma_sw;
            variables_prev.volumetric_mechanical_strain =
                variables_prev.volumetric_strain +
                identity2.transpose() * C_el.inverse() *
                    sigma_sw_prev->sigma_sw;
            return;
        }

        double const phi_m_prev = phi_prev->phi - phi_M_prev->phi;

        auto const [delta_phi_m, delta_e_sw, delta_p_L_m, delta_sigma_sw] =
            computeMicroPorosity<DisplacementDim>(
                identity2.transpose() * C_el.inverse(), rho_LR, mu,
                *micro_porosity_parameters, alpha, phi, -p_cap_ip, **p_L_m_prev,
                variables_prev, **S_L_m_prev, phi_m_prev, x_position, t, dt,
                medium.property(MPL::PropertyType::saturation_micro),
                solid_phase.property(MPL::PropertyType::swelling_stress_rate));

        phi_M.phi = phi - (phi_m_prev + delta_phi_m);
        variables_prev.transport_porosity = phi_M_prev->phi;
        variables.transport_porosity = phi_M.phi;

        *p_L_m = **p_L_m_prev + delta_p_L_m;
        {  // Update micro saturation.
            MPL::VariableArray variables_prev;
            variables_prev.capillary_pressure = -**p_L_m_prev;
            MPL::VariableArray variables;
            variables.capillary_pressure = -*p_L_m;

            *S_L_m = medium.property(MPL::PropertyType::saturation_micro)
                         .template value<double>(variables, x_position, t, dt);
        }
        sigma_sw.sigma_sw = sigma_sw_prev->sigma_sw + delta_sigma_sw;
    }
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
RichardsMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                ShapeFunctionPressure, DisplacementDim>::
    RichardsMechanicsLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const /*local_matrix_size*/,
        NumLib::GenericIntegrationMethod const& integration_method,
        bool const is_axially_symmetric,
        RichardsMechanicsProcessData<DisplacementDim>& process_data)
    : LocalAssemblerInterface<DisplacementDim>{
          e, integration_method, is_axially_symmetric, process_data}
{
    unsigned const n_integration_points =
        this->integration_method_.getNumberOfPoints();

    ip_data_.resize(n_integration_points);
    secondary_data_.N_u.resize(n_integration_points);

    auto const shape_matrices_u =
        NumLib::initShapeMatrices<ShapeFunctionDisplacement,
                                  ShapeMatricesTypeDisplacement,
                                  DisplacementDim>(e, is_axially_symmetric,
                                                   this->integration_method_);

    auto const shape_matrices_p =
        NumLib::initShapeMatrices<ShapeFunctionPressure,
                                  ShapeMatricesTypePressure, DisplacementDim>(
            e, is_axially_symmetric, this->integration_method_);

    auto const& medium =
        this->process_data_.media_map.getMedium(this->element_.getID());

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto& ip_data = ip_data_[ip];
        auto const& sm_u = shape_matrices_u[ip];
        ip_data_[ip].integration_weight =
            this->integration_method_.getWeightedPoint(ip).getWeight() *
            sm_u.integralMeasure * sm_u.detJ;

        ip_data.N_u = sm_u.N;
        ip_data.dNdx_u = sm_u.dNdx;

        ParameterLib::SpatialPosition x_position = {
            std::nullopt, this->element_.getID(),
            MathLib::Point3d(
                NumLib::interpolateCoordinates<ShapeFunctionDisplacement,
                                               ShapeMatricesTypeDisplacement>(
                    this->element_, ip_data.N_u))};

        ip_data.N_p = shape_matrices_p[ip].N;
        ip_data.dNdx_p = shape_matrices_p[ip].dNdx;

        // Initial porosity. Could be read from integration point data or mesh.
        auto& porosity =
            std::get<ProcessLib::ThermoRichardsMechanics::PorosityData>(
                this->current_states_[ip])
                .phi;
        porosity = medium->property(MPL::porosity)
                       .template initialValue<double>(
                           x_position,
                           std::numeric_limits<
                               double>::quiet_NaN() /* t independent */);

        auto& transport_porosity =
            std::get<
                ProcessLib::ThermoRichardsMechanics::TransportPorosityData>(
                this->current_states_[ip])
                .phi;
        transport_porosity = porosity;
        if (medium->hasProperty(MPL::PropertyType::transport_porosity))
        {
            transport_porosity =
                medium->property(MPL::transport_porosity)
                    .template initialValue<double>(
                        x_position,
                        std::numeric_limits<
                            double>::quiet_NaN() /* t independent */);
        }

        secondary_data_.N_u[ip] = shape_matrices_u[ip].N;
    }
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
void RichardsMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                     ShapeFunctionPressure, DisplacementDim>::
    setInitialConditionsConcrete(Eigen::VectorXd const local_x,
                                 double const t,
                                 int const /*process_id*/)
{
    assert(local_x.size() == pressure_size + displacement_size);

    auto const [p_L, u] = localDOF(local_x);

    constexpr double dt = std::numeric_limits<double>::quiet_NaN();
    auto const& medium =
        this->process_data_.media_map.getMedium(this->element_.getID());
    MPL::VariableArray variables;

    auto const& solid_phase = medium->phase("Solid");

    auto const& identity2 = MathLib::KelvinVector::Invariants<
        MathLib::KelvinVector::kelvin_vector_dimensions(
            DisplacementDim)>::identity2;

    unsigned const n_integration_points =
        this->integration_method_.getNumberOfPoints();
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto const& N_p = ip_data_[ip].N_p;

        ParameterLib::SpatialPosition x_position = {
            std::nullopt, this->element_.getID(),
            MathLib::Point3d(
                NumLib::interpolateCoordinates<ShapeFunctionPressure,
                                               ShapeMatricesTypePressure>(
                    this->element_, N_p))};

        double p_cap_ip;
        NumLib::shapeFunctionInterpolate(-p_L, N_p, p_cap_ip);

        variables.capillary_pressure = p_cap_ip;
        variables.liquid_phase_pressure = -p_cap_ip;
        // setting pG to 1 atm
        // TODO : rewrite equations s.t. p_L = pG-p_cap
        variables.gas_phase_pressure = 1.0e5;

        {
            auto& p_L_m = std::get<MicroPressure>(this->current_states_[ip]);
            auto& p_L_m_prev =
                std::get<PrevState<MicroPressure>>(this->prev_states_[ip]);
            **p_L_m_prev = -p_cap_ip;
            *p_L_m = -p_cap_ip;
        }

        auto const temperature =
            medium->property(MPL::PropertyType::reference_temperature)
                .template value<double>(variables, x_position, t, dt);
        variables.temperature = temperature;

        auto& S_L_prev =
            std::get<
                PrevState<ProcessLib::ThermoRichardsMechanics::SaturationData>>(
                this->prev_states_[ip])
                ->S_L;
        S_L_prev = medium->property(MPL::PropertyType::saturation)
                       .template value<double>(variables, x_position, t, dt);

        if (this->process_data_.initial_stress.isTotalStress())
        {
            auto const alpha_b =
                medium->property(MPL::PropertyType::biot_coefficient)
                    .template value<double>(variables, x_position, t, dt);

            variables.liquid_saturation = S_L_prev;
            double const chi_S_L =
                medium->property(MPL::PropertyType::bishops_effective_stress)
                    .template value<double>(variables, x_position, t, dt);

            // Initial stresses are total stress, which were assigned to
            // sigma_eff in
            // RichardsMechanicsLocalAssembler::initializeConcrete().
            auto& sigma_eff =
                std::get<ProcessLib::ConstitutiveRelations::EffectiveStressData<
                    DisplacementDim>>(this->current_states_[ip]);

            auto& sigma_eff_prev =
                std::get<PrevState<ProcessLib::ConstitutiveRelations::
                                       EffectiveStressData<DisplacementDim>>>(
                    this->prev_states_[ip]);

            // Reset sigma_eff to effective stress
            sigma_eff.sigma_eff.noalias() +=
                chi_S_L * alpha_b * (-p_cap_ip) * identity2;
            sigma_eff_prev->sigma_eff = sigma_eff.sigma_eff;
        }

        if (medium->hasProperty(MPL::PropertyType::saturation_micro))
        {
            MPL::VariableArray vars;
            vars.capillary_pressure = p_cap_ip;

            auto& S_L_m = std::get<MicroSaturation>(this->current_states_[ip]);
            auto& S_L_m_prev =
                std::get<PrevState<MicroSaturation>>(this->prev_states_[ip]);

            *S_L_m = medium->property(MPL::PropertyType::saturation_micro)
                         .template value<double>(vars, x_position, t, dt);
            *S_L_m_prev = S_L_m;
        }

        {
            auto& n_l = std::get<MicroWaterContent>(this->current_states_[ip]);
            auto& n_l_prev =
        std::get<PrevState<MicroWaterContent>>(this->prev_states_[ip]);
    auto& rho_lR =
        std::get<MicroLiquidDensity>(this->current_states_[ip]);
    auto& rho_lR_prev =
        std::get<PrevState<MicroLiquidDensity>>(this->prev_states_[ip]);
    auto& phi_m = std::get<MicroPorosity>(this->current_states_[ip]);
    auto& phi_m_prev =
        std::get<PrevState<MicroPorosity>>(this->prev_states_[ip]);

    // Default fallback keeps state positive for vdW algebra.
    double n_l_initial = 1e-6;
    double rho_lR_initial = 1.0;
            if (medium->hasProperty(MPL::PropertyType::saturation_micro))
            {
                auto const S_L_m_init =
                    *std::get<MicroSaturation>(this->current_states_[ip]);
                n_l_initial = std::max(1e-12, S_L_m_init);
            }

            auto const* const potential_exchange_params_ptr =
                this->getPotentialExchangeParameters();

            if (isPotentialExchangeEnabled(potential_exchange_params_ptr))
            {
                auto const porosity =
                    std::get<ProcessLib::ThermoRichardsMechanics::PorosityData>(
                        this->current_states_[ip])
                        .phi;
                auto const transport_porosity =
                    std::get<ProcessLib::ThermoRichardsMechanics::
                                 TransportPorosityData>(this->current_states_[ip])
                        .phi;
                n_l_initial = std::max(1e-12, porosity - transport_porosity);

                n_l_initial =
                    potential_exchange_params_ptr->initial_micro_water_content
                        .value_or(n_l_initial);
                rho_lR_initial = std::max(
                    1e-16,
                    potential_exchange_params_ptr
                        ->micro_liquid_density_reference);
            }

            *n_l = n_l_initial;
            **n_l_prev = n_l_initial;
            *rho_lR = rho_lR_initial;
            **rho_lR_prev = rho_lR_initial;
            *phi_m = n_l_initial;
            **phi_m_prev = n_l_initial;

            if (isPotentialExchangeEnabled(potential_exchange_params_ptr))
            {
                auto const porosity =
                    std::get<ProcessLib::ThermoRichardsMechanics::PorosityData>(
                        this->current_states_[ip])
                        .phi;
                auto const transport_porosity_init =
                    std::get<ProcessLib::ThermoRichardsMechanics::
                                 TransportPorosityData>(this->current_states_[ip])
                        .phi;
                auto const rho_LR_initial =
                    medium->phase("AqueousLiquid")
                        .property(MPL::PropertyType::density)
                        .template value<double>(variables, x_position, t, dt);
                auto const compatibility_output =
                    computeCompatibilityMicroHydraulicOutput(
                        n_l_initial, rho_LR_initial,
                        {.phi = porosity,
                         .phi_M_prev = transport_porosity_init,
                         .phi_m_prev = n_l_initial,
                         .volumetric_strain = 0.0,
                         .volumetric_strain_prev = 0.0},
                        *potential_exchange_params_ptr);
                auto& p_L_m =
                    std::get<MicroPressure>(this->current_states_[ip]);
                auto& p_L_m_prev =
                    std::get<PrevState<MicroPressure>>(this->prev_states_[ip]);
                auto& S_L_m =
                    std::get<MicroSaturation>(this->current_states_[ip]);
                auto& S_L_m_prev =
                    std::get<PrevState<MicroSaturation>>(
                        this->prev_states_[ip]);
                *p_L_m = compatibility_output.p_L_m;
                **p_L_m_prev = compatibility_output.p_L_m;
                *S_L_m = compatibility_output.S_L_m;
                **S_L_m_prev = compatibility_output.S_L_m;

                auto& transport_porosity =
                    std::get<ProcessLib::ThermoRichardsMechanics::
                                 TransportPorosityData>(
                        this->current_states_[ip])
                        .phi;
                auto& transport_porosity_prev = std::get<PrevState<
                    ProcessLib::ThermoRichardsMechanics::TransportPorosityData>>(
                    this->prev_states_[ip]);
                auto const transport_porosity_update =
                    computeTransportPorosityUpdate(
                        porosity, transport_porosity, n_l_initial, n_l_initial,
                        /*volumetric_strain=*/0.0,
                        /*volumetric_strain_prev=*/0.0,
                        potential_exchange_params_ptr
                            ->macro_porosity_update_mode);
                *phi_m = transport_porosity_update.phi_m;
                **phi_m_prev = transport_porosity_update.phi_m_prev;
                transport_porosity = transport_porosity_update.phi_M;
                transport_porosity_prev->phi =
                    transport_porosity_update.phi_M_prev;
            }
        }

        // Set eps_m_prev from potentially non-zero eps and sigma_sw from
        // restart.
        auto& state_current = this->current_states_[ip];
        variables.stress =
            std::get<ProcessLib::ConstitutiveRelations::EffectiveStressData<
                DisplacementDim>>(state_current)
                .sigma_eff;

        auto const& N_u = ip_data_[ip].N_u;
        auto const& dNdx_u = ip_data_[ip].dNdx_u;
        auto const x_coord =
            x_position.getCoordinates().value()[0];  // r for axisymetric
        auto const B =
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunctionDisplacement::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                dNdx_u, N_u, x_coord, this->is_axially_symmetric_);
        auto& eps =
            std::get<StrainData<DisplacementDim>>(this->current_states_[ip])
                .eps;
        eps.noalias() = B * u;

        // Set mechanical strain temporary to compute tangent stiffness.
        variables.mechanical_strain
            .emplace<MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
                eps);

        auto const C_el = ip_data_[ip].computeElasticTangentStiffness(
            variables, t, x_position, dt, this->solid_material_,
            *this->material_states_[ip].material_state_variables);

        auto const& sigma_sw =
            std::get<ProcessLib::ThermoRichardsMechanics::
                         ConstitutiveStress_StrainTemperature::
                             SwellingDataStateful<DisplacementDim>>(
                this->current_states_[ip])
                .sigma_sw;
        auto& eps_m_prev =
            std::get<PrevState<ProcessLib::ConstitutiveRelations::
                                   MechanicalStrainData<DisplacementDim>>>(
                this->prev_states_[ip])
                ->eps_m;

        bool const swelling_stress_active =
            solid_phase.hasProperty(MPL::PropertyType::swelling_stress_rate) ||
            isPotentialExchangeEnabled(this->getPotentialExchangeParameters());
        eps_m_prev.noalias() =
            swelling_stress_active ? eps + C_el.inverse() * sigma_sw : eps;
    }
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
void RichardsMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure,
    DisplacementDim>::assemble(double const t, double const dt,
                               std::vector<double> const& local_x,
                               std::vector<double> const& local_x_prev,
                               std::vector<double>& local_M_data,
                               std::vector<double>& local_K_data,
                               std::vector<double>& local_rhs_data)
{
    assert(local_x.size() == pressure_size + displacement_size);

    auto const [p_L, u] = localDOF(local_x);
    auto const [p_L_prev, u_prev] = localDOF(local_x_prev);

    auto K = MathLib::createZeroedMatrix<
        typename ShapeMatricesTypeDisplacement::template MatrixType<
            displacement_size + pressure_size,
            displacement_size + pressure_size>>(
        local_K_data, displacement_size + pressure_size,
        displacement_size + pressure_size);

    auto M = MathLib::createZeroedMatrix<
        typename ShapeMatricesTypeDisplacement::template MatrixType<
            displacement_size + pressure_size,
            displacement_size + pressure_size>>(
        local_M_data, displacement_size + pressure_size,
        displacement_size + pressure_size);

    auto rhs = MathLib::createZeroedVector<
        typename ShapeMatricesTypeDisplacement::template VectorType<
            displacement_size + pressure_size>>(
        local_rhs_data, displacement_size + pressure_size);

    auto const& identity2 = MathLib::KelvinVector::Invariants<
        MathLib::KelvinVector::kelvin_vector_dimensions(
            DisplacementDim)>::identity2;

    auto const& medium =
        this->process_data_.media_map.getMedium(this->element_.getID());
    auto const& liquid_phase = medium->phase("AqueousLiquid");
    auto const& solid_phase = medium->phase("Solid");
    MPL::VariableArray variables;
    MPL::VariableArray variables_prev;

    ParameterLib::SpatialPosition x_position;
    x_position.setElementID(this->element_.getID());

    unsigned const n_integration_points =
        this->integration_method_.getNumberOfPoints();
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto const& w = ip_data_[ip].integration_weight;

        auto const& N_u = ip_data_[ip].N_u;
        auto const& dNdx_u = ip_data_[ip].dNdx_u;

        auto const& N_p = ip_data_[ip].N_p;
        auto const& dNdx_p = ip_data_[ip].dNdx_p;

        x_position = {
            std::nullopt, this->element_.getID(),
            MathLib::Point3d(
                NumLib::interpolateCoordinates<ShapeFunctionDisplacement,
                                               ShapeMatricesTypeDisplacement>(
                    this->element_, N_u))};
        auto const x_coord = x_position.getCoordinates().value()[0];

        auto const B =
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunctionDisplacement::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                dNdx_u, N_u, x_coord, this->is_axially_symmetric_);

        auto& eps =
            std::get<StrainData<DisplacementDim>>(this->current_states_[ip]);
        eps.eps.noalias() = B * u;

        auto& S_L =
            std::get<ProcessLib::ThermoRichardsMechanics::SaturationData>(
                this->current_states_[ip])
                .S_L;
        auto const S_L_prev =
            std::get<
                PrevState<ProcessLib::ThermoRichardsMechanics::SaturationData>>(
                this->prev_states_[ip])
                ->S_L;

        double p_cap_ip;
        NumLib::shapeFunctionInterpolate(-p_L, N_p, p_cap_ip);

        double p_cap_prev_ip;
        NumLib::shapeFunctionInterpolate(-p_L_prev, N_p, p_cap_prev_ip);

        variables.capillary_pressure = p_cap_ip;
        variables.liquid_phase_pressure = -p_cap_ip;
        // setting pG to 1 atm
        // TODO : rewrite equations s.t. p_L = pG-p_cap
        variables.gas_phase_pressure = 1.0e5;

        auto const temperature =
            medium->property(MPL::PropertyType::reference_temperature)
                .template value<double>(variables, x_position, t, dt);
        variables.temperature = temperature;

        auto const alpha =
            medium->property(MPL::PropertyType::biot_coefficient)
                .template value<double>(variables, x_position, t, dt);
        auto& state_current = this->current_states_[ip];
        variables.stress =
            std::get<ProcessLib::ConstitutiveRelations::EffectiveStressData<
                DisplacementDim>>(state_current)
                .sigma_eff;
        // Set mechanical strain temporary to compute tangent stiffness.
        variables.mechanical_strain
            .emplace<MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
                eps.eps);
        auto const C_el = ip_data_[ip].computeElasticTangentStiffness(
            variables, t, x_position, dt, this->solid_material_,
            *this->material_states_[ip].material_state_variables);

        auto const beta_SR = (1 - alpha) / this->solid_material_.getBulkModulus(
                                               t, x_position, &C_el);
        variables.grain_compressibility = beta_SR;

        auto const rho_LR =
            liquid_phase.property(MPL::PropertyType::density)
                .template value<double>(variables, x_position, t, dt);
        variables.density = rho_LR;
        auto const& b = this->process_data_.specific_body_force;

        S_L = medium->property(MPL::PropertyType::saturation)
                  .template value<double>(variables, x_position, t, dt);
        variables.liquid_saturation = S_L;
        variables_prev.liquid_saturation = S_L_prev;

        // tangent derivative for Jacobian
        double const dS_L_dp_cap =
            medium->property(MPL::PropertyType::saturation)
                .template dValue<double>(variables,
                                         MPL::Variable::capillary_pressure,
                                         x_position, t, dt);
        // secant derivative from time discretization for storage
        // use tangent, if secant is not available
        double const DeltaS_L_Deltap_cap =
            (p_cap_ip == p_cap_prev_ip)
                ? dS_L_dp_cap
                : (S_L - S_L_prev) / (p_cap_ip - p_cap_prev_ip);

        auto const chi = [medium, x_position, t, dt](double const S_L)
        {
            MPL::VariableArray vs;
            vs.liquid_saturation = S_L;
            return medium->property(MPL::PropertyType::bishops_effective_stress)
                .template value<double>(vs, x_position, t, dt);
        };
        double const chi_S_L = chi(S_L);
        double const chi_S_L_prev = chi(S_L_prev);

        double const p_FR = -chi_S_L * p_cap_ip;
        variables.effective_pore_pressure = p_FR;
        variables_prev.effective_pore_pressure = -chi_S_L_prev * p_cap_prev_ip;

        // Set volumetric strain rate for the general case without swelling.
        variables.volumetric_strain = Invariants::trace(eps.eps);
        variables_prev.volumetric_strain = Invariants::trace(B * u_prev);

        auto& phi = std::get<ProcessLib::ThermoRichardsMechanics::PorosityData>(
                        this->current_states_[ip])
                        .phi;
        {  // Porosity update
            auto const phi_prev = std::get<PrevState<
                ProcessLib::ThermoRichardsMechanics::PorosityData>>(
                                      this->prev_states_[ip])
                                      ->phi;
            variables_prev.porosity = phi_prev;
            phi = medium->property(MPL::PropertyType::porosity)
                      .template value<double>(variables, variables_prev,
                                              x_position, t, dt);
            variables.porosity = phi;
        }

        if (alpha < phi)
        {
            OGS_FATAL(
                "RichardsMechanics: Biot-coefficient {} is smaller than "
                "porosity {} in element/integration point {}/{}.",
                alpha, phi, this->element_.getID(), ip);
        }

        // Swelling and possibly volumetric strain rate update.
        {
            auto& sigma_sw =
                std::get<ProcessLib::ThermoRichardsMechanics::
                             ConstitutiveStress_StrainTemperature::
                                 SwellingDataStateful<DisplacementDim>>(
                    this->current_states_[ip])
                    .sigma_sw;
            auto const& sigma_sw_prev = std::get<PrevState<
                ProcessLib::ThermoRichardsMechanics::
                    ConstitutiveStress_StrainTemperature::SwellingDataStateful<
                        DisplacementDim>>>(this->prev_states_[ip])
                                            ->sigma_sw;

            // If there is swelling, compute it. Update volumetric strain rate,
            // s.t. it corresponds to the mechanical part only.
            sigma_sw = sigma_sw_prev;
            if (solid_phase.hasProperty(
                    MPL::PropertyType::swelling_stress_rate))
            {
                auto const sigma_sw_dot =
                    MathLib::KelvinVector::tensorToKelvin<DisplacementDim>(
                        MPL::formEigenTensor<3>(
                            solid_phase[MPL::PropertyType::swelling_stress_rate]
                                .value(variables, variables_prev, x_position, t,
                                       dt)));
                sigma_sw += sigma_sw_dot * dt;

                variables.volumetric_mechanical_strain =
                    variables.volumetric_strain +
                    identity2.transpose() * C_el.inverse() * sigma_sw;
                variables_prev.volumetric_mechanical_strain =
                    variables_prev.volumetric_strain +
                    identity2.transpose() * C_el.inverse() * sigma_sw_prev;
            }
            else
            {
                variables.volumetric_mechanical_strain =
                    variables.volumetric_strain;
                variables_prev.volumetric_mechanical_strain =
                    variables_prev.volumetric_strain;
            }

            if (medium->hasProperty(MPL::PropertyType::transport_porosity))
            {
                auto& transport_porosity =
                    std::get<ProcessLib::ThermoRichardsMechanics::
                                 TransportPorosityData>(
                        this->current_states_[ip])
                        .phi;
                auto const transport_porosity_prev =
                    std::get<PrevState<ProcessLib::ThermoRichardsMechanics::
                                           TransportPorosityData>>(
                        this->prev_states_[ip])
                        ->phi;
                variables_prev.transport_porosity = transport_porosity_prev;

                transport_porosity =
                    medium->property(MPL::PropertyType::transport_porosity)
                        .template value<double>(variables, variables_prev,
                                                x_position, t, dt);
                variables.transport_porosity = transport_porosity;
            }
            else
            {
                variables.transport_porosity = phi;
            }
        }

        double const k_rel =
            medium->property(MPL::PropertyType::relative_permeability)
                .template value<double>(variables, x_position, t, dt);
        auto const mu =
            liquid_phase.property(MPL::PropertyType::viscosity)
                .template value<double>(variables, x_position, t, dt);

        auto const& sigma_sw =
            std::get<ProcessLib::ThermoRichardsMechanics::
                         ConstitutiveStress_StrainTemperature::
                             SwellingDataStateful<DisplacementDim>>(
                this->current_states_[ip])
                .sigma_sw;
        auto const& sigma_eff =
            std::get<ProcessLib::ConstitutiveRelations::EffectiveStressData<
                DisplacementDim>>(this->current_states_[ip])
                .sigma_eff;

        // Set mechanical variables for the intrinsic permeability model
        // For stress dependent permeability.
        {
            auto const sigma_total =
                (sigma_eff - alpha * p_FR * identity2).eval();

            // For stress dependent permeability.
            variables.total_stress.emplace<SymmetricTensor>(
                MathLib::KelvinVector::kelvinVectorToSymmetricTensor(
                    sigma_total));
        }

        variables.equivalent_plastic_strain =
            this->material_states_[ip]
                .material_state_variables->getEquivalentPlasticStrain();

        auto const K_intrinsic = MPL::formEigenTensor<DisplacementDim>(
            medium->property(MPL::PropertyType::permeability)
                .value(variables, x_position, t, dt));

        GlobalDimMatrixType const rho_K_over_mu =
            K_intrinsic * rho_LR * k_rel / mu;

        //
        // displacement equation, displacement part
        //
        {
            auto& eps_m = std::get<ProcessLib::ConstitutiveRelations::
                                       MechanicalStrainData<DisplacementDim>>(
                              this->current_states_[ip])
                              .eps_m;
            bool const swelling_stress_active =
                solid_phase.hasProperty(MPL::PropertyType::swelling_stress_rate) ||
                isPotentialExchangeEnabled(
                    this->getPotentialExchangeParameters());
            eps_m.noalias() = swelling_stress_active
                                  ? eps.eps + C_el.inverse() * sigma_sw
                                  : eps.eps;
            variables.mechanical_strain.emplace<
                MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
                eps_m);
        }

        {
            auto& state_current = this->current_states_[ip];
            auto const& state_previous = this->prev_states_[ip];
            auto& sigma_eff =
                std::get<ProcessLib::ConstitutiveRelations::EffectiveStressData<
                    DisplacementDim>>(state_current);
            auto const& sigma_eff_prev =
                std::get<PrevState<ProcessLib::ConstitutiveRelations::
                                       EffectiveStressData<DisplacementDim>>>(
                    state_previous);
            auto const& eps_m =
                std::get<ProcessLib::ConstitutiveRelations::
                             MechanicalStrainData<DisplacementDim>>(state_current);
            auto& eps_m_prev =
                std::get<PrevState<ProcessLib::ConstitutiveRelations::
                                       MechanicalStrainData<DisplacementDim>>>(
                    state_previous);

            auto const C = ip_data_[ip].updateConstitutiveRelation(
                variables, t, x_position, dt, temperature, sigma_eff,
                sigma_eff_prev, eps_m, eps_m_prev, this->solid_material_,
                this->material_states_[ip].material_state_variables);

            if (this->process_data_.use_numerical_jacobian)
            {
                K.template block<displacement_size, displacement_size>(
                     displacement_index, displacement_index)
                    .noalias() += B.transpose() * C * B * w;
            }
        }

        // p_SR
        variables.solid_grain_pressure =
            p_FR - sigma_eff.dot(identity2) / (3 * (1 - phi));
        auto const rho_SR =
            solid_phase.property(MPL::PropertyType::density)
                .template value<double>(variables, x_position, t, dt);

        //
        // displacement equation, displacement part
        //
        double const rho = rho_SR * (1 - phi) + S_L * phi * rho_LR;
        rhs.template segment<displacement_size>(displacement_index).noalias() -=
            (B.transpose() * sigma_eff - N_u_op(N_u).transpose() * rho * b) * w;

        //
        // pressure equation, pressure part.
        //
        auto const beta_LR =
            1 / rho_LR *
            liquid_phase.property(MPL::PropertyType::density)
                .template dValue<double>(variables,
                                         MPL::Variable::liquid_phase_pressure,
                                         x_position, t, dt);

        double const a0 = S_L * (alpha - phi) * beta_SR;
        // Volumetric average specific storage of the solid and fluid phases.
        double const specific_storage =
            DeltaS_L_Deltap_cap * (p_cap_ip * a0 - phi) +
            S_L * (phi * beta_LR + a0);
        M.template block<pressure_size, pressure_size>(pressure_index,
                                                       pressure_index)
            .noalias() += N_p.transpose() * rho_LR * specific_storage * N_p * w;

        K.template block<pressure_size, pressure_size>(pressure_index,
                                                       pressure_index)
            .noalias() += dNdx_p.transpose() * rho_K_over_mu * dNdx_p * w;

        rhs.template segment<pressure_size>(pressure_index).noalias() +=
            dNdx_p.transpose() * rho_LR * rho_K_over_mu * b * w;

        auto const* const potential_exchange_params_ptr =
            this->getPotentialExchangeParameters();
        bool const potential_exchange_enabled =
            isPotentialExchangeEnabled(potential_exchange_params_ptr);
        if ((medium->hasProperty(MPL::PropertyType::saturation_micro) ||
             potential_exchange_enabled) &&
            this->process_data_.micro_porosity_parameters)
        {
            double const alpha_bar =
                this->process_data_.micro_porosity_parameters
                    ->mass_exchange_coefficient;
            auto const p_L_m =
                *std::get<MicroPressure>(this->current_states_[ip]);
            double const p_L_ip = -p_cap_ip;
            double const pressure_tolerance =
                getPotentialPressureTolerance(
                    potential_exchange_params_ptr);

            bool use_vdw_micro_potential_for_active_exchange = false;
            double mu_lR_vdw = 0.0;
            double dmu_lR_vdw_drho_lR = 0.0;

            if (potential_exchange_enabled)
            {
                auto const n_l =
                    std::max(1e-16,
                             *std::get<MicroWaterContent>(
                                 this->current_states_[ip]));
                auto const transport_porosity_prev =
                    std::get<PrevState<ProcessLib::ThermoRichardsMechanics::
                                           TransportPorosityData>>(
                        this->prev_states_[ip])
                        ->phi;
                auto const n_l_prev = **std::get<PrevState<MicroWaterContent>>(
                    this->prev_states_[ip]);
                PotentialExchangeLocalSolveContext const local_solve_context{
                    .phi = phi,
                    .phi_M_prev = transport_porosity_prev,
                    .phi_m_prev = n_l_prev,
                    .volumetric_strain = variables.volumetric_strain,
                    .volumetric_strain_prev = variables_prev.volumetric_strain};
                auto const micro_potential = computeActiveMicroPotential(
                    n_l, rho_LR, local_solve_context,
                    *potential_exchange_params_ptr);
                use_vdw_micro_potential_for_active_exchange = true;
                mu_lR_vdw = micro_potential.mu_lR;
                dmu_lR_vdw_drho_lR = micro_potential.dmu_lR_drho_lR;
            }

            auto const potential_exchange_result = computePotentialExchangeUpdate(
                alpha_bar, mu, p_L_ip, p_L_m, rho_LR, beta_LR,
                pressure_tolerance, potential_exchange_enabled,
                use_vdw_micro_potential_for_active_exchange, mu_lR_vdw,
                dmu_lR_vdw_drho_lR,
                /*use_custom_dmu_lR_vdw_dpL=*/false, /*dmu_lR_vdw_dpL=*/0.0,
                /*use_fd_jacobian_for_direct_macro_derivative=*/false,
                /*fd_jacobian_perturbation=*/1e-8);
            rhs.template segment<pressure_size>(pressure_index).noalias() +=
                N_p.transpose() * potential_exchange_result.exchange.rho_L_hat * w;
        }

        //
        // displacement equation, pressure part
        //
        K.template block<displacement_size, pressure_size>(displacement_index,
                                                           pressure_index)
            .noalias() -= B.transpose() * alpha * chi_S_L * identity2 * N_p * w;

        //
        // pressure equation, displacement part.
        //
        M.template block<pressure_size, displacement_size>(pressure_index,
                                                           displacement_index)
            .noalias() += N_p.transpose() * S_L * rho_LR * alpha *
                          identity2.transpose() * B * w;
    }

    if (this->process_data_.apply_mass_lumping)
    {
        auto pressure_mass_block_diag = M.template block<pressure_size, pressure_size>(
            pressure_index, pressure_index);
        pressure_mass_block_diag = pressure_mass_block_diag.colwise().sum().eval().asDiagonal();
    }
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
void RichardsMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                     ShapeFunctionPressure, DisplacementDim>::
    assembleWithJacobianEvalConstitutiveSetting(
        double const t, double const dt,
        ParameterLib::SpatialPosition const& x_position,
        RichardsMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                        ShapeFunctionPressure,
                                        DisplacementDim>::IpData& ip_data,
        MPL::VariableArray& variables, MPL::VariableArray& variables_prev,
        MPL::Medium const* const medium, TemperatureData const T_data,
        CapillaryPressureData<DisplacementDim> const& p_cap_data,
        ConstitutiveData<DisplacementDim>& constitutive_data,
        StatefulData<DisplacementDim>& state_current,
        StatefulDataPrev<DisplacementDim> const& state_previous,
        std::optional<MicroPorosityParameters> const& micro_porosity_parameters,
        PotentialExchangeParameters const* const
            potential_exchange_parameters,
        MaterialLib::Solids::MechanicsBase<DisplacementDim> const&
            solid_material,
        ProcessLib::ThermoRichardsMechanics::MaterialStateData<DisplacementDim>&
            material_state_data)
{
    auto const& liquid_phase = medium->phase("AqueousLiquid");
    auto const& solid_phase = medium->phase("Solid");

    auto const& identity2 = MathLib::KelvinVector::Invariants<
        MathLib::KelvinVector::kelvin_vector_dimensions(
            DisplacementDim)>::identity2;

    double const temperature = T_data();
    double const p_cap_ip = p_cap_data.p_cap;
    double const p_cap_prev_ip = p_cap_data.p_cap_prev;

    auto const& eps = std::get<StrainData<DisplacementDim>>(state_current);
    auto& S_L =
        std::get<ProcessLib::ThermoRichardsMechanics::SaturationData>(state_current).S_L;
    auto const S_L_prev =
        std::get<
            PrevState<ProcessLib::ThermoRichardsMechanics::SaturationData>>(
            state_previous)
            ->S_L;
    auto const alpha =
        medium->property(MPL::PropertyType::biot_coefficient)
            .template value<double>(variables, x_position, t, dt);
    *std::get<ProcessLib::ThermoRichardsMechanics::BiotData>(constitutive_data) = alpha;

    variables.stress =
        std::get<ProcessLib::ConstitutiveRelations::EffectiveStressData<
            DisplacementDim>>(state_current)
            .sigma_eff;
    // Set mechanical strain temporary to compute tangent stiffness.
    variables.mechanical_strain
        .emplace<MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
            eps.eps);
    auto const C_el = ip_data.computeElasticTangentStiffness(
        variables, t, x_position, dt, solid_material,
        *material_state_data.material_state_variables);

    auto const beta_SR =
        (1 - alpha) / solid_material.getBulkModulus(t, x_position, &C_el);
    variables.grain_compressibility = beta_SR;
    std::get<ProcessLib::ThermoRichardsMechanics::SolidCompressibilityData>(constitutive_data)
        .beta_SR = beta_SR;

    auto const rho_LR =
        liquid_phase.property(MPL::PropertyType::density)
            .template value<double>(variables, x_position, t, dt);
    variables.density = rho_LR;
    *std::get<LiquidDensity>(constitutive_data) = rho_LR;

    S_L = medium->property(MPL::PropertyType::saturation)
              .template value<double>(variables, x_position, t, dt);
    variables.liquid_saturation = S_L;
    variables_prev.liquid_saturation = S_L_prev;

    // tangent derivative for Jacobian
    double const dS_L_dp_cap =
        medium->property(MPL::PropertyType::saturation)
            .template dValue<double>(variables,
                                     MPL::Variable::capillary_pressure,
                                     x_position, t, dt);
    std::get<ProcessLib::ThermoRichardsMechanics::SaturationDataDeriv>(constitutive_data)
        .dS_L_dp_cap = dS_L_dp_cap;
    // secant derivative from time discretization for storage
    // use tangent, if secant is not available
    double const DeltaS_L_Deltap_cap =
        (p_cap_ip == p_cap_prev_ip)
            ? dS_L_dp_cap
            : (S_L - S_L_prev) / (p_cap_ip - p_cap_prev_ip);
    std::get<SaturationSecantDerivative>(constitutive_data).DeltaS_L_Deltap_cap =
        DeltaS_L_Deltap_cap;

    auto const chi = [medium, x_position, t, dt](double const S_L)
    {
        MPL::VariableArray vs;
        vs.liquid_saturation = S_L;
        return medium->property(MPL::PropertyType::bishops_effective_stress)
            .template value<double>(vs, x_position, t, dt);
    };
    double const chi_S_L = chi(S_L);
    std::get<ProcessLib::ThermoRichardsMechanics::BishopsData>(constitutive_data).chi_S_L =
        chi_S_L;
    double const chi_S_L_prev = chi(S_L_prev);
    std::get<PrevState<ProcessLib::ThermoRichardsMechanics::BishopsData>>(constitutive_data)
        ->chi_S_L = chi_S_L_prev;

    auto const dchi_dS_L =
        medium->property(MPL::PropertyType::bishops_effective_stress)
            .template dValue<double>(
                variables, MPL::Variable::liquid_saturation, x_position, t, dt);
    std::get<ProcessLib::ThermoRichardsMechanics::BishopsData>(constitutive_data).dchi_dS_L =
        dchi_dS_L;

    double const p_FR = -chi_S_L * p_cap_ip;
    variables.effective_pore_pressure = p_FR;
    variables_prev.effective_pore_pressure = -chi_S_L_prev * p_cap_prev_ip;

    // Set volumetric strain rate for the general case without swelling.
    variables.volumetric_strain = Invariants::trace(eps.eps);
    // TODO (CL) changed that, using eps_prev for the moment, not B * u_prev
    // variables_prev.volumetric_strain = Invariants::trace(B * u_prev);
    variables_prev.volumetric_strain = Invariants::trace(
        std::get<PrevState<StrainData<DisplacementDim>>>(state_previous)->eps);

    auto& phi =
        std::get<ProcessLib::ThermoRichardsMechanics::PorosityData>(state_current).phi;
    {  // Porosity update
        auto const phi_prev =
            std::get<
                PrevState<ProcessLib::ThermoRichardsMechanics::PorosityData>>(
                state_previous)
                ->phi;
        variables_prev.porosity = phi_prev;
        phi = medium->property(MPL::PropertyType::porosity)
                  .template value<double>(variables, variables_prev, x_position,
                                          t, dt);
        variables.porosity = phi;
    }
    std::get<ProcessLib::ThermoRichardsMechanics::PorosityData>(constitutive_data).phi = phi;

    if (alpha < phi)
    {
        auto const eid =
            x_position.getElementID()
                ? static_cast<std::ptrdiff_t>(*x_position.getElementID())
                : static_cast<std::ptrdiff_t>(-1);
        OGS_FATAL(
            "RichardsMechanics: Biot-coefficient {} is smaller than porosity "
            "{} in element {}.",
            alpha, phi, eid);
    }

    auto const mu = liquid_phase.property(MPL::PropertyType::viscosity)
                        .template value<double>(variables, x_position, t, dt);
    *std::get<ProcessLib::ThermoRichardsMechanics::LiquidViscosityData>(constitutive_data) =
        mu;

    {
        // Swelling and possibly volumetric strain rate update.
        auto& sigma_sw =
            std::get<ProcessLib::ThermoRichardsMechanics::
                         ConstitutiveStress_StrainTemperature::
                             SwellingDataStateful<DisplacementDim>>(state_current);
        auto const& sigma_sw_prev =
            std::get<PrevState<ProcessLib::ThermoRichardsMechanics::
                                   ConstitutiveStress_StrainTemperature::
                                       SwellingDataStateful<DisplacementDim>>>(
                state_previous);
        auto const transport_porosity_prev = std::get<PrevState<
            ProcessLib::ThermoRichardsMechanics::TransportPorosityData>>(
            state_previous);
        auto const phi_prev = std::get<
            PrevState<ProcessLib::ThermoRichardsMechanics::PorosityData>>(
            state_previous);
        auto& transport_porosity = std::get<
            ProcessLib::ThermoRichardsMechanics::TransportPorosityData>(state_current);
        auto& p_L_m = std::get<MicroPressure>(state_current);
        auto const p_L_m_prev = std::get<PrevState<MicroPressure>>(state_previous);
        auto& S_L_m = std::get<MicroSaturation>(state_current);
        auto const S_L_m_prev = std::get<PrevState<MicroSaturation>>(state_previous);

        updateSwellingStressAndVolumetricStrain<DisplacementDim>(
            *medium, solid_phase, C_el, rho_LR, mu, micro_porosity_parameters,
            potential_exchange_parameters, alpha, phi, p_cap_ip, variables,
            variables_prev, x_position, t, dt, sigma_sw, sigma_sw_prev,
            transport_porosity_prev, phi_prev, transport_porosity, p_L_m_prev,
            S_L_m_prev, p_L_m, S_L_m);
    }

    auto const transport_porosity_prev_value = std::get<PrevState<
        ProcessLib::ThermoRichardsMechanics::TransportPorosityData>>(state_previous)
                                                    ->phi;
    auto const n_l_prev_value =
        **std::get<PrevState<MicroWaterContent>>(state_previous);

    updateMicroscaleHydraulicState<DisplacementDim>(
        state_current, state_previous, p_cap_ip, rho_LR, mu, dt, variables, variables_prev,
        {.phi = phi,
         .phi_M_prev = transport_porosity_prev_value,
         .phi_m_prev = n_l_prev_value,
         .volumetric_strain = variables.volumetric_strain,
         .volumetric_strain_prev = variables_prev.volumetric_strain},
        micro_porosity_parameters, potential_exchange_parameters);
    updatePorositySplitState<DisplacementDim>(
        state_current, state_previous, phi, variables, variables_prev,
        potential_exchange_parameters);
    updateTotalPorosityState<DisplacementDim>(
        state_current, state_previous, phi, variables, variables_prev,
        potential_exchange_parameters);
    std::get<ProcessLib::ThermoRichardsMechanics::PorosityData>(constitutive_data).phi =
        std::get<ProcessLib::ThermoRichardsMechanics::PorosityData>(state_current).phi;
    updateSwellingState<DisplacementDim>(
        solid_phase, C_el, state_current, state_previous, variables, variables_prev, x_position,
        t, dt, potential_exchange_parameters);

    if (medium->hasProperty(MPL::PropertyType::transport_porosity))
    {
        if (!medium->hasProperty(MPL::PropertyType::saturation_micro) &&
            !isPotentialExchangeEnabled(potential_exchange_parameters))
        {
            auto& transport_porosity =
                std::get<
                    ProcessLib::ThermoRichardsMechanics::TransportPorosityData>(
                    state_current)
                    .phi;
            auto const transport_porosity_prev = std::get<PrevState<
                ProcessLib::ThermoRichardsMechanics::TransportPorosityData>>(
                                                     state_previous)
                                                     ->phi;
            variables_prev.transport_porosity = transport_porosity_prev;

            transport_porosity =
                medium->property(MPL::PropertyType::transport_porosity)
                    .template value<double>(variables, variables_prev,
                                            x_position, t, dt);
            variables.transport_porosity = transport_porosity;
        }
    }
    else
    {
        variables.transport_porosity = phi;
    }

    // Set mechanical variables for the intrinsic permeability model
    // For stress dependent permeability.
    {
        // TODO mechanical constitutive relation will be evaluated afterwards
        auto const sigma_total =
            (std::get<ProcessLib::ConstitutiveRelations::EffectiveStressData<
                 DisplacementDim>>(state_current)
                 .sigma_eff +
             alpha * p_FR * identity2)
                .eval();
        // For stress dependent permeability.
        variables.total_stress.emplace<SymmetricTensor>(
            MathLib::KelvinVector::kelvinVectorToSymmetricTensor(sigma_total));
    }

    variables.equivalent_plastic_strain =
        material_state_data.material_state_variables
            ->getEquivalentPlasticStrain();

    double const k_rel =
        medium->property(MPL::PropertyType::relative_permeability)
            .template value<double>(variables, x_position, t, dt);

    auto const K_intrinsic = MPL::formEigenTensor<DisplacementDim>(
        medium->property(MPL::PropertyType::permeability)
            .value(variables, x_position, t, dt));

    std::get<
        ProcessLib::ThermoRichardsMechanics::PermeabilityData<DisplacementDim>>(
        constitutive_data)
        .k_rel = k_rel;
    std::get<
        ProcessLib::ThermoRichardsMechanics::PermeabilityData<DisplacementDim>>(
        constitutive_data)
        .Ki = K_intrinsic;

    //
    // displacement equation, displacement part
    //

    {
        auto& sigma_sw =
            std::get<ProcessLib::ThermoRichardsMechanics::
                         ConstitutiveStress_StrainTemperature::
                             SwellingDataStateful<DisplacementDim>>(state_current)
                .sigma_sw;

        auto& eps_m =
            std::get<ProcessLib::ConstitutiveRelations::MechanicalStrainData<
                DisplacementDim>>(state_current)
                .eps_m;
        bool const swelling_stress_active =
            solid_phase.hasProperty(MPL::PropertyType::swelling_stress_rate) ||
            isPotentialExchangeEnabled(potential_exchange_parameters);
        eps_m.noalias() =
            swelling_stress_active ? eps.eps + C_el.inverse() * sigma_sw
                                   : eps.eps;
        variables.mechanical_strain
            .emplace<MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
                eps_m);
    }

    {
        auto& sigma_eff =
            std::get<ProcessLib::ConstitutiveRelations::EffectiveStressData<
                DisplacementDim>>(state_current);
        auto const& sigma_eff_prev =
            std::get<PrevState<ProcessLib::ConstitutiveRelations::
                                   EffectiveStressData<DisplacementDim>>>(
                state_previous);
        auto const& eps_m =
            std::get<ProcessLib::ConstitutiveRelations::MechanicalStrainData<
                DisplacementDim>>(state_current);
        auto& eps_m_prev =
            std::get<PrevState<ProcessLib::ConstitutiveRelations::
                                   MechanicalStrainData<DisplacementDim>>>(
                state_previous);

        auto C = ip_data.updateConstitutiveRelation(
            variables, t, x_position, dt, temperature, sigma_eff,
            sigma_eff_prev, eps_m, eps_m_prev, solid_material,
            material_state_data.material_state_variables);

        *std::get<StiffnessTensor<DisplacementDim>>(constitutive_data) = std::move(C);
    }

    // p_SR
    variables.solid_grain_pressure =
        p_FR - std::get<ProcessLib::ConstitutiveRelations::EffectiveStressData<
                   DisplacementDim>>(state_current)
                       .sigma_eff.dot(identity2) /
                   (3 * (1 - phi));
    auto const rho_SR =
        solid_phase.property(MPL::PropertyType::density)
            .template value<double>(variables, x_position, t, dt);

    double const rho = rho_SR * (1 - phi) + S_L * phi * rho_LR;
    *std::get<Density>(constitutive_data) = rho;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
void RichardsMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                     ShapeFunctionPressure, DisplacementDim>::
    assembleWithJacobian(double const t, double const dt,
                         std::vector<double> const& local_x,
                         std::vector<double> const& local_x_prev,
                         std::vector<double>& local_rhs_data,
                         std::vector<double>& local_Jac_data)
{
    assert(local_x.size() == pressure_size + displacement_size);

    auto const [p_L, u] = localDOF(local_x);
    auto const [p_L_prev, u_prev] = localDOF(local_x_prev);

    auto local_Jac = MathLib::createZeroedMatrix<
        typename ShapeMatricesTypeDisplacement::template MatrixType<
            displacement_size + pressure_size,
            displacement_size + pressure_size>>(
        local_Jac_data, displacement_size + pressure_size,
        displacement_size + pressure_size);

    auto local_rhs = MathLib::createZeroedVector<
        typename ShapeMatricesTypeDisplacement::template VectorType<
            displacement_size + pressure_size>>(
        local_rhs_data, displacement_size + pressure_size);

    auto const& identity2 = MathLib::KelvinVector::Invariants<
        MathLib::KelvinVector::kelvin_vector_dimensions(
            DisplacementDim)>::identity2;

    typename ShapeMatricesTypePressure::NodalMatrixType laplace_p =
        ShapeMatricesTypePressure::NodalMatrixType::Zero(pressure_size,
                                                         pressure_size);

    typename ShapeMatricesTypePressure::NodalMatrixType storage_p_a_p =
        ShapeMatricesTypePressure::NodalMatrixType::Zero(pressure_size,
                                                         pressure_size);

    typename ShapeMatricesTypePressure::NodalMatrixType storage_p_a_S_Jpp =
        ShapeMatricesTypePressure::NodalMatrixType::Zero(pressure_size,
                                                         pressure_size);

    typename ShapeMatricesTypePressure::NodalMatrixType storage_p_a_S =
        ShapeMatricesTypePressure::NodalMatrixType::Zero(pressure_size,
                                                         pressure_size);

    typename ShapeMatricesTypeDisplacement::template MatrixType<
        displacement_size, pressure_size>
        Kup = ShapeMatricesTypeDisplacement::template MatrixType<
            displacement_size, pressure_size>::Zero(displacement_size,
                                                    pressure_size);

    typename ShapeMatricesTypeDisplacement::template MatrixType<
        pressure_size, displacement_size>
        Kpu = ShapeMatricesTypeDisplacement::template MatrixType<
            pressure_size, displacement_size>::Zero(pressure_size,
                                                    displacement_size);

    auto const& medium =
        this->process_data_.media_map.getMedium(this->element_.getID());
    auto const& liquid_phase = medium->phase("AqueousLiquid");
    auto const& solid_phase = medium->phase("Solid");
    MPL::VariableArray variables;
    MPL::VariableArray variables_prev;

    unsigned const n_integration_points =
        this->integration_method_.getNumberOfPoints();
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        ConstitutiveData<DisplacementDim> constitutive_data;
        auto& state_current = this->current_states_[ip];
        auto const& state_previous = this->prev_states_[ip];
        [[maybe_unused]] auto models = createConstitutiveModels(
            this->process_data_, this->solid_material_);

        auto const& w = ip_data_[ip].integration_weight;

        auto const& N_u = ip_data_[ip].N_u;
        auto const& dNdx_u = ip_data_[ip].dNdx_u;

        auto const& N_p = ip_data_[ip].N_p;
        auto const& dNdx_p = ip_data_[ip].dNdx_p;

        ParameterLib::SpatialPosition x_position = {
            std::nullopt, this->element_.getID(),
            MathLib::Point3d(
                NumLib::interpolateCoordinates<ShapeFunctionDisplacement,
                                               ShapeMatricesTypeDisplacement>(
                    this->element_, N_u))};
        auto const x_coord = x_position.getCoordinates().value()[0];

        auto const B =
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunctionDisplacement::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                dNdx_u, N_u, x_coord, this->is_axially_symmetric_);

        double p_cap_ip;
        NumLib::shapeFunctionInterpolate(-p_L, N_p, p_cap_ip);

        double p_cap_prev_ip;
        NumLib::shapeFunctionInterpolate(-p_L_prev, N_p, p_cap_prev_ip);

        variables.capillary_pressure = p_cap_ip;
        variables.liquid_phase_pressure = -p_cap_ip;
        // setting pG to 1 atm
        // TODO : rewrite equations s.t. p_L = pG-p_cap
        variables.gas_phase_pressure = 1.0e5;

        auto const temperature =
            medium->property(MPL::PropertyType::reference_temperature)
                .template value<double>(variables, x_position, t, dt);
        variables.temperature = temperature;

        std::get<StrainData<DisplacementDim>>(state_current).eps.noalias() = B * u;

        assembleWithJacobianEvalConstitutiveSetting(
            t, dt, x_position, ip_data_[ip], variables, variables_prev, medium,
            TemperatureData{temperature},
            CapillaryPressureData<DisplacementDim>{
                p_cap_ip, p_cap_prev_ip,
                Eigen::Vector<double, DisplacementDim>::Zero()},
            constitutive_data, state_current, state_previous, this->process_data_.micro_porosity_parameters,
            this->getPotentialExchangeParameters(),
            this->solid_material_, this->material_states_[ip]);

        {
            auto const& C = *std::get<StiffnessTensor<DisplacementDim>>(constitutive_data);
            local_Jac
                .template block<displacement_size, displacement_size>(
                    displacement_index, displacement_index)
                .noalias() += B.transpose() * C * B * w;
        }

        auto const& b = this->process_data_.specific_body_force;

        {
            auto const& sigma_eff =
                std::get<ProcessLib::ConstitutiveRelations::EffectiveStressData<
                    DisplacementDim>>(this->current_states_[ip])
                    .sigma_eff;
            double const rho = *std::get<Density>(constitutive_data);
            local_rhs.template segment<displacement_size>(displacement_index)
                .noalias() -= (B.transpose() * sigma_eff -
                               N_u_op(N_u).transpose() * rho * b) *
                              w;
        }

        //
        // displacement equation, pressure part
        //

        double const alpha =
            *std::get<ProcessLib::ThermoRichardsMechanics::BiotData>(constitutive_data);
        double const dS_L_dp_cap =
            std::get<ProcessLib::ThermoRichardsMechanics::SaturationDataDeriv>(
                constitutive_data)
                .dS_L_dp_cap;

        {
            double const chi_S_L =
                std::get<ProcessLib::ThermoRichardsMechanics::BishopsData>(constitutive_data)
                    .chi_S_L;
            Kup.noalias() +=
                B.transpose() * alpha * chi_S_L * identity2 * N_p * w;
            double const dchi_dS_L =
                std::get<ProcessLib::ThermoRichardsMechanics::BishopsData>(constitutive_data)
                    .dchi_dS_L;

            local_Jac
                .template block<displacement_size, pressure_size>(
                    displacement_index, pressure_index)
                .noalias() -= B.transpose() * alpha *
                              (chi_S_L + dchi_dS_L * p_cap_ip * dS_L_dp_cap) *
                              identity2 * N_p * w;
        }

        double const phi =
            std::get<ProcessLib::ThermoRichardsMechanics::PorosityData>(constitutive_data).phi;
        double const rho_LR = *std::get<LiquidDensity>(constitutive_data);
        local_Jac
            .template block<displacement_size, pressure_size>(
                displacement_index, pressure_index)
            .noalias() +=
            N_u_op(N_u).transpose() * phi * rho_LR * dS_L_dp_cap * b * N_p * w;

        // For the swelling stress with double structure model the corresponding
        // Jacobian u-p entry would be required, but it does not improve
        // convergence and sometimes worsens it:
        // if (medium->hasProperty(MPL::PropertyType::saturation_micro))
        // {
        //     -B.transpose() *
        //         dsigma_sw_dS_L_m* dS_L_m_dp_cap_m*(p_L_m - p_L_m_prev) /
        //         (p_cap_ip - p_cap_prev_ip) * N_p* w;
        // }
        if (!medium->hasProperty(MPL::PropertyType::saturation_micro) &&
            !isPotentialExchangeEnabled(
                this->getPotentialExchangeParameters()) &&
            solid_phase.hasProperty(MPL::PropertyType::swelling_stress_rate))
        {
            using DimMatrix = Eigen::Matrix<double, 3, 3>;
            auto const dsigma_sw_dS_L =
                MathLib::KelvinVector::tensorToKelvin<DisplacementDim>(
                    solid_phase
                        .property(MPL::PropertyType::swelling_stress_rate)
                        .template dValue<DimMatrix>(
                            variables, variables_prev,
                            MPL::Variable::liquid_saturation, x_position, t,
                            dt));
            local_Jac
                .template block<displacement_size, pressure_size>(
                    displacement_index, pressure_index)
                .noalias() +=
                B.transpose() * dsigma_sw_dS_L * dS_L_dp_cap * N_p * w;
        }
        //
        // pressure equation, displacement part.
        //
        double const S_L =
            std::get<ProcessLib::ThermoRichardsMechanics::SaturationData>(
                this->current_states_[ip])
                .S_L;
        if (this->process_data_.explicit_hm_coupling_in_unsaturated_zone)
        {
            double const chi_S_L_prev = std::get<PrevState<
                ProcessLib::ThermoRichardsMechanics::BishopsData>>(constitutive_data)
                                            ->chi_S_L;
            Kpu.noalias() += N_p.transpose() * chi_S_L_prev * rho_LR * alpha *
                             identity2.transpose() * B * w;
        }
        else
        {
            Kpu.noalias() += N_p.transpose() * S_L * rho_LR * alpha *
                             identity2.transpose() * B * w;
        }

        //
        // pressure equation, pressure part.
        //

        double const k_rel =
            std::get<ProcessLib::ThermoRichardsMechanics::PermeabilityData<
                DisplacementDim>>(constitutive_data)
                .k_rel;
        auto const& K_intrinsic =
            std::get<ProcessLib::ThermoRichardsMechanics::PermeabilityData<
                DisplacementDim>>(constitutive_data)
                .Ki;
        double const mu =
            *std::get<ProcessLib::ThermoRichardsMechanics::LiquidViscosityData>(
                constitutive_data);

        GlobalDimMatrixType const rho_Ki_over_mu = K_intrinsic * rho_LR / mu;

        laplace_p.noalias() +=
            dNdx_p.transpose() * k_rel * rho_Ki_over_mu * dNdx_p * w;

        auto const beta_LR =
            1 / rho_LR *
            liquid_phase.property(MPL::PropertyType::density)
                .template dValue<double>(variables,
                                         MPL::Variable::liquid_phase_pressure,
                                         x_position, t, dt);

        double const beta_SR =
            std::get<
                ProcessLib::ThermoRichardsMechanics::SolidCompressibilityData>(
                constitutive_data)
                .beta_SR;
        double const a0 = (alpha - phi) * beta_SR;
        double const specific_storage_a_p = S_L * (phi * beta_LR + S_L * a0);
        double const specific_storage_a_S = phi - p_cap_ip * S_L * a0;

        double const dspecific_storage_a_p_dp_cap =
            dS_L_dp_cap * (phi * beta_LR + 2 * S_L * a0);
        double const dspecific_storage_a_S_dp_cap =
            -a0 * (S_L + p_cap_ip * dS_L_dp_cap);

        storage_p_a_p.noalias() +=
            N_p.transpose() * rho_LR * specific_storage_a_p * N_p * w;

        double const DeltaS_L_Deltap_cap =
            std::get<SaturationSecantDerivative>(constitutive_data).DeltaS_L_Deltap_cap;
        storage_p_a_S.noalias() -= N_p.transpose() * rho_LR *
                                   specific_storage_a_S * DeltaS_L_Deltap_cap *
                                   N_p * w;

        local_Jac
            .template block<pressure_size, pressure_size>(pressure_index,
                                                          pressure_index)
            .noalias() += N_p.transpose() * (p_cap_ip - p_cap_prev_ip) / dt *
                          rho_LR * dspecific_storage_a_p_dp_cap * N_p * w;

        double const S_L_prev =
            std::get<
                PrevState<ProcessLib::ThermoRichardsMechanics::SaturationData>>(
                this->prev_states_[ip])
                ->S_L;
        storage_p_a_S_Jpp.noalias() -=
            N_p.transpose() * rho_LR *
            ((S_L - S_L_prev) * dspecific_storage_a_S_dp_cap +
             specific_storage_a_S * dS_L_dp_cap) /
            dt * N_p * w;

        if (!this->process_data_.explicit_hm_coupling_in_unsaturated_zone)
        {
            local_Jac
                .template block<pressure_size, pressure_size>(pressure_index,
                                                              pressure_index)
                .noalias() -= N_p.transpose() * rho_LR * dS_L_dp_cap * alpha *
                              identity2.transpose() * B * (u - u_prev) / dt *
                              N_p * w;
        }

        double const dk_rel_dS_l =
            medium->property(MPL::PropertyType::relative_permeability)
                .template dValue<double>(variables,
                                         MPL::Variable::liquid_saturation,
                                         x_position, t, dt);
        typename ShapeMatricesTypeDisplacement::GlobalDimVectorType const
            grad_p_cap = -dNdx_p * p_L;
        local_Jac
            .template block<pressure_size, pressure_size>(pressure_index,
                                                          pressure_index)
            .noalias() += dNdx_p.transpose() * rho_Ki_over_mu * grad_p_cap *
                          dk_rel_dS_l * dS_L_dp_cap * N_p * w;

        local_Jac
            .template block<pressure_size, pressure_size>(pressure_index,
                                                          pressure_index)
            .noalias() += dNdx_p.transpose() * rho_LR * rho_Ki_over_mu * b *
                          dk_rel_dS_l * dS_L_dp_cap * N_p * w;

        local_rhs.template segment<pressure_size>(pressure_index).noalias() +=
            dNdx_p.transpose() * rho_LR * k_rel * rho_Ki_over_mu * b * w;

        auto const* const potential_exchange_params_ptr =
            this->getPotentialExchangeParameters();
        bool const potential_exchange_enabled =
            isPotentialExchangeEnabled(potential_exchange_params_ptr);
        if ((medium->hasProperty(MPL::PropertyType::saturation_micro) ||
             potential_exchange_enabled) &&
            this->process_data_.micro_porosity_parameters)
        {
            double const alpha_bar =
                this->process_data_.micro_porosity_parameters
                    ->mass_exchange_coefficient;
            auto const p_L_m =
                *std::get<MicroPressure>(this->current_states_[ip]);
            double const p_L_ip = -p_cap_ip;
            double const pressure_tolerance =
                getPotentialPressureTolerance(
                    potential_exchange_params_ptr);

            bool use_vdw_micro_potential_for_active_exchange = false;
            double mu_lR_vdw = 0.0;
            double dmu_lR_vdw_drho_lR = 0.0;
            bool use_custom_dmu_lR_vdw_dpL = false;
            double dmu_lR_vdw_dpL = 0.0;
            bool use_fd_jacobian_for_direct_macro_derivative = false;
            double fd_jacobian_perturbation = 1e-8;
            if (potential_exchange_enabled)
            {
                auto const n_l =
                    std::max(1e-16,
                             *std::get<MicroWaterContent>(
                                 this->current_states_[ip]));
                auto const transport_porosity_prev =
                    std::get<PrevState<ProcessLib::ThermoRichardsMechanics::
                                           TransportPorosityData>>(
                        this->prev_states_[ip])
                        ->phi;
                auto const n_l_prev = **std::get<PrevState<MicroWaterContent>>(
                    this->prev_states_[ip]);
                PotentialExchangeLocalSolveContext const local_solve_context{
                    .phi = phi,
                    .phi_M_prev = transport_porosity_prev,
                    .phi_m_prev = n_l_prev,
                    .volumetric_strain = variables.volumetric_strain,
                    .volumetric_strain_prev = variables_prev.volumetric_strain};
                auto const micro_potential = computeActiveMicroPotential(
                    n_l, rho_LR, local_solve_context,
                    *potential_exchange_params_ptr);
                use_vdw_micro_potential_for_active_exchange = true;
                mu_lR_vdw = micro_potential.mu_lR;
                dmu_lR_vdw_drho_lR = micro_potential.dmu_lR_drho_lR;
                use_fd_jacobian_for_direct_macro_derivative =
                    potential_exchange_params_ptr
                        ->use_fd_jacobian_for_exchange;
                fd_jacobian_perturbation =
                    potential_exchange_params_ptr->fd_jacobian_perturbation;

                // In analytic mode, include implicit n_l(p_L) chain coupling
                // in dmu_lR/dp_L for the active exchange Jacobian term.
                if (!use_fd_jacobian_for_direct_macro_derivative)
                {
                    double const drho_LR_dpL = rho_LR * beta_LR;
                    auto const macro_potential = computeYoungLaplaceMacroPotential(
                        p_L_ip, rho_LR, pressure_tolerance);
                    double const alpha_M_effective = alpha_bar * rho_LR / mu;
                    auto const exchange = computePotentialDrivenMassExchange(
                        alpha_M_effective, macro_potential.mu_LR,
                        micro_potential.mu_lR);
                    double const dn_l_dpL = computeImplicitNlDpL(
                        n_l_prev, p_L_ip, dt, rho_LR, drho_LR_dpL, alpha_bar, mu,
                        macro_potential, micro_potential, exchange,
                        local_solve_context,
                        *potential_exchange_params_ptr);

                    dmu_lR_vdw_dpL = micro_potential.dmu_lR_dnl * dn_l_dpL +
                                     micro_potential.dmu_lR_drho_lR *
                                         drho_LR_dpL;
                    use_custom_dmu_lR_vdw_dpL = true;

                    auto const n_l_prev =
                        std::max(1e-16,
                                 **std::get<PrevState<MicroWaterContent>>(
                                     this->prev_states_[ip]));
                    if (potential_exchange_params_ptr->check_local_jacobian)
                    {
                        auto const fd_data =
                            computeLocalJacobianDiagnosticData(
                                n_l_prev, p_L_ip, dt, rho_LR, drho_LR_dpL,
                                alpha_bar, mu, pressure_tolerance,
                                local_solve_context,
                                *potential_exchange_params_ptr);
                        auto const analytic_drho_L_hat_dpL =
                            -exchange.drho_l_hat_dalpha_M *
                                (alpha_bar / mu * drho_LR_dpL) -
                            exchange.drho_l_hat_dmu_LR *
                                (macro_potential.dmu_LR_dpLR +
                                 macro_potential.dmu_LR_drho_LR *
                                     drho_LR_dpL) -
                            exchange.drho_l_hat_dmu_lR * dmu_lR_vdw_dpL;
                        maybeLogLocalJacobianDiagnostic(
                            p_L_ip, n_l_prev, n_l, dn_l_dpL,
                            analytic_drho_L_hat_dpL, fd_data,
                            *potential_exchange_params_ptr);
                    }
                }
            }

            auto const potential_exchange_result = computePotentialExchangeUpdate(
                alpha_bar, mu, p_L_ip, p_L_m, rho_LR, beta_LR,
                pressure_tolerance, potential_exchange_enabled,
                use_vdw_micro_potential_for_active_exchange, mu_lR_vdw,
                dmu_lR_vdw_drho_lR, use_custom_dmu_lR_vdw_dpL,
                dmu_lR_vdw_dpL,
                use_fd_jacobian_for_direct_macro_derivative,
                fd_jacobian_perturbation);
            // Phase 2C: activate potential-driven exchange source in the macro
            // balance.
            local_rhs.template segment<pressure_size>(pressure_index)
                .noalias() += N_p.transpose() * potential_exchange_result.exchange.rho_L_hat *
                              w;

            // Direct macro Jacobian term for the exchange source. In analytic
            // mode this includes the implicit n_l(p_L) chain contribution.
            local_Jac
                .template block<pressure_size, pressure_size>(pressure_index,
                                                              pressure_index)
                .noalias() -= N_p.transpose() *
                              potential_exchange_result.drho_L_hat_dpL_direct * N_p * w;

            // Keep the microscale pressure-state sensitivity lagged via the
            // secant term only in the placeholder microscale path. In the
            // vdW+ n_l opt-in path this term is intentionally omitted because
            // the active microscale potential is no longer p_L_m/rho_LR.
            if (!use_vdw_micro_potential_for_active_exchange &&
                p_cap_ip != p_cap_prev_ip)
            {
                auto const p_L_m_prev = **std::get<PrevState<MicroPressure>>(
                    this->prev_states_[ip]);
                local_Jac
                    .template block<pressure_size, pressure_size>(
                        pressure_index, pressure_index)
                    .noalias() += N_p.transpose() * alpha_bar / mu *
                                  (p_L_m - p_L_m_prev) /
                                  (p_cap_ip - p_cap_prev_ip) * N_p * w;
            }
        }
    }

    if (this->process_data_.apply_mass_lumping)
    {
        storage_p_a_p = storage_p_a_p.colwise().sum().eval().asDiagonal();
        storage_p_a_S = storage_p_a_S.colwise().sum().eval().asDiagonal();
        storage_p_a_S_Jpp =
            storage_p_a_S_Jpp.colwise().sum().eval().asDiagonal();
    }

    // pressure equation, pressure part.
    local_Jac
        .template block<pressure_size, pressure_size>(pressure_index,
                                                      pressure_index)
        .noalias() += laplace_p + storage_p_a_p / dt + storage_p_a_S_Jpp;

    // pressure equation, displacement part.
    local_Jac
        .template block<pressure_size, displacement_size>(pressure_index,
                                                          displacement_index)
        .noalias() = Kpu / dt;

    // pressure equation
    local_rhs.template segment<pressure_size>(pressure_index).noalias() -=
        laplace_p * p_L +
        (storage_p_a_p + storage_p_a_S) * (p_L - p_L_prev) / dt +
        Kpu * (u - u_prev) / dt;

    // displacement equation
    local_rhs.template segment<displacement_size>(displacement_index)
        .noalias() += Kup * p_L;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
void RichardsMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                     ShapeFunctionPressure, DisplacementDim>::
    assembleWithJacobianForPressureEquations(
        const double /*t*/, double const /*dt*/,
        Eigen::VectorXd const& /*local_x*/,
        Eigen::VectorXd const& /*local_x_prev*/,
        std::vector<double>& /*local_b_data*/,
        std::vector<double>& /*local_Jac_data*/)
{
    OGS_FATAL("RichardsMechanics; The staggered scheme is not implemented.");
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
void RichardsMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                     ShapeFunctionPressure, DisplacementDim>::
    assembleWithJacobianForDeformationEquations(
        const double /*t*/, double const /*dt*/,
        Eigen::VectorXd const& /*local_x*/,
        Eigen::VectorXd const& /*local_x_prev*/,
        std::vector<double>& /*local_b_data*/,
        std::vector<double>& /*local_Jac_data*/)
{
    OGS_FATAL("RichardsMechanics; The staggered scheme is not implemented.");
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
void RichardsMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                     ShapeFunctionPressure, DisplacementDim>::
    assembleWithJacobianForStaggeredScheme(double const t, double const dt,
                                           Eigen::VectorXd const& local_x,
                                           Eigen::VectorXd const& local_x_prev,
                                           int const process_id,
                                           std::vector<double>& local_b_data,
                                           std::vector<double>& local_Jac_data)
{
    // For the equations with pressure
    if (process_id == 0)
    {
        assembleWithJacobianForPressureEquations(t, dt, local_x, local_x_prev,
                                                 local_b_data, local_Jac_data);
        return;
    }

    // For the equations with deformation
    assembleWithJacobianForDeformationEquations(t, dt, local_x, local_x_prev,
                                                local_b_data, local_Jac_data);
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          int DisplacementDim>
void RichardsMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                     ShapeFunctionPressure, DisplacementDim>::
    computeSecondaryVariableConcrete(double const t, double const dt,
                                     Eigen::VectorXd const& local_x,
                                     Eigen::VectorXd const& local_x_prev)
{
    auto const [p_L, u] = localDOF(local_x);
    auto const [p_L_prev, u_prev] = localDOF(local_x_prev);

    auto const& identity2 = MathLib::KelvinVector::Invariants<
        MathLib::KelvinVector::kelvin_vector_dimensions(
            DisplacementDim)>::identity2;

    auto const& medium =
        this->process_data_.media_map.getMedium(this->element_.getID());
    auto const& liquid_phase = medium->phase("AqueousLiquid");
    auto const& solid_phase = medium->phase("Solid");
    MPL::VariableArray variables;
    MPL::VariableArray variables_prev;

    unsigned const n_integration_points =
        this->integration_method_.getNumberOfPoints();

    double saturation_avg = 0;
    double porosity_avg = 0;

    using KV = MathLib::KelvinVector::KelvinVectorType<DisplacementDim>;
    KV sigma_avg = KV::Zero();

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto const& N_p = ip_data_[ip].N_p;
        auto const& N_u = ip_data_[ip].N_u;
        auto const& dNdx_u = ip_data_[ip].dNdx_u;

        ParameterLib::SpatialPosition x_position = {
            std::nullopt, this->element_.getID(),
            MathLib::Point3d(
                NumLib::interpolateCoordinates<ShapeFunctionDisplacement,
                                               ShapeMatricesTypeDisplacement>(
                    this->element_, N_u))};
        auto const x_coord = x_position.getCoordinates().value()[0];

        auto const B =
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunctionDisplacement::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                dNdx_u, N_u, x_coord, this->is_axially_symmetric_);

        double p_cap_ip;
        NumLib::shapeFunctionInterpolate(-p_L, N_p, p_cap_ip);

        double p_cap_prev_ip;
        NumLib::shapeFunctionInterpolate(-p_L_prev, N_p, p_cap_prev_ip);

        variables.capillary_pressure = p_cap_ip;
        variables.liquid_phase_pressure = -p_cap_ip;
        // setting pG to 1 atm
        // TODO : rewrite equations s.t. p_L = pG-p_cap
        variables.gas_phase_pressure = 1.0e5;

        auto const temperature =
            medium->property(MPL::PropertyType::reference_temperature)
                .template value<double>(variables, x_position, t, dt);
        variables.temperature = temperature;

        auto& eps =
            std::get<StrainData<DisplacementDim>>(this->current_states_[ip])
                .eps;
        eps.noalias() = B * u;
        auto& S_L =
            std::get<ProcessLib::ThermoRichardsMechanics::SaturationData>(
                this->current_states_[ip])
                .S_L;
        auto const S_L_prev =
            std::get<
                PrevState<ProcessLib::ThermoRichardsMechanics::SaturationData>>(
                this->prev_states_[ip])
                ->S_L;
        S_L = medium->property(MPL::PropertyType::saturation)
                  .template value<double>(variables, x_position, t, dt);
        variables.liquid_saturation = S_L;
        variables_prev.liquid_saturation = S_L_prev;

        auto const chi = [medium, x_position, t, dt](double const S_L)
        {
            MPL::VariableArray vs;
            vs.liquid_saturation = S_L;
            return medium->property(MPL::PropertyType::bishops_effective_stress)
                .template value<double>(vs, x_position, t, dt);
        };
        double const chi_S_L = chi(S_L);
        double const chi_S_L_prev = chi(S_L_prev);

        auto const alpha =
            medium->property(MPL::PropertyType::biot_coefficient)
                .template value<double>(variables, x_position, t, dt);
        auto& state_current = this->current_states_[ip];
        variables.stress =
            std::get<ProcessLib::ConstitutiveRelations::EffectiveStressData<
                DisplacementDim>>(state_current)
                .sigma_eff;
        // Set mechanical strain temporary to compute tangent stiffness.
        variables.mechanical_strain
            .emplace<MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
                eps);
        auto const C_el = ip_data_[ip].computeElasticTangentStiffness(
            variables, t, x_position, dt, this->solid_material_,
            *this->material_states_[ip].material_state_variables);

        auto const beta_SR = (1 - alpha) / this->solid_material_.getBulkModulus(
                                               t, x_position, &C_el);
        variables.grain_compressibility = beta_SR;

        variables.effective_pore_pressure = -chi_S_L * p_cap_ip;
        variables_prev.effective_pore_pressure = -chi_S_L_prev * p_cap_prev_ip;

        // Set volumetric strain rate for the general case without swelling.
        variables.volumetric_strain = Invariants::trace(eps);
        variables_prev.volumetric_strain = Invariants::trace(B * u_prev);

        auto& phi = std::get<ProcessLib::ThermoRichardsMechanics::PorosityData>(
                        this->current_states_[ip])
                        .phi;
        {  // Porosity update
            auto const phi_prev = std::get<PrevState<
                ProcessLib::ThermoRichardsMechanics::PorosityData>>(
                                      this->prev_states_[ip])
                                      ->phi;
            variables_prev.porosity = phi_prev;
            phi = medium->property(MPL::PropertyType::porosity)
                      .template value<double>(variables, variables_prev,
                                              x_position, t, dt);
            variables.porosity = phi;
        }

        auto const rho_LR =
            liquid_phase.property(MPL::PropertyType::density)
                .template value<double>(variables, x_position, t, dt);
        variables.density = rho_LR;
        auto const mu =
            liquid_phase.property(MPL::PropertyType::viscosity)
                .template value<double>(variables, x_position, t, dt);

        {
            // Swelling and possibly volumetric strain rate update.
            auto& sigma_sw =
                std::get<ProcessLib::ThermoRichardsMechanics::
                             ConstitutiveStress_StrainTemperature::
                                 SwellingDataStateful<DisplacementDim>>(
                    this->current_states_[ip]);
            auto const& sigma_sw_prev = std::get<
                PrevState<ProcessLib::ThermoRichardsMechanics::
                              ConstitutiveStress_StrainTemperature::
                                  SwellingDataStateful<DisplacementDim>>>(
                this->prev_states_[ip]);
            auto const transport_porosity_prev = std::get<PrevState<
                ProcessLib::ThermoRichardsMechanics::TransportPorosityData>>(
                this->prev_states_[ip]);
            auto const phi_prev = std::get<
                PrevState<ProcessLib::ThermoRichardsMechanics::PorosityData>>(
                this->prev_states_[ip]);
            auto& transport_porosity = std::get<
                ProcessLib::ThermoRichardsMechanics::TransportPorosityData>(
                this->current_states_[ip]);
            auto& p_L_m = std::get<MicroPressure>(this->current_states_[ip]);
            auto const p_L_m_prev =
                std::get<PrevState<MicroPressure>>(this->prev_states_[ip]);
            auto& S_L_m = std::get<MicroSaturation>(this->current_states_[ip]);
            auto const S_L_m_prev =
                std::get<PrevState<MicroSaturation>>(this->prev_states_[ip]);

            updateSwellingStressAndVolumetricStrain<DisplacementDim>(
                *medium, solid_phase, C_el, rho_LR, mu,
                this->process_data_.micro_porosity_parameters,
                this->getPotentialExchangeParameters(), alpha, phi, p_cap_ip,
                variables, variables_prev, x_position, t, dt, sigma_sw,
                sigma_sw_prev, transport_porosity_prev, phi_prev,
                transport_porosity, p_L_m_prev, S_L_m_prev, p_L_m, S_L_m);
        }

        auto const transport_porosity_prev_value = std::get<PrevState<
            ProcessLib::ThermoRichardsMechanics::TransportPorosityData>>(
            this->prev_states_[ip])
                                                        ->phi;
        auto const n_l_prev_value =
            **std::get<PrevState<MicroWaterContent>>(
                this->prev_states_[ip]);

        updateMicroscaleHydraulicState<DisplacementDim>(
            this->current_states_[ip], this->prev_states_[ip], p_cap_ip,
            rho_LR, mu, dt, variables, variables_prev,
            {.phi = phi,
             .phi_M_prev = transport_porosity_prev_value,
             .phi_m_prev = n_l_prev_value,
             .volumetric_strain = variables.volumetric_strain,
             .volumetric_strain_prev = variables_prev.volumetric_strain},
            this->process_data_.micro_porosity_parameters,
            this->getPotentialExchangeParameters());
        updatePorositySplitState<DisplacementDim>(
            this->current_states_[ip], this->prev_states_[ip], phi, variables,
            variables_prev, this->getPotentialExchangeParameters());
        updateTotalPorosityState<DisplacementDim>(
            this->current_states_[ip], this->prev_states_[ip], phi, variables,
            variables_prev, this->getPotentialExchangeParameters());
        updateSwellingState<DisplacementDim>(
            solid_phase, C_el, this->current_states_[ip],
            this->prev_states_[ip], variables, variables_prev, x_position, t,
            dt, this->getPotentialExchangeParameters());

        if (medium->hasProperty(MPL::PropertyType::transport_porosity))
        {
            if (!medium->hasProperty(MPL::PropertyType::saturation_micro) &&
                !isPotentialExchangeEnabled(
                    this->getPotentialExchangeParameters()))
            {
                auto& transport_porosity =
                    std::get<ProcessLib::ThermoRichardsMechanics::
                                 TransportPorosityData>(
                        this->current_states_[ip])
                        .phi;
                auto const transport_porosity_prev =
                    std::get<PrevState<ProcessLib::ThermoRichardsMechanics::
                                           TransportPorosityData>>(
                        this->prev_states_[ip])
                        ->phi;

                variables_prev.transport_porosity = transport_porosity_prev;

                transport_porosity =
                    medium->property(MPL::PropertyType::transport_porosity)
                        .template value<double>(variables, variables_prev,
                                                x_position, t, dt);
                variables.transport_porosity = transport_porosity;
            }
        }
        else
        {
            variables.transport_porosity = phi;
        }

        auto const& sigma_eff =
            std::get<ProcessLib::ConstitutiveRelations::EffectiveStressData<
                DisplacementDim>>(this->current_states_[ip])
                .sigma_eff;

        // Set mechanical variables for the intrinsic permeability model
        // For stress dependent permeability.
        {
            auto const sigma_total =
                (sigma_eff + alpha * chi_S_L * identity2 * p_cap_ip).eval();
            // For stress dependent permeability.
            variables.total_stress.emplace<SymmetricTensor>(
                MathLib::KelvinVector::kelvinVectorToSymmetricTensor(
                    sigma_total));
        }

        variables.equivalent_plastic_strain =
            this->material_states_[ip]
                .material_state_variables->getEquivalentPlasticStrain();

        auto const K_intrinsic = MPL::formEigenTensor<DisplacementDim>(
            medium->property(MPL::PropertyType::permeability)
                .value(variables, x_position, t, dt));

        double const k_rel =
            medium->property(MPL::PropertyType::relative_permeability)
                .template value<double>(variables, x_position, t, dt);

        GlobalDimMatrixType const K_over_mu = k_rel * K_intrinsic / mu;

        double const p_FR = -chi_S_L * p_cap_ip;
        // p_SR
        variables.solid_grain_pressure =
            p_FR - sigma_eff.dot(identity2) / (3 * (1 - phi));
        auto const rho_SR =
            solid_phase.property(MPL::PropertyType::density)
                .template value<double>(variables, x_position, t, dt);
        *std::get<DrySolidDensity>(this->output_data_[ip]) = (1 - phi) * rho_SR;

        {
            auto& state_current = this->current_states_[ip];
            auto const& sigma_sw =
                std::get<ProcessLib::ThermoRichardsMechanics::
                             ConstitutiveStress_StrainTemperature::
                                 SwellingDataStateful<DisplacementDim>>(state_current)
                    .sigma_sw;
            auto& eps_m =
                std::get<ProcessLib::ConstitutiveRelations::
                             MechanicalStrainData<DisplacementDim>>(state_current)
                    .eps_m;
            bool const swelling_stress_active =
                solid_phase.hasProperty(MPL::PropertyType::swelling_stress_rate) ||
                isPotentialExchangeEnabled(
                    this->getPotentialExchangeParameters());
            eps_m.noalias() = swelling_stress_active
                                  ? eps + C_el.inverse() * sigma_sw
                                  : eps;
            variables.mechanical_strain.emplace<
                MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
                eps_m);
        }

        {
            auto& state_current = this->current_states_[ip];
            auto const& state_previous = this->prev_states_[ip];
            auto& sigma_eff =
                std::get<ProcessLib::ConstitutiveRelations::EffectiveStressData<
                    DisplacementDim>>(state_current);
            auto const& sigma_eff_prev =
                std::get<PrevState<ProcessLib::ConstitutiveRelations::
                                       EffectiveStressData<DisplacementDim>>>(
                    state_previous);
            auto const& eps_m =
                std::get<ProcessLib::ConstitutiveRelations::
                             MechanicalStrainData<DisplacementDim>>(state_current);
            auto const& eps_m_prev =
                std::get<PrevState<ProcessLib::ConstitutiveRelations::
                                       MechanicalStrainData<DisplacementDim>>>(
                    state_previous);

            ip_data_[ip].updateConstitutiveRelation(
                variables, t, x_position, dt, temperature, sigma_eff,
                sigma_eff_prev, eps_m, eps_m_prev, this->solid_material_,
                this->material_states_[ip].material_state_variables);
        }

        auto const& b = this->process_data_.specific_body_force;

        // Compute the velocity
        auto const& dNdx_p = ip_data_[ip].dNdx_p;
        std::get<
            ProcessLib::ThermoRichardsMechanics::DarcyLawData<DisplacementDim>>(
            this->output_data_[ip])
            ->noalias() = -K_over_mu * dNdx_p * p_L + rho_LR * K_over_mu * b;

        saturation_avg += S_L;
        porosity_avg += phi;
        sigma_avg += sigma_eff;
    }
    saturation_avg /= n_integration_points;
    porosity_avg /= n_integration_points;
    sigma_avg /= n_integration_points;

    (*this->process_data_.element_saturation)[this->element_.getID()] =
        saturation_avg;
    (*this->process_data_.element_porosity)[this->element_.getID()] =
        porosity_avg;

    Eigen::Map<KV>(
        &(*this->process_data_.element_stresses)[this->element_.getID() *
                                                 KV::RowsAtCompileTime]) =
        MathLib::KelvinVector::kelvinVectorToSymmetricTensor(sigma_avg);

    NumLib::interpolateToHigherOrderNodes<
        ShapeFunctionPressure, typename ShapeFunctionDisplacement::MeshElement,
        DisplacementDim>(this->element_, this->is_axially_symmetric_, p_L,
                         *this->process_data_.pressure_interpolated);
}
}  // namespace RichardsMechanics
}  // namespace ProcessLib
