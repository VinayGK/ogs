// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <optional>

namespace ProcessLib::RichardsMechanics
{
enum class VKPotentialExchangeMode
{
    FullPotential
};

enum class VKMicroPotentialConvention
{
    PositiveReduced,
    NegativeAttractive
};

enum class VKLocalNonlinearSolveMode
{
    ScalarExchange,
    ScalarNotebookStorage
};

enum class VKMacroPorosityUpdateMode
{
    AlgebraicSplit,
    NotebookAdditiveRate
};

enum class VKMicroSolidVolumeFractionMode
{
    Reference,
    CurrentPorositySplit
};

inline constexpr char const* toString(
    VKMicroPotentialConvention const convention)
{
    switch (convention)
    {
        case VKMicroPotentialConvention::PositiveReduced:
            return "positive_reduced";
        case VKMicroPotentialConvention::NegativeAttractive:
            return "negative_attractive";
    }
    return "unknown";
}

inline constexpr double microPotentialSignFactor(
    VKMicroPotentialConvention const convention)
{
    return convention == VKMicroPotentialConvention::NegativeAttractive ? -1.0
                                                                       : 1.0;
}

inline constexpr char const* toString(VKLocalNonlinearSolveMode const mode)
{
    switch (mode)
    {
        case VKLocalNonlinearSolveMode::ScalarExchange:
            return "scalar_exchange";
        case VKLocalNonlinearSolveMode::ScalarNotebookStorage:
            return "scalar_notebook_storage";
    }
    return "unknown";
}

inline constexpr char const* toString(VKMacroPorosityUpdateMode const mode)
{
    switch (mode)
    {
        case VKMacroPorosityUpdateMode::AlgebraicSplit:
            return "algebraic_split";
        case VKMacroPorosityUpdateMode::NotebookAdditiveRate:
            return "notebook_additive_rate";
    }
    return "unknown";
}

inline constexpr char const* toString(
    VKMicroSolidVolumeFractionMode const mode)
{
    switch (mode)
    {
        case VKMicroSolidVolumeFractionMode::Reference:
            return "reference";
        case VKMicroSolidVolumeFractionMode::CurrentPorositySplit:
            return "current_porosity_split";
    }
    return "unknown";
}

struct VKPotentialExchangeParameters
{
    bool enabled = false;
    VKPotentialExchangeMode mode = VKPotentialExchangeMode::FullPotential;

    // Young-Laplace macro potential branch tolerance.
    double pressure_tolerance = 0.0;

    // vdW microscale potential parameters / reference state constants.
    double hamaker_constant = 0.0;
    double specific_surface = 0.0;
    double micro_solid_density_reference = 0.0;          // rho_SR
    double micro_solid_volume_fraction_reference = 0.0;  // n_S
    VKMicroPotentialConvention micro_potential_convention =
        VKMicroPotentialConvention::PositiveReduced;
    VKLocalNonlinearSolveMode local_nonlinear_solve_mode =
        VKLocalNonlinearSolveMode::ScalarExchange;
    VKMacroPorosityUpdateMode macro_porosity_update_mode =
        VKMacroPorosityUpdateMode::AlgebraicSplit;
    VKMicroSolidVolumeFractionMode micro_solid_volume_fraction_mode =
        VKMicroSolidVolumeFractionMode::Reference;

    // Optional GP-local n_l initialization (future full 2C path).
    std::optional<double> initial_micro_water_content;

    // Optional Jacobian approximation for VK exchange contribution only.
    // If true, drho_L_hat/dp_L is computed by finite difference in the local
    // helper path.
    bool use_fd_jacobian_for_exchange = false;
    double fd_jacobian_perturbation = 1e-8;

    // Optional diagnostic self-check for the local implicit n_l chain rule.
    // If enabled, the code compares analytic and FD local derivatives once.
    bool check_local_jacobian = false;
    double local_jacobian_perturbation = 1e-8;
    double local_jacobian_relative_tolerance = 1e-3;

    // Optional VK-only mechanical gain on relaxation of the vdW-derived
    // compatibility pressure p_L_m = -rho_LR * mu_lR. Zero preserves the
    // current committed behavior.
    double vdw_relaxation_stress_gain = 0.0;

    // Optional VK-only mechanical gain on positive microscale water-content
    // increments. This is intended for the notebook-consistent fully saturated
    // microscale interpretation, where swelling is driven by n_l growth rather
    // than by a change in saturation.
    double micro_water_content_stress_gain = 0.0;

    // Optional VK-only notebook-style reversible swelling-strain slope driven
    // by signed microscale water-content increments Delta n_l. If enabled, the
    // VK branch uses Delta eps_sw = slope * Delta n_l and converts that
    // isotropic swelling-strain increment into a stress increment via the
    // elastic stiffness, instead of using the legacy saturation-driven
    // swelling_stress_rate path.
    double micro_water_content_swelling_slope = 0.0;
};
}  // namespace ProcessLib::RichardsMechanics
