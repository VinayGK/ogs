// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <optional>

namespace ProcessLib::RichardsMechanics
{
enum class PotentialExchangeMode
{
    FullPotential
};

enum class MicroPotentialConvention
{
    PositiveReduced,
    NegativeAttractive
};

enum class LocalNonlinearSolveMode
{
    ScalarExchange,
    ScalarReferenceStorage,
    ScalarReferenceMassStorage
};

enum class MacroPorosityUpdateMode
{
    AlgebraicSplit,
    ReferenceAdditiveRate
};

enum class MicroSolidVolumeFractionMode
{
    Reference,
    CurrentPorositySplit
};

enum class PotentialExchangeRoleMapping
{
    CurrentOgs,
    MathematicaReferenceRoles
};

inline constexpr char const* toString(
    MicroPotentialConvention const convention)
{
    switch (convention)
    {
        case MicroPotentialConvention::PositiveReduced:
            return "positive_reduced";
        case MicroPotentialConvention::NegativeAttractive:
            return "negative_attractive";
    }
    return "unknown";
}

inline constexpr double microPotentialSignFactor(
    MicroPotentialConvention const convention)
{
    return convention == MicroPotentialConvention::NegativeAttractive ? -1.0
                                                                       : 1.0;
}

inline constexpr char const* toString(LocalNonlinearSolveMode const mode)
{
    switch (mode)
    {
        case LocalNonlinearSolveMode::ScalarExchange:
            return "scalar_exchange";
        case LocalNonlinearSolveMode::ScalarReferenceStorage:
            return "scalar_notebook_storage";
        case LocalNonlinearSolveMode::ScalarReferenceMassStorage:
            return "scalar_notebook_mass_storage";
    }
    return "unknown";
}

inline constexpr char const* toString(MacroPorosityUpdateMode const mode)
{
    switch (mode)
    {
        case MacroPorosityUpdateMode::AlgebraicSplit:
            return "algebraic_split";
        case MacroPorosityUpdateMode::ReferenceAdditiveRate:
            return "notebook_additive_rate";
    }
    return "unknown";
}

inline constexpr char const* toString(
    MicroSolidVolumeFractionMode const mode)
{
    switch (mode)
    {
        case MicroSolidVolumeFractionMode::Reference:
            return "reference";
        case MicroSolidVolumeFractionMode::CurrentPorositySplit:
            return "current_porosity_split";
    }
    return "unknown";
}

inline constexpr char const* toString(
    PotentialExchangeRoleMapping const mapping)
{
    switch (mapping)
    {
        case PotentialExchangeRoleMapping::CurrentOgs:
            return "current_ogs";
        case PotentialExchangeRoleMapping::MathematicaReferenceRoles:
            return "notebook_roles";
    }
    return "unknown";
}

struct PotentialExchangeParameters
{
    bool enabled = false;
    PotentialExchangeMode mode = PotentialExchangeMode::FullPotential;

    // Young-Laplace macro potential branch tolerance.
    double pressure_tolerance = 0.0;

    // vdW microscale potential parameters / reference state constants.
    double hamaker_constant = 0.0;
    double specific_surface = 0.0;
    double micro_solid_density_reference = 0.0;          // rho_SR
    double micro_solid_volume_fraction_reference = 0.0;  // n_S
    double micro_liquid_density_reference = 0.0;         // rho_l0
    double micro_liquid_density_a = 0.0;                 // a_rho
    double micro_liquid_density_b = 0.0;                 // b_rho
    MicroPotentialConvention micro_potential_convention =
        MicroPotentialConvention::PositiveReduced;
    LocalNonlinearSolveMode local_nonlinear_solve_mode =
        LocalNonlinearSolveMode::ScalarExchange;
    MacroPorosityUpdateMode macro_porosity_update_mode =
        MacroPorosityUpdateMode::AlgebraicSplit;
    MicroSolidVolumeFractionMode micro_solid_volume_fraction_mode =
        MicroSolidVolumeFractionMode::Reference;
    PotentialExchangeRoleMapping potential_role_mapping =
        PotentialExchangeRoleMapping::CurrentOgs;

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
