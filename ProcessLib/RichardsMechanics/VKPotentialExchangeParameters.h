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

    // Optional GP-local n_l initialization (future full 2C path).
    std::optional<double> initial_micro_water_content;

    // Optional Jacobian approximation for VK exchange contribution only.
    // If true, drho_L_hat/dp_L is computed by finite difference in the local
    // helper path.
    bool use_fd_jacobian_for_exchange = false;
    double fd_jacobian_perturbation = 1e-8;
};
}  // namespace ProcessLib::RichardsMechanics
