/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "CreateMFrontRichardsMechanics.h"

#include "CreateMFrontGeneric.h"
#include "MFrontRichardsMechanics.h"

namespace MaterialLib::Solids::MFront
{
template <int DisplacementDim>
std::unique_ptr<MechanicsBase<DisplacementDim>> createMFrontRichardsMechanics(
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    BaseLib::ConfigTree const& config)
{
    auto conf = createMFrontConfig(DisplacementDim, parameters, config);

    return std::make_unique<MFrontRichardsMechanics<DisplacementDim>>(
        std::move(conf.behaviour), std::move(conf.material_properties),
        std::move(conf.state_variables_initial_properties),
        local_coordinate_system);
}
}  // namespace MaterialLib::Solids::MFront

namespace MaterialLib::Solids::MFront
{
template std::unique_ptr<MechanicsBase<2>> createMFrontRichardsMechanics<2>(
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    BaseLib::ConfigTree const& config);

template std::unique_ptr<MechanicsBase<3>> createMFrontRichardsMechanics<3>(
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    BaseLib::ConfigTree const& config);
}  // namespace MaterialLib::Solids::MFront
