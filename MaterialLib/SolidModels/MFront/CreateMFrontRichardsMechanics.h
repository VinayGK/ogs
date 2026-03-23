/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <memory>
#include <vector>

#include "BaseLib/ConfigTree-fwd.h"
#include "MaterialLib/SolidModels/MechanicsBase.h"
#include "ParameterLib/Parameter.h"

namespace MaterialLib::Solids::MFront
{
template <int DisplacementDim>
std::unique_ptr<MechanicsBase<DisplacementDim>> createMFrontRichardsMechanics(
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    BaseLib::ConfigTree const& config);

extern template std::unique_ptr<MechanicsBase<2>>
createMFrontRichardsMechanics<2>(
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    BaseLib::ConfigTree const& config);

extern template std::unique_ptr<MechanicsBase<3>>
createMFrontRichardsMechanics<3>(
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    BaseLib::ConfigTree const& config);
}  // namespace MaterialLib::Solids::MFront
