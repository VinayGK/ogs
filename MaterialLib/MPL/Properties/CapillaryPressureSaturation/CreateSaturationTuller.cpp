/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "BaseLib/ConfigTree.h"
#include "SaturationTuller.h"

namespace MaterialPropertyLib
{
std::unique_ptr<SaturationTuller> createSaturationTuller(
    BaseLib::ConfigTree const& config)
{
    config.checkConfigParameter("type", "SaturationTuller");

    auto property_name = config.peekConfigParameter<std::string>("name");

    DBUG("Create SaturationTuller medium property {:s}.", property_name);

    auto const area_factor_tuller =
        config.getConfigParameter<double>("area_factor_tuller");
    auto const pore_area_shape_factor_tuller =
        config.getConfigParameter<double>("pore_area_shape_factor_tuller");
    auto const characteristic_pore_size =
        config.getConfigParameter<double>("characteristic_pore_size");
    auto const surface_tension =
        config.getConfigParameter<double>("surface_tension");

    return std::make_unique<SaturationTuller>(
        std::move(property_name), area_factor_tuller,
        pore_area_shape_factor_tuller, characteristic_pore_size,
        surface_tension);
}
}  // namespace MaterialPropertyLib
