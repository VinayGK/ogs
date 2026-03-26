/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#include <gtest/gtest.h>

#include <cmath>

#include "MaterialLib/MPL/Properties/CapillaryPressureSaturation/SaturationTuller.h"

namespace MPL = MaterialPropertyLib;

TEST(MaterialPropertyLib, SaturationTuller)
{
    double constexpr area_factor_tuller = 1.0;
    double constexpr pore_area_shape_factor_tuller = 0.8584073464102069;
    double constexpr characteristic_pore_size = 1e-5;
    double constexpr surface_tension = 0.0715;

    MPL::Property const& saturation = MPL::SaturationTuller{
        "saturation", area_factor_tuller, pore_area_shape_factor_tuller,
        characteristic_pore_size, surface_tension};

    MPL::VariableArray variable_array;
    ParameterLib::SpatialPosition const pos;
    double const t = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();

    double const capillary_prefactor =
        4.0 * pore_area_shape_factor_tuller * surface_tension *
        surface_tension /
        (area_factor_tuller * characteristic_pore_size *
         characteristic_pore_size);

    variable_array.capillary_pressure = 0.0;
    EXPECT_DOUBLE_EQ(saturation.template value<double>(variable_array, pos, t, dt),
                     1.0);
    EXPECT_DOUBLE_EQ(
        saturation.template dValue<double>(
            variable_array, MPL::Variable::capillary_pressure, pos, t, dt),
        0.0);
    EXPECT_DOUBLE_EQ(
        saturation.template d2Value<double>(
            variable_array, MPL::Variable::capillary_pressure,
            MPL::Variable::capillary_pressure, pos, t, dt),
        0.0);

    for (double p_cap : {1e5, 2.5e5, 5e5, 7.5e5, 1e6})
    {
        variable_array.capillary_pressure = p_cap;

        auto const expected_saturation =
            1.0 - std::exp(-capillary_prefactor / (p_cap * p_cap));
        auto const expected_dS =
            -2.0 * capillary_prefactor *
            std::exp(-capillary_prefactor / (p_cap * p_cap)) /
            (p_cap * p_cap * p_cap);
        auto const expected_d2S =
            std::exp(-capillary_prefactor / (p_cap * p_cap)) *
            (6.0 * capillary_prefactor / std::pow(p_cap, 4) -
             4.0 * capillary_prefactor * capillary_prefactor /
                 std::pow(p_cap, 6));

        EXPECT_NEAR(saturation.template value<double>(variable_array, pos, t, dt),
                    expected_saturation, 1e-15);
        EXPECT_NEAR(
            saturation.template dValue<double>(
                variable_array, MPL::Variable::capillary_pressure, pos, t, dt),
            expected_dS, 1e-20);
        EXPECT_NEAR(
            saturation.template d2Value<double>(
                variable_array, MPL::Variable::capillary_pressure,
                MPL::Variable::capillary_pressure, pos, t, dt),
            expected_d2S, 1e-25);
    }
}
