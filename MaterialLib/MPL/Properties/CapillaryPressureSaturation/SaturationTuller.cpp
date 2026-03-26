/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "SaturationTuller.h"

#include <cmath>

namespace MaterialPropertyLib
{
SaturationTuller::SaturationTuller(std::string name,
                                   double const area_factor_tuller,
                                   double const pore_area_shape_factor_tuller,
                                   double const characteristic_pore_size,
                                   double const surface_tension)
    : area_factor_tuller_(area_factor_tuller),
      pore_area_shape_factor_tuller_(pore_area_shape_factor_tuller),
      characteristic_pore_size_(characteristic_pore_size),
      surface_tension_(surface_tension),
      capillary_prefactor_(4.0 * pore_area_shape_factor_tuller_ *
                           surface_tension_ * surface_tension_ /
                           (area_factor_tuller_ * characteristic_pore_size_ *
                            characteristic_pore_size_))
{
    name_ = std::move(name);

    if (area_factor_tuller_ <= 0.0)
    {
        OGS_FATAL("SaturationTuller requires area_factor_tuller > 0.");
    }
    if (pore_area_shape_factor_tuller_ <= 0.0)
    {
        OGS_FATAL(
            "SaturationTuller requires pore_area_shape_factor_tuller > 0.");
    }
    if (characteristic_pore_size_ <= 0.0)
    {
        OGS_FATAL("SaturationTuller requires characteristic_pore_size > 0.");
    }
    if (surface_tension_ <= 0.0)
    {
        OGS_FATAL("SaturationTuller requires surface_tension > 0.");
    }
}

PropertyDataType SaturationTuller::value(VariableArray const& variable_array,
                                         ParameterLib::SpatialPosition const&,
                                         double const,
                                         double const) const
{
    double const p_cap = variable_array.capillary_pressure;

    if (p_cap <= pressure_tolerance_)
    {
        return 1.0;
    }

    auto const capillary_term = capillary_prefactor_ / (p_cap * p_cap);
    return 1.0 - std::exp(-capillary_term);
}

PropertyDataType SaturationTuller::dValue(
    VariableArray const& variable_array,
    Variable const variable,
    ParameterLib::SpatialPosition const&,
    double const,
    double const) const
{
    if (variable != Variable::capillary_pressure)
    {
        OGS_FATAL(
            "SaturationTuller::dValue is implemented for derivatives with "
            "respect to capillary pressure only.");
    }

    double const p_cap = variable_array.capillary_pressure;

    if (p_cap <= pressure_tolerance_)
    {
        return 0.0;
    }

    auto const inv_p_cap = 1.0 / p_cap;
    auto const exp_term =
        std::exp(-capillary_prefactor_ * inv_p_cap * inv_p_cap);
    return -2.0 * capillary_prefactor_ * exp_term * inv_p_cap * inv_p_cap *
           inv_p_cap;
}

PropertyDataType SaturationTuller::d2Value(
    VariableArray const& variable_array,
    Variable const variable1,
    Variable const variable2,
    ParameterLib::SpatialPosition const&,
    double const,
    double const) const
{
    if (variable1 != Variable::capillary_pressure ||
        variable2 != Variable::capillary_pressure)
    {
        OGS_FATAL(
            "SaturationTuller::d2Value is implemented for second derivatives "
            "with respect to capillary pressure only.");
    }

    double const p_cap = variable_array.capillary_pressure;

    if (p_cap <= pressure_tolerance_)
    {
        return 0.0;
    }

    auto const p_cap2 = p_cap * p_cap;
    auto const p_cap4 = p_cap2 * p_cap2;
    auto const p_cap6 = p_cap4 * p_cap2;
    auto const exp_term = std::exp(-capillary_prefactor_ / p_cap2);

    return exp_term *
           (6.0 * capillary_prefactor_ / p_cap4 -
            4.0 * capillary_prefactor_ * capillary_prefactor_ / p_cap6);
}
}  // namespace MaterialPropertyLib
