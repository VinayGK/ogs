/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#pragma once

#include "MaterialLib/MPL/Property.h"

namespace MaterialPropertyLib
{
/**
 * \brief Tuller-style saturation law used by DSM micro-macro coupling.
 *
 * \details This property must be a medium property. It computes the liquid
 * saturation from capillary pressure using the same scalar law as the
 * DSM micro-scale unsaturated branch.
 *
 * For capillary pressure \f$p_\mathrm{cap}\f$ the saturation is
 *
 * \f[
 * S_L(p_\mathrm{cap}) =
 * \begin{cases}
 * 1 - \exp\!\left(-\dfrac{A_T}{p_\mathrm{cap}^2}\right),
 *     & p_\mathrm{cap} > p_\mathrm{tol},\\[1ex]
 * 1, & p_\mathrm{cap} \le p_\mathrm{tol},
 * \end{cases}
 * \f]
 *
 * where
 *
 * \f[
 * A_T =
 * \frac{4 \beta_T \gamma^2}{a_T r_c^2}.
 * \f]
 */
class SaturationTuller final : public Property
{
public:
    SaturationTuller(std::string name,
                     double area_factor_tuller,
                     double pore_area_shape_factor_tuller,
                     double characteristic_pore_size,
                     double surface_tension);

    void checkScale() const override
    {
        if (!std::holds_alternative<Medium*>(scale_))
        {
            OGS_FATAL(
                "The property 'SaturationTuller' is implemented on the "
                "'media' scale only.");
        }
    }

    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& pos,
                           double const t,
                           double const dt) const override;

    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const variable,
                            ParameterLib::SpatialPosition const& pos,
                            double const t,
                            double const dt) const override;

    PropertyDataType d2Value(VariableArray const& variable_array,
                             Variable const variable1,
                             Variable const variable2,
                             ParameterLib::SpatialPosition const& pos,
                             double const t,
                             double const dt) const override;

private:
    double const area_factor_tuller_;
    double const pore_area_shape_factor_tuller_;
    double const characteristic_pore_size_;
    double const surface_tension_;
    double const capillary_prefactor_;
    static constexpr double pressure_tolerance_ = 1e-12;
};
}  // namespace MaterialPropertyLib
