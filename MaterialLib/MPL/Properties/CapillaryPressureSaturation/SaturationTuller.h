// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "MaterialLib/MPL/Property.h"

namespace MaterialPropertyLib
{
class Medium;
class Phase;
class Component;

/**
 * \brief Tuller/Young-Laplace-inspired saturation model used in the VK notebook
 * macro retention law.
 *
 * The implemented saturation relation is (for capillary pressure p_c):
 * \f[
 * S_L(p_c)=
 * \begin{cases}
 * S_{L,\max}, & p_c \le p_\mathrm{tol},\\
 * S_{L,\mathrm{res}} +
 * (S_{L,\max}-S_{L,\mathrm{res}})
 * \left(1-\exp\left(-\frac{C_T}{p_c^2}\right)\right), & p_c > p_\mathrm{tol},
 * \end{cases}
 * \f]
 * with
 * \f[
 * C_T = \frac{4 F_\gamma \sigma^2}{A_n L^2}.
 * \f]
 *
 * The parameters \f$A_n\f$, \f$F_\gamma\f$, \f$L\f$, \f$\sigma\f$ correspond
 * to the VK notebook names `AreaFactorTuller`, `PoreAreaShapefactorTuller`,
 * `CharacteristicPoreSize`, and `SurfaceTension`.
 */
class SaturationTuller final : public Property
{
public:
    SaturationTuller(std::string name, double residual_liquid_saturation,
                     double maximum_liquid_saturation,
                     double area_factor_tuller,
                     double pore_area_shapefactor_tuller,
                     double characteristic_pore_size, double surface_tension,
                     double pressure_tolerance);

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
                           double const t, double const dt) const override;
    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const variable,
                            ParameterLib::SpatialPosition const& pos,
                            double const t, double const dt) const override;
    PropertyDataType d2Value(VariableArray const& variable_array,
                             Variable const variable1, Variable const variable2,
                             ParameterLib::SpatialPosition const& pos,
                             double const t, double const dt) const override;

private:
    double const S_L_res_;
    double const S_L_max_;
    double const coefficient_;
    double const pressure_tolerance_;
};
}  // namespace MaterialPropertyLib

