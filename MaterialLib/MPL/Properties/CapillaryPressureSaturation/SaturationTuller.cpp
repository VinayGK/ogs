// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "SaturationTuller.h"

#include <algorithm>
#include <cassert>
#include <cmath>

namespace MaterialPropertyLib
{
namespace
{
// Option B film coefficient C_film [Pa^(1/3)].
// From the vdW disjoining-pressure balance p_c = A_H/(6 pi h^3):
//   h(p_c)   = (A_H/(6 pi p_c))^(1/3)                                  [m]
//   S_film   = (a_v/phi_M) * h(p_c) = C_film * p_c^(-1/3)              [-]
//   C_film   = (a_v/phi_M) * (A_H/(6 pi))^(1/3)
// Units: [1/m] * (J/[1])^(1/3) = [1/m] * (Pa*m^3)^(1/3)
//        = [1/m] * Pa^(1/3) * m = Pa^(1/3)  -> S_film dimensionless. ok
// The 6*pi prefactor matches the micro vdW convention in
// PotentialExchange.h (disjoining pressure Pi = A/(6 pi h^3)).
double computeFilmCoefficient(double const macro_specific_surface,
                              double const hamaker_constant,
                              double const macro_porosity)
{
    constexpr double pi = 3.141592653589793238462643383279502884;
    return (macro_specific_surface / macro_porosity) *
           std::cbrt(hamaker_constant / (6.0 * pi));  // [Pa^(1/3)]
}
}  // namespace

SaturationTuller::SaturationTuller(
    std::string name, double const residual_liquid_saturation,
    double const maximum_liquid_saturation, double const area_factor_tuller,
    double const pore_area_shapefactor_tuller,
    double const characteristic_pore_size, double const surface_tension,
    double const pressure_tolerance, double const macro_specific_surface,
    double const hamaker_constant, double const macro_porosity)
    : S_L_res_(residual_liquid_saturation),
      S_L_max_(maximum_liquid_saturation),
      coefficient_(4.0 * pore_area_shapefactor_tuller * surface_tension *
                   surface_tension /
                   (area_factor_tuller * characteristic_pore_size *
                    characteristic_pore_size)),
      pressure_tolerance_(pressure_tolerance),
      film_active_(macro_specific_surface > 0.0),
      film_coefficient_(film_active_
                            ? computeFilmCoefficient(macro_specific_surface,
                                                     hamaker_constant,
                                                     macro_porosity)
                            : 0.0)
{
    name_ = std::move(name);

    if (!(0.0 <= S_L_res_ && S_L_res_ <= S_L_max_ && S_L_max_ <= 1.0))
    {
        OGS_FATAL(
            "SaturationTuller bounds must satisfy 0 <= S_L_res <= S_L_max <= "
            "1, but got S_L_res = {:g}, S_L_max = {:g}.",
            S_L_res_, S_L_max_);
    }

    if (!(area_factor_tuller > 0.0))
    {
        OGS_FATAL("SaturationTuller requires area_factor_tuller > 0, got {:g}.",
                  area_factor_tuller);
    }
    if (!(pore_area_shapefactor_tuller > 0.0))
    {
        OGS_FATAL(
            "SaturationTuller requires pore_area_shapefactor_tuller > 0, got "
            "{:g}.",
            pore_area_shapefactor_tuller);
    }
    if (!(characteristic_pore_size > 0.0))
    {
        OGS_FATAL(
            "SaturationTuller requires characteristic_pore_size > 0, got "
            "{:g}.",
            characteristic_pore_size);
    }
    if (!(surface_tension > 0.0))
    {
        OGS_FATAL("SaturationTuller requires surface_tension > 0, got {:g}.",
                  surface_tension);
    }
    if (!(pressure_tolerance_ >= 0.0))
    {
        OGS_FATAL(
            "SaturationTuller requires pressure_tolerance >= 0, got {:g}.",
            pressure_tolerance_);
    }
    if (!(coefficient_ > 0.0))
    {
        OGS_FATAL("SaturationTuller internal coefficient must be > 0, got {:g}.",
                  coefficient_);
    }

    // Option B film branch: when the macro specific surface is supplied
    // (film_active_), the Hamaker constant and macro porosity must be valid.
    if (film_active_)
    {
        if (!(hamaker_constant > 0.0))
        {
            OGS_FATAL(
                "SaturationTuller film branch requires hamaker_constant > 0 "
                "when macro_specific_surface > 0, got {:g}.",
                hamaker_constant);
        }
        if (!(macro_porosity > 0.0 && macro_porosity <= 1.0))
        {
            OGS_FATAL(
                "SaturationTuller film branch requires 0 < macro_porosity <= "
                "1 when macro_specific_surface > 0, got {:g}.",
                macro_porosity);
        }
    }
}

PropertyDataType SaturationTuller::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    double const p_cap = variable_array.capillary_pressure;
    if (p_cap <= pressure_tolerance_)
    {
        return S_L_max_;
    }

    double const e = std::exp(-coefficient_ / (p_cap * p_cap));
    double const S_eff = 1.0 - e;
    double S = S_L_res_ + (S_L_max_ - S_L_res_) * S_eff;

    if (film_active_)
    {
        // Adsorptive-film contribution S_film = C_film * p_cap^(-1/3) [-],
        // coexisting with the capillary corner at every p_cap (no junction).
        S += film_coefficient_ * std::cbrt(1.0 / p_cap);  // [-]
    }

    return std::clamp(S, S_L_res_, S_L_max_);
}

PropertyDataType SaturationTuller::dValue(
    VariableArray const& variable_array, Variable const variable,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
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

    double const e = std::exp(-coefficient_ / (p_cap * p_cap));
    double const dS_eff_dp_cap = -(2.0 * coefficient_ / (p_cap * p_cap * p_cap)) * e;
    double dS = dS_eff_dp_cap * (S_L_max_ - S_L_res_);

    if (film_active_)
    {
        // d/dp_cap [C_film * p_cap^(-1/3)] = -(1/3) C_film * p_cap^(-4/3) [1/Pa]
        dS += -(film_coefficient_ / 3.0) * std::pow(p_cap, -4.0 / 3.0);  // [1/Pa]
    }
    return dS;
}

PropertyDataType SaturationTuller::d2Value(
    VariableArray const& variable_array, Variable const variable1,
    Variable const variable2, ParameterLib::SpatialPosition const& /*pos*/,
    double const /*t*/, double const /*dt*/) const
{
    (void)variable1;
    (void)variable2;
    assert((variable1 == Variable::capillary_pressure) &&
           (variable2 == Variable::capillary_pressure) &&
           "SaturationTuller::d2Value is implemented for derivatives with "
           "respect to capillary pressure only.");

    double const p_cap = variable_array.capillary_pressure;
    if (p_cap <= pressure_tolerance_)
    {
        return 0.0;
    }

    double const p2 = p_cap * p_cap;
    double const p4 = p2 * p2;
    double const p6 = p4 * p2;
    double const e = std::exp(-coefficient_ / p2);
    double const d2S_eff_dp_cap2 =
        2.0 * coefficient_ * e * (3.0 * p2 - 2.0 * coefficient_) / p6;
    double d2S = d2S_eff_dp_cap2 * (S_L_max_ - S_L_res_);

    if (film_active_)
    {
        // d2/dp_cap^2 [C_film * p_cap^(-1/3)] = (4/9) C_film * p_cap^(-7/3) [1/Pa^2]
        d2S += (4.0 / 9.0) * film_coefficient_ * std::pow(p_cap, -7.0 / 3.0);  // [1/Pa^2]
    }
    return d2S;
}
}  // namespace MaterialPropertyLib

