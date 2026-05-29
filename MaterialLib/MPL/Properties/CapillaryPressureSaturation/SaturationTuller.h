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
 * \brief Tuller/Young-Laplace-inspired saturation model used in the DSM
 *        micro-macro reference branch
 * macro retention law.
 *
 * The implemented saturation relation is (for capillary pressure p_c):
 * \f[
 * S_L(p_c)=
 * \begin{cases}
 * S_{L,\max}, & p_c \le p_\mathrm{tol},\\
 * S_{L,\mathrm{res}} +
 * (S_{L,\max}-S_{L,\mathrm{res}})
 * \left(1-\exp\left(-\frac{C_T}{p_c^2}\right)\right)
 * + S_\mathrm{film}(p_c), & p_c > p_\mathrm{tol},
 * \end{cases}
 * \f]
 * with
 * \f[
 * C_T = \frac{4 F_\gamma \sigma^2}{A_n L^2}.
 * \f]
 *
 * The parameters \f$A_n\f$, \f$F_\gamma\f$, \f$L\f$, \f$\sigma\f$ correspond
 * to the DSM parameter names `AreaFactorTuller`, `PoreAreaShapefactorTuller`,
 * `CharacteristicPoreSize`, and `SurfaceTension`.
 *
 * \par Option B film branch (coexistence closure)
 * When the macro specific surface \f$a_v\f$ is supplied (and positive) an
 * adsorptive-film contribution is added at every capillary pressure, with
 * no air-entry threshold and no junction/blend. The film thickness follows
 * the van der Waals disjoining-pressure balance
 * \f$p_c = A_H/(6\pi h^3)\f$ (same prefactor and physics as the micro vdW
 * potential in PotentialExchange.h), giving
 * \f[
 * h(p_c) = \left(\frac{A_H}{6\pi p_c}\right)^{1/3},\qquad
 * S_\mathrm{film}(p_c) = \frac{a_v}{\phi^\mathrm{Macro}}\,h(p_c)
 *        = C_\mathrm{film}\,p_c^{-1/3},\quad
 * C_\mathrm{film} = \frac{a_v}{\phi^\mathrm{Macro}}
 *        \left(\frac{A_H}{6\pi}\right)^{1/3}.
 * \f]
 * The closure is parameterized through water content (\f$S_\mathrm{film}\f$ is
 * a saturation); \f$h\f$ is an internal pointwise intermediate only, never a
 * state variable or calibration target. Setting \f$a_v = 0\f$ (the default)
 * disables the film branch and recovers the pure-Tuller relation exactly.
 */
class SaturationTuller final : public Property
{
public:
    SaturationTuller(std::string name, double residual_liquid_saturation,
                     double maximum_liquid_saturation,
                     double area_factor_tuller,
                     double pore_area_shapefactor_tuller,
                     double characteristic_pore_size, double surface_tension,
                     double pressure_tolerance,
                     double macro_specific_surface = 0.0,
                     double hamaker_constant = 0.0,
                     double macro_porosity = 0.0);

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

    // Option B film branch. film_active_ gates the adsorptive-film term;
    // when false the relation reduces to pure Tuller exactly.
    bool const film_active_;
    double const film_coefficient_;  // C_film [Pa^(1/3)], = 0 when inactive
};
}  // namespace MaterialPropertyLib
