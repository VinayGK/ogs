/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <limits>

#include "MaterialLib/SolidModels/MechanicsBase.h"

namespace ProcessLib::RichardsMechanics
{
template <int DisplacementDim>
struct PressureCoupledSolidData
{
    using SaturationStrainJacobian =
        typename MaterialLib::Solids::MechanicsBase<
            DisplacementDim>::SaturationStrainJacobian;

    bool is_active = false;
    double saturation = std::numeric_limits<double>::quiet_NaN();
    double dS_L_dp_cap = std::numeric_limits<double>::quiet_NaN();
    SaturationStrainJacobian dS_L_dStrain =
        SaturationStrainJacobian::Zero();
    MathLib::KelvinVector::KelvinVectorType<DisplacementDim>
        dSigma_dLiquidPressure =
            MathLib::KelvinVector::KelvinVectorType<DisplacementDim>::Zero();
};

template <int DisplacementDim>
PressureCoupledSolidData<DisplacementDim> makePressureCoupledSolidData(
    typename MaterialLib::Solids::MechanicsBase<
        DisplacementDim>::PressureCoupledResponse const& response)
{
    return {
        true,
        response.saturation,
        // RichardsMechanics assembles in capillary pressure p_cap = -p_L.
        -response.dSaturation_dLiquidPressure,
        response.dSaturation_dStrain,
        response.dStress_dLiquidPressure,
    };
}

}  // namespace ProcessLib::RichardsMechanics
