/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <gtest/gtest.h>

#include "MaterialLib/SolidModels/MechanicsBase.h"
#include "ProcessLib/RichardsMechanics/ConstitutiveRelations/PressureCoupledSolidData.h"

TEST(RichardsMechanics, makePressureCoupledSolidData)
{
    constexpr int DisplacementDim = 2;
    using MechanicsBase =
        MaterialLib::Solids::MechanicsBase<DisplacementDim>;
    using KelvinMatrix = typename MechanicsBase::KelvinMatrix;
    using KelvinVector = typename MechanicsBase::KelvinVector;
    using SaturationStrainJacobian =
        typename MechanicsBase::SaturationStrainJacobian;

    typename MechanicsBase::PressureCoupledResponse response{
        KelvinVector::Zero(),
        0.625,
        KelvinMatrix::Identity(),
        KelvinVector::Constant(2.5),
        SaturationStrainJacobian::Constant(1.5),
        4.0,
        std::make_unique<typename MechanicsBase::MaterialStateVariables>()};

    auto const mapped =
        ProcessLib::RichardsMechanics::makePressureCoupledSolidData<
            DisplacementDim>(response);

    EXPECT_TRUE(mapped.is_active);
    EXPECT_DOUBLE_EQ(mapped.saturation, 0.625);
    EXPECT_DOUBLE_EQ(mapped.dS_L_dp_cap, -4.0);
    EXPECT_TRUE(
        mapped.dS_L_dStrain.isApprox(SaturationStrainJacobian::Constant(1.5)));
    EXPECT_TRUE(mapped.dSigma_dLiquidPressure.isApprox(
        KelvinVector::Constant(2.5)));
}
