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

#include <memory>
#include <optional>

#include "ConstitutiveRelations/Base.h"
#include "ConstitutiveRelations/PressureCoupledSolidData.h"
#include "MaterialLib/SolidModels/MechanicsBase.h"
#include "MathLib/KelvinVector.h"
#include "ProcessLib/ConstitutiveRelations/EffectiveStressData.h"
#include "ProcessLib/ConstitutiveRelations/MechanicalStrainData.h"

namespace ProcessLib
{
namespace RichardsMechanics
{
template <typename BMatricesType, typename ShapeMatrixTypeDisplacement,
          typename ShapeMatricesTypePressure, int DisplacementDim, int NPoints>
struct IntegrationPointData final
{
    typename ShapeMatrixTypeDisplacement::NodalRowVectorType N_u;
    typename ShapeMatrixTypeDisplacement::GlobalDimNodalMatrixType dNdx_u;

    struct ConstitutiveRelationUpdate final
    {
        typename BMatricesType::KelvinMatrixType C;
        std::optional<PressureCoupledSolidData<DisplacementDim>>
            pressure_coupled_data;
    };

    typename ShapeMatricesTypePressure::NodalRowVectorType N_p;
    typename ShapeMatricesTypePressure::GlobalDimNodalMatrixType dNdx_p;

    double integration_weight = std::numeric_limits<double>::quiet_NaN();

    MathLib::KelvinVector::
        KelvinMatrixType<DisplacementDim> static computeElasticTangentStiffness(
            MaterialPropertyLib::VariableArray const& variable_array,
            double const t,
            ParameterLib::SpatialPosition const& x_position,
            double const dt,
            MaterialLib::Solids::MechanicsBase<DisplacementDim> const&
                solid_material,
            typename MaterialLib::Solids::MechanicsBase<DisplacementDim>::
                MaterialStateVariables const& material_state_variables)
    {
        namespace MPL = MaterialPropertyLib;

        MPL::VariableArray variable_array_prev = variable_array;

        auto&& solution = solid_material.integrateStress(
            variable_array_prev, variable_array, t, x_position, dt,
            material_state_variables);

        if (!solution)
        {
            OGS_FATAL("Computation of elastic tangent stiffness failed.");
        }

        MathLib::KelvinVector::KelvinMatrixType<DisplacementDim> C =
            std::move(std::get<2>(*solution));

        return C;
    }

    static ConstitutiveRelationUpdate updateConstitutiveRelation(
        MaterialPropertyLib::VariableArray const& variable_array,
        MaterialPropertyLib::VariableArray const& variable_array_prev,
        double const t,
        ParameterLib::SpatialPosition const& x_position,
        double const dt,
        double const temperature,
        ProcessLib::ConstitutiveRelations::EffectiveStressData<DisplacementDim>&
            sigma_eff,
        PrevState<ProcessLib::ConstitutiveRelations::EffectiveStressData<
            DisplacementDim>> const& sigma_eff_prev,
        ProcessLib::ConstitutiveRelations::MechanicalStrainData<
            DisplacementDim> const&
        /*eps_m*/,
        PrevState<ProcessLib::ConstitutiveRelations::MechanicalStrainData<
            DisplacementDim>> const& eps_m_prev,
        MaterialLib::Solids::MechanicsBase<DisplacementDim> const&
            solid_material,
        std::unique_ptr<typename MaterialLib::Solids::MechanicsBase<
            DisplacementDim>::MaterialStateVariables>& material_state_variables)
    {
        auto variable_array_prev_local = variable_array_prev;
        variable_array_prev_local.stress = sigma_eff_prev->sigma_eff;
        variable_array_prev_local.mechanical_strain
            .emplace<MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
                eps_m_prev->eps_m);
        variable_array_prev_local.temperature = temperature;

        if (auto pressure_coupled_response =
                solid_material.integrateStressPressureCoupled(
                    variable_array_prev_local, variable_array, t, x_position,
                    dt, *material_state_variables))
        {
            auto pressure_coupled_data =
                makePressureCoupledSolidData<DisplacementDim>(
                    *pressure_coupled_response);
            sigma_eff.sigma_eff = pressure_coupled_response->stress;
            material_state_variables =
                std::move(pressure_coupled_response->state);
            return {std::move(pressure_coupled_response->dStress_dStrain),
                    std::move(pressure_coupled_data)};
        }

        auto&& solution = solid_material.integrateStress(
            variable_array_prev_local, variable_array, t, x_position, dt,
            *material_state_variables);

        if (!solution)
        {
            OGS_FATAL("Computation of local constitutive relation failed.");
        }

        typename BMatricesType::KelvinMatrixType C;
        std::tie(sigma_eff.sigma_eff, material_state_variables, C) =
            std::move(*solution);

        return {std::move(C), std::nullopt};
    }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

}  // namespace RichardsMechanics
}  // namespace ProcessLib
