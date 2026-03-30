/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <algorithm>

#include "MFrontGeneric.h"
#include "Variable.h"

namespace MaterialLib
{
namespace Solids
{
namespace MFront
{
template <int DisplacementDim>
class MFrontRichardsMechanics
    : private MFrontGeneric<DisplacementDim,
                            boost::mp11::mp_list<Strain, LiquidPressure>,
                            boost::mp11::mp_list<Stress, Saturation>,
                            boost::mp11::mp_list<Temperature>>,
      public MechanicsBase<DisplacementDim>
{
    using Base = MFrontGeneric<DisplacementDim,
                               boost::mp11::mp_list<Strain, LiquidPressure>,
                               boost::mp11::mp_list<Stress, Saturation>,
                               boost::mp11::mp_list<Temperature>>;
    using KelvinVector = typename Base::KelvinVector;
    using KelvinMatrix = typename Base::KelvinMatrix;

public:
    using PressureCoupledResponse =
        typename MechanicsBase<DisplacementDim>::PressureCoupledResponse;

    using Base::Base;

    std::unique_ptr<
        typename MechanicsBase<DisplacementDim>::MaterialStateVariables>
    createMaterialStateVariables() const override
    {
        return Base::createMaterialStateVariables();
    }

    void initializeInternalStateVariables(
        double const t,
        ParameterLib::SpatialPosition const& x,
        typename MechanicsBase<DisplacementDim>::MaterialStateVariables&
            material_state_variables) const override
    {
        Base::initializeInternalStateVariables(t, x, material_state_variables);
    }

    std::optional<PressureCoupledResponse> integrateStressPressureCoupled(
        MaterialPropertyLib::VariableArray const& variable_array_prev,
        MaterialPropertyLib::VariableArray const& variable_array,
        double const t,
        ParameterLib::SpatialPosition const& x,
        double const dt,
        typename MechanicsBase<DisplacementDim>::MaterialStateVariables const&
            material_state_variables) const override
    {
        auto res = Base::integrateStress(variable_array_prev,
                                         variable_array,
                                         t,
                                         x,
                                         dt,
                                         material_state_variables);

        if (!res)
        {
            return std::nullopt;
        }

        auto& [forces_data, state, tangent_operator_data] = *res;
        auto const view = this->createThermodynamicForcesView();
        auto const swelling_stress =
            getInternalVectorByName(*state, "sigma_S");
        auto const liquid_mass_exchange_source =
            -getInternalScalarByName(*state, "rho_l_hat");

        PressureCoupledResponse response{
            view.block(stress, forces_data),
            view.block(saturation, forces_data),
            blocks_view_.block(stress, strain, tangent_operator_data),
            blocks_view_.block(stress, liquid_pressure, tangent_operator_data),
            blocks_view_.block(saturation, strain, tangent_operator_data),
            blocks_view_.block(saturation, liquid_pressure,
                               tangent_operator_data),
            swelling_stress,
            liquid_mass_exchange_source,
            std::move(state)};
        return response;
    }

    std::optional<std::tuple<KelvinVector,
                             std::unique_ptr<typename MechanicsBase<
                                 DisplacementDim>::MaterialStateVariables>,
                             KelvinMatrix>>
    integrateStress(
        MaterialPropertyLib::VariableArray const& variable_array_prev,
        MaterialPropertyLib::VariableArray const& variable_array,
        double const t,
        ParameterLib::SpatialPosition const& x,
        double const dt,
        typename MechanicsBase<DisplacementDim>::MaterialStateVariables const&
            material_state_variables) const override
    {
        auto response = integrateStressPressureCoupled(variable_array_prev,
                                                       variable_array,
                                                       t,
                                                       x,
                                                       dt,
                                                       material_state_variables);

        if (!response)
        {
            return std::nullopt;
        }

        return std::optional<
            std::tuple<KelvinVector,
                       std::unique_ptr<typename MechanicsBase<
                           DisplacementDim>::MaterialStateVariables>,
                       KelvinMatrix>>{std::in_place,
                                      response->stress,
                                      std::move(response->state),
                                      response->dStress_dStrain};
    }

    std::vector<typename MechanicsBase<DisplacementDim>::InternalVariable>
    getInternalVariables() const override
    {
        return Base::getInternalVariables();
    }

    double getBulkModulus(double const t,
                          ParameterLib::SpatialPosition const& x,
                          KelvinMatrix const* const C) const override
    {
        return Base::getBulkModulus(t, x, C);
    }

    double computeFreeEnergyDensity(
        double const t,
        ParameterLib::SpatialPosition const& x,
        double const dt,
        KelvinVector const& eps,
        KelvinVector const& sigma,
        typename MechanicsBase<DisplacementDim>::MaterialStateVariables const&
            material_state_variables) const override
    {
        return Base::computeFreeEnergyDensity(
            t, x, dt, eps, sigma, material_state_variables);
    }

private:
    KelvinVector getInternalVectorByName(
        typename MechanicsBase<DisplacementDim>::MaterialStateVariables const&
            material_state_variables,
        std::string const& name) const
    {
        auto const internal_variables = Base::getInternalVariables();
        auto const it =
            std::find_if(internal_variables.begin(), internal_variables.end(),
                         [&name](auto const& internal_variable)
                         { return internal_variable.name == name; });

        KelvinVector values = KelvinVector::Zero();
        if (it == internal_variables.end() ||
            it->num_components != values.rows())
        {
            return values;
        }

        std::vector<double> cache;
        auto const& raw_values = it->getter(material_state_variables, cache);
        for (Eigen::Index i = 0;
             i < values.rows() &&
             i < static_cast<Eigen::Index>(raw_values.size());
             ++i)
        {
            values[i] = raw_values[static_cast<std::size_t>(i)];
        }
        return values;
    }

    double getInternalScalarByName(
        typename MechanicsBase<DisplacementDim>::MaterialStateVariables const&
            material_state_variables,
        std::string const& name) const
    {
        auto const internal_variables = Base::getInternalVariables();
        auto const it =
            std::find_if(internal_variables.begin(), internal_variables.end(),
                         [&name](auto const& internal_variable)
                         { return internal_variable.name == name; });

        if (it == internal_variables.end() || it->num_components != 1)
        {
            return 0.0;
        }

        std::vector<double> cache;
        auto const& values = it->getter(material_state_variables, cache);
        return values.empty() ? 0.0 : values[0];
    }

    OGSMFrontTangentOperatorBlocksView<
        DisplacementDim,
        ForcesGradsCombinations<boost::mp11::mp_list<Strain, LiquidPressure>,
                                boost::mp11::mp_list<Stress, Saturation>,
                                boost::mp11::mp_list<Temperature>>::type>
        blocks_view_ = this->createTangentOperatorBlocksView();
};

extern template class MFrontRichardsMechanics<2>;
extern template class MFrontRichardsMechanics<3>;

}  // namespace MFront
}  // namespace Solids
}  // namespace MaterialLib
