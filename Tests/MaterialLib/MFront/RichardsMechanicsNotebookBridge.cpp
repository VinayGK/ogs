/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#ifdef OGS_USE_MFRONT

#include <algorithm>
#include <cmath>
#include <memory>
#include <optional>
#include <vector>

#include <gtest/gtest.h>
#include <boost/mp11.hpp>

#include "BaseLib/ConfigTree.h"
#include "MaterialLib/SolidModels/MFront/CreateMFrontGeneric.h"
#include "MaterialLib/SolidModels/MFront/Variable.h"
#include "ParameterLib/ConstantParameter.h"
#include "Tests/TestTools.h"

namespace MSM = MaterialLib::Solids::MFront;
namespace MPL = MaterialPropertyLib;
using KV = MathLib::KelvinVector::KelvinVectorType<3>;

static auto createParameters()
{
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> parameters;

    auto add_param = [&parameters](char const* name, double const value)
    {
        parameters.push_back(
            std::make_unique<ParameterLib::ConstantParameter<double>>(name,
                                                                      value));
    };

    add_param("E", 1e9);
    add_param("nu", 0.3);
    add_param("alpha_b", 1.0);
    add_param("p_sw", 2.0);
    add_param("lambda_v", 0.5);
    add_param("pc0", 2e5);
    add_param("v0", 1.7857142857142857);
    add_param("p_sat_scale", 1e3);

    return parameters;
}

static auto createBridgeModel(
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters)
{
    auto local_coordinate_system = std::nullopt;

    const char* xml = R"XML(
        <type>MFrontRichardsMechanics</type>
        <behaviour>RichardsMechanicsNotebookBridge</behaviour>
        <library path_is_relative_to_prj_file="false">libOgsMFrontBehaviour</library>
        <material_properties>
            <material_property name="YoungModulus" parameter="E"/>
            <material_property name="PoissonRatio" parameter="nu"/>
            <material_property name="CriticalStateLineSlope" parameter="alpha_b"/>
            <material_property name="SwellingLineSlope" parameter="p_sw"/>
            <material_property name="VirginConsolidationLineSlope" parameter="lambda_v"/>
            <material_property name="CharacteristicPreConsolidationPressure" parameter="pc0"/>
            <material_property name="InitialVolumeRatio" parameter="v0"/>
            <material_property name="SaturationPressureScale" parameter="p_sat_scale"/>
        </material_properties>
        <initial_values>
            <state_variable name="PreConsolidationPressure" parameter="pc0"/>
            <state_variable name="VolumeRatio" parameter="v0"/>
        </initial_values>
        )XML";

    auto ptree = Tests::readXml(xml);
    BaseLib::ConfigTree config_tree(std::move(ptree), "FILENAME",
                                    &BaseLib::ConfigTree::onerror,
                                    &BaseLib::ConfigTree::onwarning);

    return MSM::createMFrontGeneric<
        3, boost::mp11::mp_list<MSM::Strain, MSM::LiquidPressure>,
        boost::mp11::mp_list<MSM::Stress, MSM::Saturation>,
        boost::mp11::mp_list<MSM::Temperature>>(
        parameters, local_coordinate_system, config_tree);
}

TEST(MaterialLib_RichardsMechanicsNotebookBridgeMFront,
     OneStepPressureResponse)
{
    auto const parameters = createParameters();
    auto model = createBridgeModel(parameters);
    ASSERT_TRUE(model != nullptr);

    auto state = model->createMaterialStateVariables();
    ASSERT_TRUE(state != nullptr);

    MPL::VariableArray variable_array_prev;
    variable_array_prev.stress.template emplace<KV>(KV::Zero());
    variable_array_prev.mechanical_strain.template emplace<KV>(KV::Zero());
    variable_array_prev.liquid_phase_pressure = 0.0;
    variable_array_prev.liquid_saturation = 0.0;
    variable_array_prev.temperature = 0.0;

    MPL::VariableArray variable_array = variable_array_prev;
    variable_array.liquid_phase_pressure = 1e3;

    ParameterLib::SpatialPosition x{};
    auto solution = model->integrateStress(variable_array_prev,
                                           variable_array,
                                           1.0,
                                           x,
                                           1.0,
                                           *state);
    ASSERT_TRUE(solution);

    auto& [forces_data, new_state, tangent_matrix] = *solution;
    ASSERT_TRUE(new_state != nullptr);

    MSM::OGSMFrontThermodynamicForcesView<
        3, boost::mp11::mp_list<MSM::Stress, MSM::Saturation>>
        view;
    auto const stress = view.block(MSM::stress, forces_data);
    auto const saturation = view.block(MSM::saturation, forces_data);

    double const expected_saturation = 1.0 / (1.0 + std::exp(-1.0));
    double const expected_stress_component = -expected_saturation * 1e3;

    EXPECT_NEAR(saturation, expected_saturation, 1e-12);
    EXPECT_NEAR(stress[0], expected_stress_component, 1e-12);
    EXPECT_NEAR(stress[1], expected_stress_component, 1e-12);
    EXPECT_NEAR(stress[2], expected_stress_component, 1e-12);
    EXPECT_NEAR(stress[3], 0.0, 1e-12);
    EXPECT_NEAR(stress[4], 0.0, 1e-12);
    EXPECT_NEAR(stress[5], 0.0, 1e-12);

    ASSERT_FALSE(tangent_matrix.data.empty());
    EXPECT_TRUE(std::all_of(
        tangent_matrix.data.begin(),
        tangent_matrix.data.end(),
        [](double const v) { return std::isfinite(v); }));
}

#endif  // OGS_USE_MFRONT
