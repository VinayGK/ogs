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
#include <fstream>
#include <limits>
#include <memory>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

#include <gtest/gtest.h>
#include <boost/mp11.hpp>

#include "BaseLib/ConfigTree.h"
#include "InfoLib/TestInfo.h"
#include "MaterialLib/SolidModels/CreateConstitutiveRelation.h"
#include "MaterialLib/SolidModels/MFront/CreateMFrontGeneric.h"
#include "MaterialLib/SolidModels/MFront/Variable.h"
#include "ParameterLib/ConstantParameter.h"
#include "Tests/TestTools.h"

namespace MSM = MaterialLib::Solids::MFront;
namespace MPL = MaterialPropertyLib;
using KV = MathLib::KelvinVector::KelvinVectorType<3>;
using MB = MaterialLib::Solids::MechanicsBase<3>;

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
    add_param("n_l_ref", 0.2);
    add_param("n_l0", 0.1);

    return parameters;
}


static double saturationFromPressure(double const pressure,
                                     double const p_sat_scale)
{
    return 1.0 / (1.0 + std::exp(-pressure / p_sat_scale));
}

static double dSaturationDp(double const pressure, double const p_sat_scale)
{
    auto const saturation = saturationFromPressure(pressure, p_sat_scale);
    return saturation * (1.0 - saturation) / p_sat_scale;
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
            <material_property name="n_l_ref" parameter="n_l_ref"/>
        </material_properties>
        <initial_values>
            <state_variable name="PreConsolidationPressure" parameter="pc0"/>
            <state_variable name="VolumeRatio" parameter="v0"/>
            <state_variable name="n_l" parameter="n_l0"/>
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

static auto createBridgeModelThroughFactory(
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
            <material_property name="n_l_ref" parameter="n_l_ref"/>
        </material_properties>
        <initial_values>
            <state_variable name="PreConsolidationPressure" parameter="pc0"/>
            <state_variable name="VolumeRatio" parameter="v0"/>
            <state_variable name="n_l" parameter="n_l0"/>
        </initial_values>
        )XML";

    auto ptree = Tests::readXml(xml);
    BaseLib::ConfigTree config_tree(std::move(ptree), "FILENAME",
                                    &BaseLib::ConfigTree::onerror,
                                    &BaseLib::ConfigTree::onwarning);

    return MaterialLib::Solids::createConstitutiveRelation<3>(
        parameters, local_coordinate_system, config_tree);
}

template <typename Model>
static void initializeState(Model const& model,
                            MB::MaterialStateVariables& state)
{
    ParameterLib::SpatialPosition x{};
    model.initializeInternalStateVariables(0.0, x, state);
}

template <typename Model>
static double getInternalVariable(Model const& model,
                                  MB::MaterialStateVariables const& state,
                                  std::string const& name)
{
    auto const internal_variables = model.getInternalVariables();
    auto const it = std::find_if(
        internal_variables.begin(), internal_variables.end(),
        [&name](auto const& internal_variable)
        { return internal_variable.name == name; });
    if (it == internal_variables.end())
    {
        return std::numeric_limits<double>::quiet_NaN();
    }

    std::vector<double> cache;
    auto const& values = it->getter(state, cache);
    if (values.size() != 1)
    {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return values[0];
}

struct BaselineRow
{
    int step = 0;
    double pressure = 0.0;
    double saturation = 0.0;
    double stress_xx = 0.0;
    double n_l = 0.0;
    double phi_m = 0.0;
    double dS_dp = 0.0;
    double dSigma_dp = 0.0;
};

static std::vector<BaselineRow> loadBaselineRows(std::string const& filename)
{
    std::ifstream in(filename);
    EXPECT_TRUE(in.good()) << filename;

    std::vector<BaselineRow> rows;
    std::string line;

    std::getline(in, line);  // header
    while (std::getline(in, line))
    {
        if (line.empty())
        {
            continue;
        }

        std::stringstream ss(line);
        std::string field;
        BaselineRow row;

        std::getline(ss, field, ',');
        row.step = std::stoi(field);
        std::getline(ss, field, ',');
        row.pressure = std::stod(field);
        std::getline(ss, field, ',');
        row.saturation = std::stod(field);
        std::getline(ss, field, ',');
        row.stress_xx = std::stod(field);
        std::getline(ss, field, ',');
        row.n_l = std::stod(field);
        std::getline(ss, field, ',');
        row.phi_m = std::stod(field);
        std::getline(ss, field, ',');
        row.dS_dp = std::stod(field);
        std::getline(ss, field, ',');
        row.dSigma_dp = std::stod(field);

        rows.push_back(row);
    }

    return rows;
}

TEST(MaterialLib_RichardsMechanicsNotebookBridgeMFront,
     OneStepPressureResponse)
{
    auto const parameters = createParameters();
    auto model = createBridgeModel(parameters);
    ASSERT_TRUE(model != nullptr);

    auto state = model->createMaterialStateVariables();
    ASSERT_TRUE(state != nullptr);
    initializeState(*model, *state);

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

TEST(MaterialLib_RichardsMechanicsNotebookBridgeMFront,
     FactoryPathOneStepPressureResponse)
{
    auto const parameters = createParameters();
    auto model = createBridgeModelThroughFactory(parameters);
    ASSERT_TRUE(model != nullptr);

    auto state = model->createMaterialStateVariables();
    ASSERT_TRUE(state != nullptr);
    initializeState(*model, *state);

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

    auto& [stress, new_state, tangent_matrix] = *solution;
    ASSERT_TRUE(new_state != nullptr);

    double const expected_saturation = saturationFromPressure(1e3, 1e3);
    double const expected_stress_component = -expected_saturation * 1e3;

    EXPECT_NEAR(stress[0], expected_stress_component, 1e-12);
    EXPECT_NEAR(stress[1], expected_stress_component, 1e-12);
    EXPECT_NEAR(stress[2], expected_stress_component, 1e-12);
    EXPECT_NEAR(stress[3], 0.0, 1e-12);
    EXPECT_NEAR(stress[4], 0.0, 1e-12);
    EXPECT_NEAR(stress[5], 0.0, 1e-12);

    EXPECT_TRUE(tangent_matrix.allFinite());
}


TEST(MaterialLib_RichardsMechanicsNotebookBridgeMFront,
     FactoryPathPressureCoupledBlocks)
{
    auto const parameters = createParameters();
    auto model = createBridgeModelThroughFactory(parameters);
    ASSERT_TRUE(model != nullptr);

    auto state = model->createMaterialStateVariables();
    ASSERT_TRUE(state != nullptr);
    initializeState(*model, *state);

    MPL::VariableArray variable_array_prev;
    variable_array_prev.stress.template emplace<KV>(KV::Zero());
    variable_array_prev.mechanical_strain.template emplace<KV>(KV::Zero());
    variable_array_prev.liquid_phase_pressure = 0.0;
    variable_array_prev.liquid_saturation = 0.0;
    variable_array_prev.temperature = 0.0;

    MPL::VariableArray variable_array = variable_array_prev;
    variable_array.liquid_phase_pressure = 1e3;

    ParameterLib::SpatialPosition x{};
    auto response = model->integrateStressPressureCoupled(
        variable_array_prev, variable_array, 1.0, x, 1.0, *state);
    ASSERT_TRUE(response);
    ASSERT_TRUE(response->state != nullptr);

    double const expected_saturation = saturationFromPressure(1e3, 1e3);
    double const expected_dS_dp = dSaturationDp(1e3, 1e3);
    double const expected_dsigma_dp = -(expected_saturation + 1e3 * expected_dS_dp);

    EXPECT_NEAR(response->saturation, expected_saturation, 1e-12);
    EXPECT_NEAR(response->stress[0], -expected_saturation * 1e3, 1e-12);
    EXPECT_NEAR(response->stress[1], -expected_saturation * 1e3, 1e-12);
    EXPECT_NEAR(response->stress[2], -expected_saturation * 1e3, 1e-12);
    EXPECT_NEAR(response->dStress_dLiquidPressure[0], expected_dsigma_dp, 1e-12);
    EXPECT_NEAR(response->dStress_dLiquidPressure[1], expected_dsigma_dp, 1e-12);
    EXPECT_NEAR(response->dStress_dLiquidPressure[2], expected_dsigma_dp, 1e-12);
    EXPECT_NEAR(response->dStress_dLiquidPressure[3], 0.0, 1e-12);
    EXPECT_NEAR(response->dStress_dLiquidPressure[4], 0.0, 1e-12);
    EXPECT_NEAR(response->dStress_dLiquidPressure[5], 0.0, 1e-12);
    EXPECT_NEAR(response->dSaturation_dLiquidPressure, expected_dS_dp, 1e-12);
    EXPECT_TRUE(response->dStress_dStrain.allFinite());
    for (int i = 0; i < response->dSaturation_dStrain.size(); ++i)
    {
        EXPECT_NEAR(response->dSaturation_dStrain[i], 0.0, 1e-12);
    }
}

TEST(MaterialLib_RichardsMechanicsNotebookBridgeMFront,
     PressureHistoryResponse)
{
    auto const parameters = createParameters();
    auto model = createBridgeModel(parameters);
    ASSERT_TRUE(model != nullptr);

    auto state = model->createMaterialStateVariables();
    ASSERT_TRUE(state != nullptr);
    initializeState(*model, *state);

    MPL::VariableArray variable_array_prev;
    variable_array_prev.stress.template emplace<KV>(KV::Zero());
    variable_array_prev.mechanical_strain.template emplace<KV>(KV::Zero());
    variable_array_prev.liquid_phase_pressure = 0.0;
    variable_array_prev.liquid_saturation = saturationFromPressure(0.0, 1e3);
    variable_array_prev.temperature = 0.0;

    std::vector<double> const pressures{0.0, 250.0, 500.0, 750.0, 1000.0};
    ParameterLib::SpatialPosition x{};

    double previous_pressure = 0.0;
    double previous_saturation = saturationFromPressure(previous_pressure, 1e3);
    KV previous_stress = KV::Zero();

    MSM::OGSMFrontThermodynamicForcesView<
        3, boost::mp11::mp_list<MSM::Stress, MSM::Saturation>>
        view;

    for (std::size_t step = 0; step < pressures.size(); ++step)
    {
        auto variable_array = variable_array_prev;
        variable_array.liquid_phase_pressure = pressures[step];

        auto solution = model->integrateStress(variable_array_prev,
                                               variable_array,
                                               static_cast<double>(step + 1),
                                               x,
                                               1.0,
                                               *state);
        ASSERT_TRUE(solution);

        auto& [forces_data, new_state, tangent_matrix] = *solution;
        ASSERT_TRUE(new_state != nullptr);

        auto const stress = view.block(MSM::stress, forces_data);
        auto const saturation = view.block(MSM::saturation, forces_data);

        auto const expected_saturation =
            saturationFromPressure(pressures[step], 1e3);
        auto expected_stress = previous_stress;
        auto const increment = -(expected_saturation * pressures[step] -
                                 previous_saturation * previous_pressure);
        expected_stress[0] += increment;
        expected_stress[1] += increment;
        expected_stress[2] += increment;

        EXPECT_NEAR(saturation, expected_saturation, 1e-12);
        EXPECT_NEAR(stress[0], expected_stress[0], 1e-12);
        EXPECT_NEAR(stress[1], expected_stress[1], 1e-12);
        EXPECT_NEAR(stress[2], expected_stress[2], 1e-12);
        EXPECT_NEAR(stress[3], 0.0, 1e-12);
        EXPECT_NEAR(stress[4], 0.0, 1e-12);
        EXPECT_NEAR(stress[5], 0.0, 1e-12);

        ASSERT_FALSE(tangent_matrix.data.empty());
        EXPECT_TRUE(std::all_of(
            tangent_matrix.data.begin(), tangent_matrix.data.end(),
            [](double const v) { return std::isfinite(v); }));

        previous_pressure = pressures[step];
        previous_saturation = saturation;
        previous_stress = stress;
        variable_array_prev = variable_array;
        variable_array_prev.stress.template emplace<KV>(stress);
        variable_array_prev.liquid_saturation = saturation;
        state = std::move(new_state);
        state->pushBackState();
    }
}


TEST(MaterialLib_RichardsMechanicsNotebookBridgeMFront,
     MicrostateHistoryResponse)
{
    auto const parameters = createParameters();
    auto model = createBridgeModelThroughFactory(parameters);
    ASSERT_TRUE(model != nullptr);

    auto state = model->createMaterialStateVariables();
    ASSERT_TRUE(state != nullptr);
    initializeState(*model, *state);

    EXPECT_NEAR(getInternalVariable(*model, *state, "n_l"), 0.1, 1e-12);

    MPL::VariableArray variable_array_prev;
    variable_array_prev.stress.template emplace<KV>(KV::Zero());
    variable_array_prev.mechanical_strain.template emplace<KV>(KV::Zero());
    variable_array_prev.liquid_phase_pressure = 0.0;
    variable_array_prev.liquid_saturation = saturationFromPressure(0.0, 1e3);
    variable_array_prev.temperature = 0.0;

    std::vector<double> const pressures{250.0, 500.0, 750.0, 1000.0};
    ParameterLib::SpatialPosition x{};

    double previous_n_l = 0.1;

    for (std::size_t step = 0; step < pressures.size(); ++step)
    {
        auto variable_array = variable_array_prev;
        variable_array.liquid_phase_pressure = pressures[step];

        auto response = model->integrateStressPressureCoupled(
            variable_array_prev,
            variable_array,
            static_cast<double>(step + 1),
            x,
            1.0,
            *state);
        ASSERT_TRUE(response);
        ASSERT_TRUE(response->state != nullptr);

        double const expected_saturation =
            saturationFromPressure(pressures[step], 1e3);
        double const expected_n_l = 0.2 * expected_saturation;

        EXPECT_NEAR(response->saturation, expected_saturation, 1e-12);

        auto const n_l_value =
            getInternalVariable(*model, *response->state, "n_l");
        auto const phi_m_value =
            getInternalVariable(*model, *response->state, "phi_m");

        ASSERT_TRUE(std::isfinite(n_l_value));
        ASSERT_TRUE(std::isfinite(phi_m_value));
        EXPECT_NEAR(n_l_value, expected_n_l, 1e-12);
        EXPECT_NEAR(phi_m_value, expected_n_l, 1e-12);
        EXPECT_NEAR(phi_m_value, n_l_value, 1e-12);
        EXPECT_GE(n_l_value, previous_n_l);

        previous_n_l = n_l_value;
        state = std::move(response->state);
        state->pushBackState();
        variable_array_prev = variable_array;
        variable_array_prev.stress.template emplace<KV>(response->stress);
        variable_array_prev.liquid_saturation = response->saturation;
    }
}

TEST(MaterialLib_RichardsMechanicsNotebookBridgeMFront,
     ReducedNotebookBaselineHistory)
{
    auto const baseline_rows = loadBaselineRows(
        TestInfoLib::TestInfo::data_path +
        "/MaterialLib/MFront/RichardsMechanicsNotebookBridge_history_baseline.csv");
    ASSERT_EQ(baseline_rows.size(), 5);

    auto const parameters = createParameters();
    auto model = createBridgeModelThroughFactory(parameters);
    ASSERT_TRUE(model != nullptr);

    auto state = model->createMaterialStateVariables();
    ASSERT_TRUE(state != nullptr);
    initializeState(*model, *state);

    MPL::VariableArray variable_array_prev;
    variable_array_prev.stress.template emplace<KV>(KV::Zero());
    variable_array_prev.mechanical_strain.template emplace<KV>(KV::Zero());
    variable_array_prev.liquid_phase_pressure = 0.0;
    variable_array_prev.liquid_saturation = saturationFromPressure(0.0, 1e3);
    variable_array_prev.temperature = 0.0;

    ParameterLib::SpatialPosition x{};

    for (auto const& row : baseline_rows)
    {
        auto variable_array = variable_array_prev;
        variable_array.liquid_phase_pressure = row.pressure;

        auto response = model->integrateStressPressureCoupled(
            variable_array_prev,
            variable_array,
            static_cast<double>(row.step + 1),
            x,
            1.0,
            *state);
        ASSERT_TRUE(response);
        ASSERT_TRUE(response->state != nullptr);

        auto const n_l_value =
            getInternalVariable(*model, *response->state, "n_l");
        auto const phi_m_value =
            getInternalVariable(*model, *response->state, "phi_m");

        EXPECT_NEAR(response->saturation, row.saturation, 1e-12);
        EXPECT_NEAR(response->stress[0], row.stress_xx, 1e-12);
        EXPECT_NEAR(response->stress[1], row.stress_xx, 1e-12);
        EXPECT_NEAR(response->stress[2], row.stress_xx, 1e-12);
        EXPECT_NEAR(response->stress[3], 0.0, 1e-12);
        EXPECT_NEAR(response->stress[4], 0.0, 1e-12);
        EXPECT_NEAR(response->stress[5], 0.0, 1e-12);
        EXPECT_NEAR(response->dSaturation_dLiquidPressure, row.dS_dp, 1e-12);
        EXPECT_NEAR(response->dStress_dLiquidPressure[0], row.dSigma_dp, 1e-12);
        EXPECT_NEAR(response->dStress_dLiquidPressure[1], row.dSigma_dp, 1e-12);
        EXPECT_NEAR(response->dStress_dLiquidPressure[2], row.dSigma_dp, 1e-12);
        EXPECT_NEAR(response->dStress_dLiquidPressure[3], 0.0, 1e-12);
        EXPECT_NEAR(response->dStress_dLiquidPressure[4], 0.0, 1e-12);
        EXPECT_NEAR(response->dStress_dLiquidPressure[5], 0.0, 1e-12);
        EXPECT_NEAR(n_l_value, row.n_l, 1e-12);
        EXPECT_NEAR(phi_m_value, row.phi_m, 1e-12);
        EXPECT_NEAR(phi_m_value, n_l_value, 1e-12);

        state = std::move(response->state);
        state->pushBackState();
        variable_array_prev = variable_array;
        variable_array_prev.stress.template emplace<KV>(response->stress);
        variable_array_prev.liquid_saturation = response->saturation;
    }
}

#endif  // OGS_USE_MFRONT
