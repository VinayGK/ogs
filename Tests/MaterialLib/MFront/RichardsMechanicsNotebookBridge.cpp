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
#include <unordered_map>
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

    add_param("E", 1e10);
    add_param("nu", 0.25);
    add_param("swelling_slope", 0.1);
    add_param("mass_exchange_coefficient", 1.0);
    add_param("rho_LR_ref", 1000.0);
    add_param("rho_l0", 1300.0);
    add_param("rho_lR0", 1300.0);
    add_param("rho_SR", 2470.0);
    add_param("density_a", 1.3);
    add_param("density_b", 1.0);
    add_param("hamaker_constant", -6e-20);
    add_param("specific_surface", 100.0);
    add_param("phi0", 0.2);
    add_param("area_factor_tuller", 1.0);
    add_param("pore_area_shape_factor_tuller", 0.8584073464102069);
    add_param("characteristic_pore_size", 1e-5);
    add_param("surface_tension", 0.0715);
    add_param("n_l0", 0.1);
    add_param("epsilon_sw0", 0.0);

    return parameters;
}

static std::string stripQuotes(std::string value)
{
    value.erase(0, value.find_first_not_of(" \t\r\n\"") );
    value.erase(value.find_last_not_of(" \t\r\n\"") + 1);
    return value;
}

static std::vector<std::string> splitCommaLine(std::string const& line)
{
    std::stringstream ss(line);
    std::string field;
    std::vector<std::string> fields;
    while (std::getline(ss, field, ','))
    {
        fields.push_back(stripQuotes(field));
    }
    return fields;
}

struct BaselineRow
{
    int step = 0;
    double pressure = 0.0;
    double saturation = 0.0;
    double mu_LR = 0.0;
    double n_l = 0.0;
    double phi_m = 0.0;
    double phi_M = 0.0;
    double rho_lR = 0.0;
    double mu_lR = 0.0;
    double rho_l_hat = 0.0;
    double epsilon_sw = 0.0;
    double stress_xx = 0.0;
};

static std::vector<BaselineRow> loadBaselineRows(std::string const& filename)
{
    std::ifstream in(filename);
    EXPECT_TRUE(in.good()) << filename;

    std::string header_line;
    std::getline(in, header_line);
    auto const headers = splitCommaLine(header_line);

    std::unordered_map<std::string, std::size_t> column;
    for (std::size_t i = 0; i < headers.size(); ++i)
    {
        column[headers[i]] = i;
    }

    auto const get_value = [&](std::vector<std::string> const& fields,
                               std::string const& key) -> double
    {
        return std::stod(fields.at(column.at(key)));
    };

    std::vector<BaselineRow> rows;
    std::string line;
    while (std::getline(in, line))
    {
        if (line.empty())
        {
            continue;
        }

        auto const fields = splitCommaLine(line);
        BaselineRow row;
        row.step = static_cast<int>(get_value(fields, "step"));
        row.pressure = get_value(fields, "pressure");
        row.saturation = get_value(fields, "S_L");
        row.mu_LR = get_value(fields, "mu_LR");
        row.n_l = get_value(fields, "n_l");
        row.phi_m = get_value(fields, "phi_m");
        row.phi_M = get_value(fields, "phi_M");
        row.rho_lR = get_value(fields, "rho_lR");
        row.mu_lR = get_value(fields, "mu_lR");
        row.rho_l_hat = get_value(fields, "rho_l_hat");
        row.epsilon_sw = get_value(fields, "epsilon_sw");
        row.stress_xx = get_value(fields, "sigma_S_xx");
        rows.push_back(row);
    }

    return rows;
}

static double scaledTolerance(double const value,
                              double const relative = 1e-8,
                              double const absolute = 1e-12)
{
    return absolute + relative * std::max(1.0, std::abs(value));
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
            <material_property name="SwellingSlope" parameter="swelling_slope"/>
            <material_property name="MassExchangeCoefficient" parameter="mass_exchange_coefficient"/>
            <material_property name="ReferenceLiquidDensityMacro" parameter="rho_LR_ref"/>
            <material_property name="ReferenceLiquidDensityMicro" parameter="rho_l0"/>
            <material_property name="ReferenceDensitySolid" parameter="rho_SR"/>
            <material_property name="MicroLiquidDensityA" parameter="density_a"/>
            <material_property name="MicroLiquidDensityB" parameter="density_b"/>
            <material_property name="HamakerConstant" parameter="hamaker_constant"/>
            <material_property name="SpecificSurface" parameter="specific_surface"/>
            <material_property name="InitialPorosity" parameter="phi0"/>
            <material_property name="AreaFactorTuller" parameter="area_factor_tuller"/>
            <material_property name="PoreAreaShapeFactorTuller" parameter="pore_area_shape_factor_tuller"/>
            <material_property name="CharacteristicPoreSize" parameter="characteristic_pore_size"/>
            <material_property name="SurfaceTension" parameter="surface_tension"/>
        </material_properties>
        <initial_values>
            <state_variable name="n_l" parameter="n_l0"/>
            <state_variable name="rho_lR" parameter="rho_lR0"/>
            <state_variable name="epsilon_sw" parameter="epsilon_sw0"/>
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
            <material_property name="SwellingSlope" parameter="swelling_slope"/>
            <material_property name="MassExchangeCoefficient" parameter="mass_exchange_coefficient"/>
            <material_property name="ReferenceLiquidDensityMacro" parameter="rho_LR_ref"/>
            <material_property name="ReferenceLiquidDensityMicro" parameter="rho_l0"/>
            <material_property name="ReferenceDensitySolid" parameter="rho_SR"/>
            <material_property name="MicroLiquidDensityA" parameter="density_a"/>
            <material_property name="MicroLiquidDensityB" parameter="density_b"/>
            <material_property name="HamakerConstant" parameter="hamaker_constant"/>
            <material_property name="SpecificSurface" parameter="specific_surface"/>
            <material_property name="InitialPorosity" parameter="phi0"/>
            <material_property name="AreaFactorTuller" parameter="area_factor_tuller"/>
            <material_property name="PoreAreaShapeFactorTuller" parameter="pore_area_shape_factor_tuller"/>
            <material_property name="CharacteristicPoreSize" parameter="characteristic_pore_size"/>
            <material_property name="SurfaceTension" parameter="surface_tension"/>
        </material_properties>
        <initial_values>
            <state_variable name="n_l" parameter="n_l0"/>
            <state_variable name="rho_lR" parameter="rho_lR0"/>
            <state_variable name="epsilon_sw" parameter="epsilon_sw0"/>
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

TEST(MaterialLib_RichardsMechanicsNotebookBridgeMFront,
     OneStepPressureResponse)
{
    auto const baseline_rows = loadBaselineRows(
        TestInfoLib::TestInfo::data_path +
        "/MaterialLib/MFront/RichardsMechanicsNotebookBridge_overlap_transfer_baseline.csv");
    ASSERT_EQ(baseline_rows.size(), 5);
    auto const& saturated_anchor = baseline_rows.front();

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
    variable_array_prev.liquid_saturation = 1.0;
    variable_array_prev.temperature = 0.0;

    MPL::VariableArray variable_array = variable_array_prev;
    variable_array.liquid_phase_pressure = 1000.0;

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

    EXPECT_NEAR(saturation, saturated_anchor.saturation,
                scaledTolerance(saturated_anchor.saturation));
    EXPECT_NEAR(stress[0], saturated_anchor.stress_xx,
                scaledTolerance(saturated_anchor.stress_xx, 5e-3, 1e-6));
    EXPECT_NEAR(stress[1], saturated_anchor.stress_xx,
                scaledTolerance(saturated_anchor.stress_xx, 5e-3, 1e-6));
    EXPECT_NEAR(stress[2], saturated_anchor.stress_xx,
                scaledTolerance(saturated_anchor.stress_xx, 5e-3, 1e-6));
    EXPECT_NEAR(stress[3], 0.0, 1e-12);
    EXPECT_NEAR(stress[4], 0.0, 1e-12);
    EXPECT_NEAR(stress[5], 0.0, 1e-12);

    ASSERT_FALSE(tangent_matrix.data.empty());
    EXPECT_TRUE(std::all_of(
        tangent_matrix.data.begin(), tangent_matrix.data.end(),
        [](double const v) { return std::isfinite(v); }));
}

TEST(MaterialLib_RichardsMechanicsNotebookBridgeMFront,
     FactoryPathOneStepPressureResponse)
{
    auto const baseline_rows = loadBaselineRows(
        TestInfoLib::TestInfo::data_path +
        "/MaterialLib/MFront/RichardsMechanicsNotebookBridge_overlap_transfer_baseline.csv");
    ASSERT_EQ(baseline_rows.size(), 5);
    auto const& saturated_anchor = baseline_rows.front();

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
    variable_array_prev.liquid_saturation = 1.0;
    variable_array_prev.temperature = 0.0;

    MPL::VariableArray variable_array = variable_array_prev;
    variable_array.liquid_phase_pressure = 1000.0;

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

    EXPECT_NEAR(stress[0], saturated_anchor.stress_xx,
                scaledTolerance(saturated_anchor.stress_xx, 5e-3, 1e-6));
    EXPECT_NEAR(stress[1], saturated_anchor.stress_xx,
                scaledTolerance(saturated_anchor.stress_xx, 5e-3, 1e-6));
    EXPECT_NEAR(stress[2], saturated_anchor.stress_xx,
                scaledTolerance(saturated_anchor.stress_xx, 5e-3, 1e-6));
    EXPECT_NEAR(stress[3], 0.0, 1e-12);
    EXPECT_NEAR(stress[4], 0.0, 1e-12);
    EXPECT_NEAR(stress[5], 0.0, 1e-12);
    EXPECT_TRUE(tangent_matrix.allFinite());
}

TEST(MaterialLib_RichardsMechanicsNotebookBridgeMFront,
     FactoryPathPressureCoupledBlocks)
{
    auto const baseline_rows = loadBaselineRows(
        TestInfoLib::TestInfo::data_path +
        "/MaterialLib/MFront/RichardsMechanicsNotebookBridge_overlap_transfer_baseline.csv");
    ASSERT_EQ(baseline_rows.size(), 5);
    auto const& saturated_anchor = baseline_rows.front();

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
    variable_array_prev.liquid_saturation = 1.0;
    variable_array_prev.temperature = 0.0;

    MPL::VariableArray variable_array = variable_array_prev;
    variable_array.liquid_phase_pressure = 1000.0;

    ParameterLib::SpatialPosition x{};
    auto response = model->integrateStressPressureCoupled(
        variable_array_prev, variable_array, 1.0, x, 1.0, *state);
    ASSERT_TRUE(response);
    ASSERT_TRUE(response->state != nullptr);

    EXPECT_NEAR(response->saturation, saturated_anchor.saturation,
                scaledTolerance(saturated_anchor.saturation));
    EXPECT_NEAR(response->stress[0], saturated_anchor.stress_xx,
                scaledTolerance(saturated_anchor.stress_xx, 5e-3, 1e-6));
    EXPECT_NEAR(response->stress[1], saturated_anchor.stress_xx,
                scaledTolerance(saturated_anchor.stress_xx, 5e-3, 1e-6));
    EXPECT_NEAR(response->stress[2], saturated_anchor.stress_xx,
                scaledTolerance(saturated_anchor.stress_xx, 5e-3, 1e-6));
    EXPECT_NEAR(response->stress[3], 0.0, 1e-12);
    EXPECT_NEAR(response->stress[4], 0.0, 1e-12);
    EXPECT_NEAR(response->stress[5], 0.0, 1e-12);
    EXPECT_NEAR(response->dStress_dLiquidPressure[0], 0.0, 1e-12);
    EXPECT_NEAR(response->dStress_dLiquidPressure[1], 0.0, 1e-12);
    EXPECT_NEAR(response->dStress_dLiquidPressure[2], 0.0, 1e-12);
    EXPECT_NEAR(response->dStress_dLiquidPressure[3], 0.0, 1e-12);
    EXPECT_NEAR(response->dStress_dLiquidPressure[4], 0.0, 1e-12);
    EXPECT_NEAR(response->dStress_dLiquidPressure[5], 0.0, 1e-12);
    EXPECT_NEAR(response->dSaturation_dLiquidPressure, 0.0, 1e-12);
    EXPECT_TRUE(response->dStress_dStrain.allFinite());
    for (int i = 0; i < response->dSaturation_dStrain.size(); ++i)
    {
        EXPECT_NEAR(response->dSaturation_dStrain[i], 0.0, 1e-12);
    }
}

TEST(MaterialLib_RichardsMechanicsNotebookBridgeMFront,
     PressureHistoryResponse)
{
    auto const baseline_rows = loadBaselineRows(
        TestInfoLib::TestInfo::data_path +
        "/MaterialLib/MFront/RichardsMechanicsNotebookBridge_overlap_transfer_baseline.csv");
    ASSERT_EQ(baseline_rows.size(), 5);

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
    variable_array_prev.liquid_saturation = 1.0;
    variable_array_prev.temperature = 0.0;

    ParameterLib::SpatialPosition x{};

    MSM::OGSMFrontThermodynamicForcesView<
        3, boost::mp11::mp_list<MSM::Stress, MSM::Saturation>>
        view;

    for (auto const& row : baseline_rows)
    {
        auto variable_array = variable_array_prev;
        variable_array.liquid_phase_pressure = row.pressure;

        auto solution = model->integrateStress(variable_array_prev,
                                               variable_array,
                                               static_cast<double>(row.step + 1),
                                               x,
                                               1.0,
                                               *state);
        ASSERT_TRUE(solution);

        auto& [forces_data, new_state, tangent_matrix] = *solution;
        ASSERT_TRUE(new_state != nullptr);

        auto const stress = view.block(MSM::stress, forces_data);
        auto const saturation = view.block(MSM::saturation, forces_data);

        EXPECT_NEAR(saturation, row.saturation, scaledTolerance(row.saturation));
        EXPECT_NEAR(stress[0], row.stress_xx, scaledTolerance(row.stress_xx, 5e-3, 1e-6));
        EXPECT_NEAR(stress[1], row.stress_xx, scaledTolerance(row.stress_xx, 5e-3, 1e-6));
        EXPECT_NEAR(stress[2], row.stress_xx, scaledTolerance(row.stress_xx, 5e-3, 1e-6));
        EXPECT_NEAR(stress[3], 0.0, 1e-12);
        EXPECT_NEAR(stress[4], 0.0, 1e-12);
        EXPECT_NEAR(stress[5], 0.0, 1e-12);

        ASSERT_FALSE(tangent_matrix.data.empty());
        EXPECT_TRUE(std::all_of(
            tangent_matrix.data.begin(), tangent_matrix.data.end(),
            [](double const v) { return std::isfinite(v); }));

        state = std::move(new_state);
        state->pushBackState();
        variable_array_prev = variable_array;
        variable_array_prev.stress.template emplace<KV>(stress);
        variable_array_prev.liquid_saturation = saturation;
    }
}

TEST(MaterialLib_RichardsMechanicsNotebookBridgeMFront,
     MicrostateHistoryResponse)
{
    auto const baseline_rows = loadBaselineRows(
        TestInfoLib::TestInfo::data_path +
        "/MaterialLib/MFront/RichardsMechanicsNotebookBridge_overlap_transfer_baseline.csv");
    ASSERT_EQ(baseline_rows.size(), 5);

    auto const parameters = createParameters();
    auto model = createBridgeModelThroughFactory(parameters);
    ASSERT_TRUE(model != nullptr);

    auto state = model->createMaterialStateVariables();
    ASSERT_TRUE(state != nullptr);
    initializeState(*model, *state);

    EXPECT_NEAR(getInternalVariable(*model, *state, "n_l"), 0.1, 1e-12);
    EXPECT_NEAR(getInternalVariable(*model, *state, "rho_lR"), 1300.0, 1e-12);
    EXPECT_NEAR(getInternalVariable(*model, *state, "epsilon_sw"), 0.0, 1e-12);

    MPL::VariableArray variable_array_prev;
    variable_array_prev.stress.template emplace<KV>(KV::Zero());
    variable_array_prev.mechanical_strain.template emplace<KV>(KV::Zero());
    variable_array_prev.liquid_phase_pressure = 0.0;
    variable_array_prev.liquid_saturation = 1.0;
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

        auto const n_l_value = getInternalVariable(*model, *response->state, "n_l");
        auto const phi_m_value = getInternalVariable(*model, *response->state, "phi_m");
        auto const phi_M_value = getInternalVariable(*model, *response->state, "phi_M");
        auto const rho_lR_value = getInternalVariable(*model, *response->state, "rho_lR");
        auto const mu_lR_value = getInternalVariable(*model, *response->state, "mu_lR");
        auto const rho_l_hat_value = getInternalVariable(*model, *response->state, "rho_l_hat");
        auto const epsilon_sw_value =
            getInternalVariable(*model, *response->state, "epsilon_sw");

        EXPECT_NEAR(response->saturation, row.saturation, scaledTolerance(row.saturation));
        EXPECT_NEAR(n_l_value, row.n_l, scaledTolerance(row.n_l, 5e-3, 1e-8));
        EXPECT_NEAR(phi_m_value, row.phi_m, scaledTolerance(row.phi_m, 5e-3, 1e-8));
        EXPECT_NEAR(phi_M_value, row.phi_M, scaledTolerance(row.phi_M, 5e-3, 1e-8));
        EXPECT_NEAR(rho_lR_value, row.rho_lR, scaledTolerance(row.rho_lR, 1e-3, 1e-6));
        EXPECT_NEAR(mu_lR_value, row.mu_lR, scaledTolerance(row.mu_lR, 1e-2, 1e-8));
        EXPECT_NEAR(rho_l_hat_value, row.rho_l_hat, scaledTolerance(row.rho_l_hat, 1e-2, 1e-8));
        EXPECT_NEAR(epsilon_sw_value, row.epsilon_sw, scaledTolerance(row.epsilon_sw, 5e-3, 1e-8));
        EXPECT_NEAR(phi_m_value, n_l_value, scaledTolerance(row.n_l, 5e-3, 1e-8));

        state = std::move(response->state);
        state->pushBackState();
        variable_array_prev = variable_array;
        variable_array_prev.stress.template emplace<KV>(response->stress);
        variable_array_prev.liquid_saturation = response->saturation;
    }
}

TEST(MaterialLib_RichardsMechanicsNotebookBridgeMFront,
     NotebookOverlapTransferBaselineHistory)
{
    auto const baseline_rows = loadBaselineRows(
        TestInfoLib::TestInfo::data_path +
        "/MaterialLib/MFront/RichardsMechanicsNotebookBridge_overlap_transfer_baseline.csv");
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
    variable_array_prev.liquid_saturation = 1.0;
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

        auto const n_l_value = getInternalVariable(*model, *response->state, "n_l");
        auto const phi_m_value = getInternalVariable(*model, *response->state, "phi_m");
        auto const phi_M_value = getInternalVariable(*model, *response->state, "phi_M");
        auto const rho_lR_value = getInternalVariable(*model, *response->state, "rho_lR");
        auto const mu_lR_value = getInternalVariable(*model, *response->state, "mu_lR");
        auto const rho_l_hat_value = getInternalVariable(*model, *response->state, "rho_l_hat");
        auto const epsilon_sw_value =
            getInternalVariable(*model, *response->state, "epsilon_sw");

        EXPECT_NEAR(response->saturation, row.saturation, scaledTolerance(row.saturation));
        EXPECT_NEAR(response->stress[0], row.stress_xx, scaledTolerance(row.stress_xx, 5e-3, 1e-6));
        EXPECT_NEAR(response->stress[1], row.stress_xx, scaledTolerance(row.stress_xx, 5e-3, 1e-6));
        EXPECT_NEAR(response->stress[2], row.stress_xx, scaledTolerance(row.stress_xx, 5e-3, 1e-6));
        EXPECT_NEAR(response->stress[3], 0.0, 1e-12);
        EXPECT_NEAR(response->stress[4], 0.0, 1e-12);
        EXPECT_NEAR(response->stress[5], 0.0, 1e-12);
        EXPECT_NEAR(response->dSaturation_dLiquidPressure, 0.0, 1e-12);
        EXPECT_NEAR(response->dStress_dLiquidPressure[0], 0.0, 1e-12);
        EXPECT_NEAR(response->dStress_dLiquidPressure[1], 0.0, 1e-12);
        EXPECT_NEAR(response->dStress_dLiquidPressure[2], 0.0, 1e-12);
        EXPECT_NEAR(response->dStress_dLiquidPressure[3], 0.0, 1e-12);
        EXPECT_NEAR(response->dStress_dLiquidPressure[4], 0.0, 1e-12);
        EXPECT_NEAR(response->dStress_dLiquidPressure[5], 0.0, 1e-12);
        EXPECT_NEAR(n_l_value, row.n_l, scaledTolerance(row.n_l, 5e-3, 1e-8));
        EXPECT_NEAR(phi_m_value, row.phi_m, scaledTolerance(row.phi_m, 5e-3, 1e-8));
        EXPECT_NEAR(phi_M_value, row.phi_M, scaledTolerance(row.phi_M, 5e-3, 1e-8));
        EXPECT_NEAR(rho_lR_value, row.rho_lR, scaledTolerance(row.rho_lR, 1e-3, 1e-6));
        EXPECT_NEAR(mu_lR_value, row.mu_lR, scaledTolerance(row.mu_lR, 1e-2, 1e-8));
        EXPECT_NEAR(rho_l_hat_value, row.rho_l_hat, scaledTolerance(row.rho_l_hat, 1e-2, 1e-8));
        EXPECT_NEAR(epsilon_sw_value, row.epsilon_sw, scaledTolerance(row.epsilon_sw, 5e-3, 1e-8));

        state = std::move(response->state);
        state->pushBackState();
        variable_array_prev = variable_array;
        variable_array_prev.stress.template emplace<KV>(response->stress);
        variable_array_prev.liquid_saturation = response->saturation;
    }
}

#endif  // OGS_USE_MFRONT
