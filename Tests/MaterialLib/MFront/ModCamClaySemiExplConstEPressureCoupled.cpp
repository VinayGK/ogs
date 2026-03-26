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
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "BaseLib/ConfigTree.h"
#include "MaterialLib/SolidModels/CreateConstitutiveRelation.h"
#include "ParameterLib/ConstantParameter.h"
#include "Tests/TestTools.h"

namespace MPL = MaterialPropertyLib;
using KV = MathLib::KelvinVector::KelvinVectorType<3>;
using KM = MathLib::KelvinVector::KelvinMatrixType<3>;
using MB = MaterialLib::Solids::MechanicsBase<3>;

namespace
{
std::vector<std::unique_ptr<ParameterLib::ParameterBase>>
createMCCBridgeParameters()
{
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> parameters;

    auto add_param = [&parameters](char const* name, double const value)
    {
        parameters.push_back(
            std::make_unique<ParameterLib::ConstantParameter<double>>(name,
                                                                      value));
    };

    add_param("YoungModulus", 52e6);
    add_param("PoissonRatio", 0.3);
    add_param("CriticalStateLineSlope", 1.2);
    add_param("SwellingLineSlope", 6.6e-3);
    add_param("VirginConsolidationLineSlope", 7.7e-2);
    add_param("InitialPreConsolidationPressure", 2e5);
    add_param("InitialVolumeRatio", 1.7857142857142858);
    add_param("ResidualLiquidSaturation", 0.0);
    add_param("ResidualGasSaturation", 0.0);
    add_param("BubblePressure", 1e4);
    add_param("VanGenuchtenExponent_m", 0.4);

    return parameters;
}

std::unique_ptr<MB> createMCCBridgeModelFromXml(
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    char const* const xml)
{
    auto local_coordinate_system = std::nullopt;

    auto ptree = Tests::readXml(xml);
    BaseLib::ConfigTree config_tree(std::move(ptree), "FILENAME",
                                    &BaseLib::ConfigTree::onerror,
                                    &BaseLib::ConfigTree::onwarning);

    return MaterialLib::Solids::createConstitutiveRelation<3>(
        parameters, local_coordinate_system, config_tree);
}

std::unique_ptr<MB> createMCCBridgeNativeModel(
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters)
{
    char const* xml = R"XML(
        <type>MFront</type>
        <behaviour>ModCamClay_semiExpl_constE</behaviour>
        <library path_is_relative_to_prj_file="false">libOgsMFrontBehaviour</library>
        <material_properties>
            <material_property name="YoungModulus" parameter="YoungModulus"/>
            <material_property name="PoissonRatio" parameter="PoissonRatio"/>
            <material_property name="CriticalStateLineSlope" parameter="CriticalStateLineSlope"/>
            <material_property name="SwellingLineSlope" parameter="SwellingLineSlope"/>
            <material_property name="VirginConsolidationLineSlope" parameter="VirginConsolidationLineSlope"/>
            <material_property name="CharacteristicPreConsolidationPressure" parameter="InitialPreConsolidationPressure"/>
        </material_properties>
        <initial_values>
            <state_variable name="PreConsolidationPressure" parameter="InitialPreConsolidationPressure"/>
            <state_variable name="VolumeRatio" parameter="InitialVolumeRatio"/>
        </initial_values>
        )XML";

    return createMCCBridgeModelFromXml(parameters, xml);
}

std::unique_ptr<MB> createMCCBridgePressureCoupledModel(
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters)
{
    char const* xml = R"XML(
        <type>MFrontRichardsMechanics</type>
        <behaviour>ModCamClay_semiExpl_constE_pressureCoupled</behaviour>
        <library path_is_relative_to_prj_file="false">libOgsMFrontBehaviour</library>
        <material_properties>
            <material_property name="YoungModulus" parameter="YoungModulus"/>
            <material_property name="PoissonRatio" parameter="PoissonRatio"/>
            <material_property name="CriticalStateLineSlope" parameter="CriticalStateLineSlope"/>
            <material_property name="SwellingLineSlope" parameter="SwellingLineSlope"/>
            <material_property name="VirginConsolidationLineSlope" parameter="VirginConsolidationLineSlope"/>
            <material_property name="CharacteristicPreConsolidationPressure" parameter="InitialPreConsolidationPressure"/>
            <material_property name="ResidualLiquidSaturation" parameter="ResidualLiquidSaturation"/>
            <material_property name="ResidualGasSaturation" parameter="ResidualGasSaturation"/>
            <material_property name="BubblePressure" parameter="BubblePressure"/>
            <material_property name="VanGenuchtenExponent_m" parameter="VanGenuchtenExponent_m"/>
        </material_properties>
        <initial_values>
            <state_variable name="PreConsolidationPressure" parameter="InitialPreConsolidationPressure"/>
            <state_variable name="VolumeRatio" parameter="InitialVolumeRatio"/>
        </initial_values>
        )XML";

    return createMCCBridgeModelFromXml(parameters, xml);
}

void initializeMCCBridgeState(MB const& model,
                              MB::MaterialStateVariables& state)
{
    ParameterLib::SpatialPosition x{};
    model.initializeInternalStateVariables(0.0, x, state);
}

double getMCCBridgeInternalScalar(MB const& model,
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

std::vector<double> getMCCBridgeInternalVector(
    MB const& model,
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
        return {};
    }

    std::vector<double> cache;
    auto const& values = it->getter(state, cache);
    return {values.begin(), values.end()};
}

KV mccBridgeIsotropicStress(double const value)
{
    KV sigma = KV::Zero();
    sigma[0] = value;
    sigma[1] = value;
    sigma[2] = value;
    return sigma;
}

KV mccBridgeIsotropicStrainFromVolumetric(double const epsilon_v)
{
    KV eps = KV::Zero();
    eps[0] = epsilon_v / 3.0;
    eps[1] = epsilon_v / 3.0;
    eps[2] = epsilon_v / 3.0;
    return eps;
}

double mccBridgeSaturationVG(double const p_cap,
                             double const residual_liquid_saturation,
                             double const residual_gas_saturation,
                             double const exponent_m,
                             double const bubble_pressure)
{
    double const saturation_max = 1.0 - residual_gas_saturation;
    if (p_cap <= 0.0)
    {
        return saturation_max;
    }

    double const exponent_n = 1.0 / (1.0 - exponent_m);
    double const p = p_cap / bubble_pressure;
    double const p_to_n = std::pow(p, exponent_n);
    double const s_eff = std::pow(p_to_n + 1.0, -exponent_m);
    double const s =
        s_eff * (saturation_max - residual_liquid_saturation) +
        residual_liquid_saturation;
    return std::clamp(s, residual_liquid_saturation, saturation_max);
}

double mccBridgeDSaturationVGdLiquidPressure(
    double const liquid_pressure,
    double const residual_liquid_saturation,
    double const residual_gas_saturation,
    double const exponent_m,
    double const bubble_pressure)
{
    double const p_cap = -liquid_pressure;
    if (p_cap <= 0.0)
    {
        return 0.0;
    }

    double const exponent_n = 1.0 / (1.0 - exponent_m);
    double const saturation_max = 1.0 - residual_gas_saturation;
    double const p = p_cap / bubble_pressure;
    double const p_to_n = std::pow(p, exponent_n);
    double const s_eff = std::pow(p_to_n + 1.0, -exponent_m);
    double const ds_eff_dp_cap =
        -exponent_m * exponent_n * p_to_n * s_eff /
        (p_cap * (p_to_n + 1.0));

    return -ds_eff_dp_cap *
           (saturation_max - residual_liquid_saturation);
}
}  // namespace

TEST(MaterialLib_RMBridgeMFront_MCC,
     PressureCoupledConstEMatchesNativeMechanicalPath)
{
    auto parameters = createMCCBridgeParameters();
    auto native = createMCCBridgeNativeModel(parameters);
    auto bridge = createMCCBridgePressureCoupledModel(parameters);

    auto native_state = native->createMaterialStateVariables();
    auto bridge_state = bridge->createMaterialStateVariables();
    initializeMCCBridgeState(*native, *native_state);
    initializeMCCBridgeState(*bridge, *bridge_state);

    ParameterLib::SpatialPosition x{};
    double const liquid_pressure_prev = -5e3;
    double const residual_liquid_saturation = 0.0;
    double const residual_gas_saturation = 0.0;
    double const exponent_m = 0.4;
    double const bubble_pressure = 1e4;

    MPL::VariableArray native_prev;
    native_prev.mechanical_strain.emplace<KV>(KV::Zero());
    native_prev.stress.emplace<KV>(mccBridgeIsotropicStress(-5e3));
    native_prev.liquid_phase_pressure = liquid_pressure_prev;
    native_prev.temperature = 293.15;

    MPL::VariableArray bridge_prev = native_prev;

    struct Step
    {
        double t;
        double eps_v;
        double liquid_pressure;
    };

    std::vector<Step> const steps = {
        {1.0, -5e-4, -2e4},
        {2.0, -1.5e-3, -4e4},
        {3.0, -3.5e-3, -7e4},
        {4.0, -6.5e-3, -1e5},
    };

    for (std::size_t i = 0; i < steps.size(); ++i)
    {
        MPL::VariableArray native_current = native_prev;
        native_current.mechanical_strain.emplace<KV>(
            mccBridgeIsotropicStrainFromVolumetric(steps[i].eps_v));
        native_current.liquid_phase_pressure = steps[i].liquid_pressure;

        MPL::VariableArray bridge_current = bridge_prev;
        bridge_current.mechanical_strain.emplace<KV>(
            mccBridgeIsotropicStrainFromVolumetric(steps[i].eps_v));
        bridge_current.liquid_phase_pressure = steps[i].liquid_pressure;

        double const dt = i == 0 ? steps[i].t : steps[i].t - steps[i - 1].t;

        auto native_solution = native->integrateStress(
            native_prev, native_current, steps[i].t, x, dt, *native_state);
        ASSERT_TRUE(native_solution);

        auto bridge_solution = bridge->integrateStressPressureCoupled(
            bridge_prev, bridge_current, steps[i].t, x, dt, *bridge_state);
        ASSERT_TRUE(bridge_solution);

        auto& [native_stress, native_state_new, native_tangent] =
            *native_solution;
        auto& bridge_response = *bridge_solution;

        EXPECT_TRUE(
            (native_stress - bridge_response.stress).isZero(5e-11));
        EXPECT_TRUE((native_tangent - bridge_response.dStress_dStrain)
                        .isZero(1e-8));

        EXPECT_TRUE(bridge_response.dStress_dLiquidPressure.isZero(1e-12));
        EXPECT_TRUE(bridge_response.dSaturation_dStrain.isZero(1e-12));

        auto const expected_saturation = mccBridgeSaturationVG(
            -steps[i].liquid_pressure, residual_liquid_saturation,
            residual_gas_saturation, exponent_m, bubble_pressure);
        auto const expected_dS_dp = mccBridgeDSaturationVGdLiquidPressure(
            steps[i].liquid_pressure, residual_liquid_saturation,
            residual_gas_saturation, exponent_m, bubble_pressure);

        EXPECT_NEAR(expected_saturation, bridge_response.saturation, 1e-15);
        EXPECT_NEAR(expected_dS_dp, bridge_response.dSaturation_dLiquidPressure,
                    1e-15);

        EXPECT_NEAR(
            getMCCBridgeInternalScalar(*native, *native_state_new,
                                       "EquivalentPlasticStrain"),
            getMCCBridgeInternalScalar(*bridge, *bridge_response.state,
                                       "EquivalentPlasticStrain"),
            1e-14);
        EXPECT_NEAR(
            getMCCBridgeInternalScalar(*native, *native_state_new,
                                       "PreConsolidationPressure"),
            getMCCBridgeInternalScalar(*bridge, *bridge_response.state,
                                       "PreConsolidationPressure"),
            1e-8);
        EXPECT_NEAR(
            getMCCBridgeInternalScalar(*native, *native_state_new,
                                       "PlasticVolumetricStrain"),
            getMCCBridgeInternalScalar(*bridge, *bridge_response.state,
                                       "PlasticVolumetricStrain"),
            1e-14);
        EXPECT_NEAR(getMCCBridgeInternalScalar(*native, *native_state_new,
                                               "VolumeRatio"),
                    getMCCBridgeInternalScalar(*bridge,
                                               *bridge_response.state,
                                               "VolumeRatio"),
                    1e-14);

        auto const native_elastic_strain = getMCCBridgeInternalVector(
            *native, *native_state_new, "ElasticStrain");
        auto const bridge_elastic_strain = getMCCBridgeInternalVector(
            *bridge, *bridge_response.state, "ElasticStrain");
        ASSERT_EQ(native_elastic_strain.size(), bridge_elastic_strain.size());
        for (std::size_t k = 0; k < native_elastic_strain.size(); ++k)
        {
            EXPECT_NEAR(native_elastic_strain[k], bridge_elastic_strain[k],
                        1e-14);
        }

        native_state = std::move(native_state_new);
        bridge_state = std::move(bridge_response.state);

        native_prev = native_current;
        native_prev.stress.emplace<KV>(native_stress);
        bridge_prev = bridge_current;
        bridge_prev.stress.emplace<KV>(bridge_response.stress);
    }
}

TEST(MaterialLib_RMBridgeMFront_MCC, PressureCoupledConstEFactoryCreationWorks)
{
    auto parameters = createMCCBridgeParameters();
    auto bridge = createMCCBridgePressureCoupledModel(parameters);
    ASSERT_TRUE(bridge);

    auto state = bridge->createMaterialStateVariables();
    initializeMCCBridgeState(*bridge, *state);

    MPL::VariableArray previous;
    previous.mechanical_strain.emplace<KV>(KV::Zero());
    previous.stress.emplace<KV>(mccBridgeIsotropicStress(-5e3));
    previous.liquid_phase_pressure = -5e3;
    previous.temperature = 293.15;

    MPL::VariableArray current = previous;
    current.mechanical_strain.emplace<KV>(
        mccBridgeIsotropicStrainFromVolumetric(-1e-3));
    current.liquid_phase_pressure = -2e4;

    ParameterLib::SpatialPosition x{};
    auto response = bridge->integrateStressPressureCoupled(previous, current,
                                                           1.0, x, 1.0,
                                                           *state);
    ASSERT_TRUE(response);
    EXPECT_TRUE(response->stress.allFinite());
    EXPECT_TRUE(response->dStress_dStrain.allFinite());
    EXPECT_TRUE(response->dStress_dLiquidPressure.allFinite());
    EXPECT_TRUE(std::isfinite(response->saturation));
    EXPECT_TRUE(std::isfinite(response->dSaturation_dLiquidPressure));
    EXPECT_TRUE(response->dSaturation_dStrain.allFinite());
}

#endif  // OGS_USE_MFRONT
