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
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "BaseLib/ConfigTree.h"
#include "MaterialLib/SolidModels/CreateConstitutiveRelation.h"
#include "MaterialLib/SolidModels/MFront/MFrontGeneric.h"
#include "ParameterLib/ConstantParameter.h"
#include "Tests/TestTools.h"

namespace MPL = MaterialPropertyLib;
namespace MSM = MaterialLib::Solids::MFront;
using KV = MathLib::KelvinVector::KelvinVectorType<3>;
using KM = MathLib::KelvinVector::KelvinMatrixType<3>;
using MB = MaterialLib::Solids::MechanicsBase<3>;

namespace
{
std::vector<std::unique_ptr<ParameterLib::ParameterBase>>
createMCCBridgeParameters(
    double const initial_volume_ratio = 1.7857142857142858)
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
    add_param("InitialVolumeRatio", initial_volume_ratio);
    add_param("ResidualLiquidSaturation", 0.0);
    add_param("ResidualGasSaturation", 0.0);
    add_param("BubblePressure", 1e4);
    add_param("VanGenuchtenExponent_m", 0.4);

    return parameters;
}

double mccBridgeInitialVolumeRatioFromDryDensity(
    double const dry_density,
    double const solid_density = 2470.0)
{
    return solid_density / dry_density;
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

bool copyMCCBridgeInternalVariables(MB const& source_model,
                                    MB::MaterialStateVariables const& source,
                                    MB const& target_model,
                                    MB::MaterialStateVariables& target)
{
    auto const source_internal_variables = source_model.getInternalVariables();
    auto const target_internal_variables = target_model.getInternalVariables();

    for (auto const& source_variable : source_internal_variables)
    {
        auto const target_it = std::find_if(
            target_internal_variables.begin(), target_internal_variables.end(),
            [&source_variable](auto const& target_variable)
            { return target_variable.name == source_variable.name; });

        if (target_it == target_internal_variables.end())
        {
            return false;
        }

        std::vector<double> cache;
        auto const& source_values = source_variable.getter(source, cache);
        auto target_values = target_it->reference(target);
        if (target_values.size() != source_values.size())
        {
            return false;
        }

        std::copy(source_values.begin(), source_values.end(),
                  target_values.begin());
    }

    return true;
}

void setMCCBridgeThermodynamicForces(MB::MaterialStateVariables& state,
                                     KV const& stress,
                                     double const saturation)
{
    auto* const mfront_state =
        dynamic_cast<MSM::MaterialStateVariablesMFront<3>*>(&state);
    if (mfront_state == nullptr)
    {
        return;
    }

    auto& thermodynamic_forces =
        mfront_state->_behaviour_data.s1.thermodynamic_forces;
    thermodynamic_forces[0] = stress[0];
    thermodynamic_forces[1] = stress[1];
    thermodynamic_forces[2] = stress[2];
    thermodynamic_forces[3] = stress[3];
    thermodynamic_forces[4] = stress[4];
    thermodynamic_forces[5] = stress[5];
    if (thermodynamic_forces.size() > 6u)
    {
        thermodynamic_forces[6] = saturation;
    }
}

struct NativeMCCOneStepResult
{
    KV stress = KV::Zero();
    KM tangent = KM::Zero();
};

std::optional<NativeMCCOneStepResult> runMCCBridgeNativeOneStep(
    double const epsilon_v,
    double const liquid_pressure)
{
    auto parameters = createMCCBridgeParameters();
    auto native = createMCCBridgeNativeModel(parameters);
    auto state = native->createMaterialStateVariables();
    initializeMCCBridgeState(*native, *state);

    MPL::VariableArray previous;
    previous.mechanical_strain.emplace<KV>(KV::Zero());
    previous.stress.emplace<KV>(mccBridgeIsotropicStress(-5e3));
    previous.liquid_phase_pressure = -5e3;
    previous.temperature = 293.15;

    MPL::VariableArray current = previous;
    current.mechanical_strain.emplace<KV>(
        mccBridgeIsotropicStrainFromVolumetric(epsilon_v));
    current.liquid_phase_pressure = liquid_pressure;

    ParameterLib::SpatialPosition x{};
    auto const solution =
        native->integrateStress(previous, current, 1.0, x, 1.0, *state);
    if (!solution)
    {
        return std::nullopt;
    }

    auto const& [stress, state_new, tangent] = *solution;
    (void)state_new;
    return NativeMCCOneStepResult{stress, tangent};
}

struct PressureCoupledMCCOneStepResult
{
    KV stress = KV::Zero();
    KM dStress_dStrain = KM::Zero();
    KV dStress_dLiquidPressure = KV::Zero();
};

std::optional<PressureCoupledMCCOneStepResult> runMCCBridgePressureCoupledOneStep(
    double const epsilon_v,
    double const liquid_pressure)
{
    auto parameters = createMCCBridgeParameters();
    auto bridge = createMCCBridgePressureCoupledModel(parameters);
    auto state = bridge->createMaterialStateVariables();
    initializeMCCBridgeState(*bridge, *state);

    MPL::VariableArray previous;
    previous.mechanical_strain.emplace<KV>(KV::Zero());
    previous.stress.emplace<KV>(mccBridgeIsotropicStress(-5e3));
    previous.liquid_phase_pressure = -5e3;
    previous.temperature = 293.15;

    MPL::VariableArray current = previous;
    current.mechanical_strain.emplace<KV>(
        mccBridgeIsotropicStrainFromVolumetric(epsilon_v));
    current.liquid_phase_pressure = liquid_pressure;

    ParameterLib::SpatialPosition x{};
    auto const response = bridge->integrateStressPressureCoupled(
        previous, current, 1.0, x, 1.0, *state);
    if (!response)
    {
        return std::nullopt;
    }

    return PressureCoupledMCCOneStepResult{response->stress,
                                           response->dStress_dStrain,
                                           response->dStress_dLiquidPressure};
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

TEST(MaterialLib_RMBridgeMFront_MCC,
     PressureCoupledConstEMatchesNativeAcrossDryDensitySweep)
{
    struct DryDensityCase
    {
        double dry_density;
    };

    std::array<DryDensityCase, 4> const dry_density_cases{{
        {1350.0},
        {1450.0},
        {1550.0},
        {1650.0},
    }};

    ParameterLib::SpatialPosition x{};
    double const residual_liquid_saturation = 0.0;
    double const residual_gas_saturation = 0.0;
    double const exponent_m = 0.4;
    double const bubble_pressure = 1e4;

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
    };

    for (auto const& c : dry_density_cases)
    {
        SCOPED_TRACE(::testing::Message()
                     << "dry_density=" << c.dry_density << " kg/m^3");

        auto parameters = createMCCBridgeParameters(
            mccBridgeInitialVolumeRatioFromDryDensity(c.dry_density));
        auto native = createMCCBridgeNativeModel(parameters);
        auto bridge = createMCCBridgePressureCoupledModel(parameters);

        auto native_state = native->createMaterialStateVariables();
        auto bridge_state = bridge->createMaterialStateVariables();
        initializeMCCBridgeState(*native, *native_state);
        initializeMCCBridgeState(*bridge, *bridge_state);

        MPL::VariableArray native_prev;
        native_prev.mechanical_strain.emplace<KV>(KV::Zero());
        native_prev.stress.emplace<KV>(mccBridgeIsotropicStress(-5e3));
        native_prev.liquid_phase_pressure = -5e3;
        native_prev.temperature = 293.15;

        MPL::VariableArray bridge_prev = native_prev;

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

            double const dt =
                i == 0 ? steps[i].t : steps[i].t - steps[i - 1].t;

            auto native_solution = native->integrateStress(
                native_prev, native_current, steps[i].t, x, dt, *native_state);
            auto bridge_solution = bridge->integrateStressPressureCoupled(
                bridge_prev, bridge_current, steps[i].t, x, dt, *bridge_state);

            ASSERT_TRUE(native_solution);
            ASSERT_TRUE(bridge_solution);

            auto& [native_stress, native_state_new, native_tangent] =
                *native_solution;
            auto& bridge_response = *bridge_solution;

            EXPECT_TRUE((native_stress - bridge_response.stress).isZero(5e-11));
            EXPECT_TRUE((native_tangent - bridge_response.dStress_dStrain)
                            .isZero(1e-8));
            EXPECT_TRUE(bridge_response.dStress_dLiquidPressure.isZero(1e-12));

            auto const expected_saturation = mccBridgeSaturationVG(
                -steps[i].liquid_pressure, residual_liquid_saturation,
                residual_gas_saturation, exponent_m, bubble_pressure);
            auto const expected_dS_dp = mccBridgeDSaturationVGdLiquidPressure(
                steps[i].liquid_pressure, residual_liquid_saturation,
                residual_gas_saturation, exponent_m, bubble_pressure);

            EXPECT_NEAR(expected_saturation, bridge_response.saturation, 1e-15);
            EXPECT_NEAR(expected_dS_dp,
                        bridge_response.dSaturation_dLiquidPressure, 1e-15);

            for (auto const& name :
                 {"EquivalentPlasticStrain", "PreConsolidationPressure",
                  "PlasticVolumetricStrain", "VolumeRatio"})
            {
                double const native_value =
                    getMCCBridgeInternalScalar(*native, *native_state_new, name);
                double const bridge_value = getMCCBridgeInternalScalar(
                    *bridge, *bridge_response.state, name);
                EXPECT_NEAR(native_value, bridge_value,
                            name == std::string("PreConsolidationPressure")
                                ? 1e-8
                                : 1e-14)
                    << name;
            }

            native_state = std::move(native_state_new);
            bridge_state = std::move(bridge_response.state);
            native_prev = native_current;
            native_prev.stress.emplace<KV>(native_stress);
            bridge_prev = bridge_current;
            bridge_prev.stress.emplace<KV>(bridge_response.stress);
        }
    }
}

TEST(MaterialLib_RMBridgeMFront_MCC,
     PressureCoupledConstETangentConsistentNearPlasticOnset)
{
    constexpr double eps_v = -1.5e-3;
    constexpr double liquid_pressure = -4e4;
    constexpr double deps_v = 1e-7;

    auto const native_base = runMCCBridgeNativeOneStep(eps_v, liquid_pressure);
    auto const native_plus =
        runMCCBridgeNativeOneStep(eps_v + deps_v, liquid_pressure);
    auto const native_minus =
        runMCCBridgeNativeOneStep(eps_v - deps_v, liquid_pressure);
    auto const bridge_base =
        runMCCBridgePressureCoupledOneStep(eps_v, liquid_pressure);
    auto const bridge_plus =
        runMCCBridgePressureCoupledOneStep(eps_v + deps_v, liquid_pressure);
    auto const bridge_minus =
        runMCCBridgePressureCoupledOneStep(eps_v - deps_v, liquid_pressure);

    ASSERT_TRUE(native_base);
    ASSERT_TRUE(native_plus);
    ASSERT_TRUE(native_minus);
    ASSERT_TRUE(bridge_base);
    ASSERT_TRUE(bridge_plus);
    ASSERT_TRUE(bridge_minus);

    auto const dSigma_dEps_v_native_fd =
        (native_plus->stress - native_minus->stress) / (2.0 * deps_v);
    auto const dSigma_dEps_v_bridge_fd =
        (bridge_plus->stress - bridge_minus->stress) / (2.0 * deps_v);

    auto const dEps_dEps_v = mccBridgeIsotropicStrainFromVolumetric(1.0);
    auto const dSigma_dEps_v_native_tangent =
        native_base->tangent * dEps_dEps_v;
    auto const dSigma_dEps_v_bridge_tangent =
        bridge_base->dStress_dStrain * dEps_dEps_v;

    for (Eigen::Index i = 0; i < KV::RowsAtCompileTime; ++i)
    {
        auto const native_fd = dSigma_dEps_v_native_fd[i];
        auto const native_tangent = dSigma_dEps_v_native_tangent[i];
        auto const bridge_fd = dSigma_dEps_v_bridge_fd[i];
        auto const bridge_tangent = dSigma_dEps_v_bridge_tangent[i];

        double const native_tol =
            1e2 + 5e-3 * std::max(std::abs(native_fd), std::abs(native_tangent));
        double const bridge_tol =
            1e2 + 5e-3 * std::max(std::abs(bridge_fd), std::abs(bridge_tangent));
        double const parity_tol = 1e2 + 1e-4 * std::max(
                                            std::abs(native_tangent),
                                            std::abs(bridge_tangent));

        EXPECT_NEAR(native_fd, native_tangent, native_tol);
        EXPECT_NEAR(bridge_fd, bridge_tangent, bridge_tol);
        EXPECT_NEAR(native_tangent, bridge_tangent, parity_tol);
    }

    EXPECT_TRUE((native_base->stress - bridge_base->stress).isZero(5e-11));
    EXPECT_TRUE((native_base->tangent - bridge_base->dStress_dStrain)
                    .isZero(1e-8));
    EXPECT_TRUE(bridge_base->dStress_dLiquidPressure.isZero(1e-12));
}

TEST(MaterialLib_RMBridgeMFront_MCC,
     PressureCoupledConstERestartStateTransferMatchesContinuousPath)
{
    auto parameters = createMCCBridgeParameters();
    auto native = createMCCBridgeNativeModel(parameters);
    auto bridge = createMCCBridgePressureCoupledModel(parameters);

    auto native_state = native->createMaterialStateVariables();
    auto bridge_state = bridge->createMaterialStateVariables();
    initializeMCCBridgeState(*native, *native_state);
    initializeMCCBridgeState(*bridge, *bridge_state);

    ParameterLib::SpatialPosition x{};

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
    };

    MPL::VariableArray native_prev;
    native_prev.mechanical_strain.emplace<KV>(KV::Zero());
    native_prev.stress.emplace<KV>(mccBridgeIsotropicStress(-5e3));
    native_prev.liquid_phase_pressure = -5e3;
    native_prev.temperature = 293.15;

    MPL::VariableArray bridge_prev = native_prev;

    // Advance to a checkpoint (step index 1).
    for (std::size_t i = 0; i < 2; ++i)
    {
        MPL::VariableArray native_current = native_prev;
        native_current.mechanical_strain.emplace<KV>(
            mccBridgeIsotropicStrainFromVolumetric(steps[i].eps_v));
        native_current.liquid_phase_pressure = steps[i].liquid_pressure;

        MPL::VariableArray bridge_current = bridge_prev;
        bridge_current.mechanical_strain.emplace<KV>(
            mccBridgeIsotropicStrainFromVolumetric(steps[i].eps_v));
        bridge_current.liquid_phase_pressure = steps[i].liquid_pressure;

        double const dt =
            i == 0 ? steps[i].t : steps[i].t - steps[i - 1].t;

        auto native_solution = native->integrateStress(
            native_prev, native_current, steps[i].t, x, dt, *native_state);
        auto bridge_solution = bridge->integrateStressPressureCoupled(
            bridge_prev, bridge_current, steps[i].t, x, dt, *bridge_state);

        ASSERT_TRUE(native_solution);
        ASSERT_TRUE(bridge_solution);

        auto& [native_stress, native_state_new, native_tangent] =
            *native_solution;
        (void)native_tangent;
        auto& bridge_response = *bridge_solution;

        native_state = std::move(native_state_new);
        bridge_state = std::move(bridge_response.state);

        native_prev = native_current;
        native_prev.stress.emplace<KV>(native_stress);
        bridge_prev = bridge_current;
        bridge_prev.stress.emplace<KV>(bridge_response.stress);
    }

    // Continuous continuation for step index 2.
    MPL::VariableArray native_current_cont = native_prev;
    native_current_cont.mechanical_strain.emplace<KV>(
        mccBridgeIsotropicStrainFromVolumetric(steps[2].eps_v));
    native_current_cont.liquid_phase_pressure = steps[2].liquid_pressure;

    MPL::VariableArray bridge_current_cont = bridge_prev;
    bridge_current_cont.mechanical_strain.emplace<KV>(
        mccBridgeIsotropicStrainFromVolumetric(steps[2].eps_v));
    bridge_current_cont.liquid_phase_pressure = steps[2].liquid_pressure;

    double const dt_cont = steps[2].t - steps[1].t;

    auto native_continuous = native->integrateStress(
        native_prev, native_current_cont, steps[2].t, x, dt_cont, *native_state);
    auto bridge_continuous = bridge->integrateStressPressureCoupled(
        bridge_prev, bridge_current_cont, steps[2].t, x, dt_cont, *bridge_state);
    ASSERT_TRUE(native_continuous);
    ASSERT_TRUE(bridge_continuous);

    // Restart from checkpointed internal variables.
    auto restart_native = createMCCBridgeNativeModel(parameters);
    auto restart_bridge = createMCCBridgePressureCoupledModel(parameters);
    auto restart_native_state = restart_native->createMaterialStateVariables();
    auto restart_bridge_state = restart_bridge->createMaterialStateVariables();
    initializeMCCBridgeState(*restart_native, *restart_native_state);
    initializeMCCBridgeState(*restart_bridge, *restart_bridge_state);

    ASSERT_TRUE(copyMCCBridgeInternalVariables(
        *native, *native_state, *restart_native, *restart_native_state));
    ASSERT_TRUE(copyMCCBridgeInternalVariables(
        *bridge, *bridge_state, *restart_bridge, *restart_bridge_state));

    auto const checkpoint_native_stress = std::get<KV>(native_prev.stress);
    auto const checkpoint_bridge_stress = std::get<KV>(bridge_prev.stress);
    double const checkpoint_saturation = mccBridgeSaturationVG(
        -bridge_prev.liquid_phase_pressure, 0.0, 0.0, 0.4, 1e4);

    setMCCBridgeThermodynamicForces(*restart_native_state,
                                    checkpoint_native_stress,
                                    checkpoint_saturation);
    setMCCBridgeThermodynamicForces(*restart_bridge_state,
                                    checkpoint_bridge_stress,
                                    checkpoint_saturation);
    restart_native_state->pushBackState();
    restart_bridge_state->pushBackState();

    MPL::VariableArray native_prev_restart = native_prev;
    MPL::VariableArray bridge_prev_restart = bridge_prev;

    MPL::VariableArray native_current_restart = native_prev_restart;
    native_current_restart.mechanical_strain.emplace<KV>(
        mccBridgeIsotropicStrainFromVolumetric(steps[2].eps_v));
    native_current_restart.liquid_phase_pressure = steps[2].liquid_pressure;

    MPL::VariableArray bridge_current_restart = bridge_prev_restart;
    bridge_current_restart.mechanical_strain.emplace<KV>(
        mccBridgeIsotropicStrainFromVolumetric(steps[2].eps_v));
    bridge_current_restart.liquid_phase_pressure = steps[2].liquid_pressure;

    auto native_restarted = restart_native->integrateStress(
        native_prev_restart, native_current_restart, steps[2].t, x, dt_cont,
        *restart_native_state);
    auto bridge_restarted = restart_bridge->integrateStressPressureCoupled(
        bridge_prev_restart, bridge_current_restart, steps[2].t, x, dt_cont,
        *restart_bridge_state);
    ASSERT_TRUE(native_restarted);
    ASSERT_TRUE(bridge_restarted);

    auto const& [native_continuous_stress, native_continuous_state,
                 native_continuous_tangent] = *native_continuous;
    auto const& [native_restarted_stress, native_restarted_state,
                 native_restarted_tangent] = *native_restarted;
    auto const& bridge_continuous_response = *bridge_continuous;
    auto const& bridge_restarted_response = *bridge_restarted;

    EXPECT_TRUE((native_continuous_stress - native_restarted_stress)
                    .isZero(1e-8));
    EXPECT_TRUE((native_continuous_tangent - native_restarted_tangent)
                    .isZero(1e-7));
    EXPECT_NEAR(getMCCBridgeInternalScalar(*native, *native_continuous_state,
                                           "EquivalentPlasticStrain"),
                getMCCBridgeInternalScalar(*restart_native,
                                           *native_restarted_state,
                                           "EquivalentPlasticStrain"),
                1e-12);
    EXPECT_NEAR(getMCCBridgeInternalScalar(*native, *native_continuous_state,
                                           "VolumeRatio"),
                getMCCBridgeInternalScalar(*restart_native,
                                           *native_restarted_state,
                                           "VolumeRatio"),
                2e-3);

    EXPECT_TRUE((bridge_continuous_response.stress -
                 bridge_restarted_response.stress)
                    .isZero(1e-8));
    EXPECT_TRUE((bridge_continuous_response.dStress_dStrain -
                 bridge_restarted_response.dStress_dStrain)
                    .isZero(1e-7));
    EXPECT_TRUE((bridge_continuous_response.dStress_dLiquidPressure -
                 bridge_restarted_response.dStress_dLiquidPressure)
                    .isZero(1e-9));
    EXPECT_NEAR(bridge_continuous_response.saturation,
                bridge_restarted_response.saturation, 1e-14);
    EXPECT_NEAR(bridge_continuous_response.dSaturation_dLiquidPressure,
                bridge_restarted_response.dSaturation_dLiquidPressure, 1e-14);

    EXPECT_TRUE((native_restarted_stress - bridge_restarted_response.stress)
                    .isZero(5e-10));
    EXPECT_TRUE((native_restarted_tangent -
                 bridge_restarted_response.dStress_dStrain)
                    .isZero(1e-7));
}

#endif  // OGS_USE_MFRONT
