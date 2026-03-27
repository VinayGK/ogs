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
using MB = MaterialLib::Solids::MechanicsBase<3>;

namespace
{
std::vector<std::unique_ptr<ParameterLib::ParameterBase>>
createNotebookMCCParameters(double const swelling_slope = 0.0,
                            double const mass_exchange_coefficient = 0.0,
                            double const n_l0 = 0.1,
                            double const rho_lR0 = 1300.0,
                            double const epsilon_sw0 = 0.0,
                            double const notebook_saturation_mode = 0.0)
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
    add_param("NotebookSaturationMode", notebook_saturation_mode);

    // Neutral first-step notebook coupling: state is updated and visible, but
    // the verified MCC stress/saturation surface stays unchanged.
    add_param("SwellingSlope", swelling_slope);
    add_param("MassExchangeCoefficient", mass_exchange_coefficient);
    add_param("ReferenceLiquidDensityMacro", 1000.0);
    add_param("ReferenceLiquidDensityMicro", 1300.0);
    add_param("ReferenceDensitySolid", 2470.0);
    add_param("MicroLiquidDensityA", 1.3);
    add_param("MicroLiquidDensityB", 1.0);
    add_param("HamakerConstant", -6e-20);
    add_param("SpecificSurface", 100.0);
    add_param("AreaFactorTuller", 1.0);
    add_param("PoreAreaShapeFactorTuller", 0.8584073464102069);
    add_param("CharacteristicPoreSize", 1e-5);
    add_param("SurfaceTension", 0.0715);
    add_param("InitialPorosity", 0.432);
    add_param("n_l0", n_l0);
    add_param("rho_lR0", rho_lR0);
    add_param("epsilon_sw0", epsilon_sw0);

    return parameters;
}

std::unique_ptr<MB> createModelFromXml(
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

std::unique_ptr<MB> createReferenceMCCModel(
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

    return createModelFromXml(parameters, xml);
}

std::unique_ptr<MB> createNotebookMCCModel(
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters)
{
    char const* xml = R"XML(
        <type>MFrontRichardsMechanics</type>
        <behaviour>RichardsMechanicsNotebookBridge_MCC</behaviour>
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
            <material_property name="NotebookSaturationMode" parameter="NotebookSaturationMode"/>
            <material_property name="SwellingSlope" parameter="SwellingSlope"/>
            <material_property name="MassExchangeCoefficient" parameter="MassExchangeCoefficient"/>
            <material_property name="ReferenceLiquidDensityMacro" parameter="ReferenceLiquidDensityMacro"/>
            <material_property name="ReferenceLiquidDensityMicro" parameter="ReferenceLiquidDensityMicro"/>
            <material_property name="ReferenceDensitySolid" parameter="ReferenceDensitySolid"/>
            <material_property name="MicroLiquidDensityA" parameter="MicroLiquidDensityA"/>
            <material_property name="MicroLiquidDensityB" parameter="MicroLiquidDensityB"/>
            <material_property name="HamakerConstant" parameter="HamakerConstant"/>
            <material_property name="SpecificSurface" parameter="SpecificSurface"/>
            <material_property name="AreaFactorTuller" parameter="AreaFactorTuller"/>
            <material_property name="PoreAreaShapeFactorTuller" parameter="PoreAreaShapeFactorTuller"/>
            <material_property name="CharacteristicPoreSize" parameter="CharacteristicPoreSize"/>
            <material_property name="SurfaceTension" parameter="SurfaceTension"/>
            <material_property name="InitialPorosity" parameter="InitialPorosity"/>
        </material_properties>
        <initial_values>
            <state_variable name="PreConsolidationPressure" parameter="InitialPreConsolidationPressure"/>
            <state_variable name="VolumeRatio" parameter="InitialVolumeRatio"/>
            <state_variable name="n_l" parameter="n_l0"/>
            <state_variable name="rho_lR" parameter="rho_lR0"/>
            <state_variable name="epsilon_sw" parameter="epsilon_sw0"/>
        </initial_values>
        )XML";

    return createModelFromXml(parameters, xml);
}

void initializeState(MB const& model, MB::MaterialStateVariables& state)
{
    ParameterLib::SpatialPosition x{};
    model.initializeInternalStateVariables(0.0, x, state);
}

double getInternalScalar(MB const& model,
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

std::vector<double> getInternalVector(MB const& model,
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

KV notebookMCCIsotropicStress(double const value)
{
    KV sigma = KV::Zero();
    sigma[0] = value;
    sigma[1] = value;
    sigma[2] = value;
    return sigma;
}

KV notebookMCCIsotropicStrainFromVolumetric(double const epsilon_v)
{
    KV eps = KV::Zero();
    eps[0] = epsilon_v / 3.0;
    eps[1] = epsilon_v / 3.0;
    eps[2] = epsilon_v / 3.0;
    return eps;
}

double notebookTullerSaturation(double const liquid_pressure)
{
    if (liquid_pressure >= -1e-12)
    {
        return 1.0;
    }

    constexpr double area_factor_tuller = 1.0;
    constexpr double pore_area_shape_factor_tuller = 0.8584073464102069;
    constexpr double characteristic_pore_size = 1e-5;
    constexpr double surface_tension = 0.0715;

    double const prefactor =
        4.0 * pore_area_shape_factor_tuller * surface_tension * surface_tension /
        (area_factor_tuller * characteristic_pore_size *
         characteristic_pore_size);
    return 1.0 - std::exp(-prefactor / (liquid_pressure * liquid_pressure));
}

double notebookTullerDSaturationDLiquidPressure(double const liquid_pressure)
{
    if (liquid_pressure >= -1e-12)
    {
        return 0.0;
    }

    constexpr double area_factor_tuller = 1.0;
    constexpr double pore_area_shape_factor_tuller = 0.8584073464102069;
    constexpr double characteristic_pore_size = 1e-5;
    constexpr double surface_tension = 0.0715;

    double const prefactor =
        4.0 * pore_area_shape_factor_tuller * surface_tension * surface_tension /
        (area_factor_tuller * characteristic_pore_size *
         characteristic_pore_size);
    double const inv_p = 1.0 / liquid_pressure;
    double const exp_term = std::exp(-prefactor * inv_p * inv_p);
    return -2.0 * prefactor * exp_term * inv_p * inv_p * inv_p;
}
}  // namespace

TEST(MaterialLib_RMBridgeMFront_NotebookMCC,
     NeutralNotebookStateMatchesVerifiedMCCBridge)
{
    auto parameters = createNotebookMCCParameters();
    auto reference = createReferenceMCCModel(parameters);
    auto notebook_mcc = createNotebookMCCModel(parameters);

    auto reference_state = reference->createMaterialStateVariables();
    auto notebook_state = notebook_mcc->createMaterialStateVariables();
    initializeState(*reference, *reference_state);
    initializeState(*notebook_mcc, *notebook_state);

    ParameterLib::SpatialPosition x{};

    MPL::VariableArray previous;
    previous.mechanical_strain.emplace<KV>(KV::Zero());
    previous.stress.emplace<KV>(notebookMCCIsotropicStress(-5e3));
    previous.liquid_phase_pressure = -5e3;
    previous.temperature = 293.15;

    auto reference_previous = previous;
    auto notebook_previous = previous;

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
        MPL::VariableArray reference_current = reference_previous;
        reference_current.mechanical_strain.emplace<KV>(
            notebookMCCIsotropicStrainFromVolumetric(steps[i].eps_v));
        reference_current.liquid_phase_pressure = steps[i].liquid_pressure;

        MPL::VariableArray notebook_current = notebook_previous;
        notebook_current.mechanical_strain.emplace<KV>(
            notebookMCCIsotropicStrainFromVolumetric(steps[i].eps_v));
        notebook_current.liquid_phase_pressure = steps[i].liquid_pressure;

        double const dt = i == 0 ? steps[i].t : steps[i].t - steps[i - 1].t;

        auto reference_response = reference->integrateStressPressureCoupled(
            reference_previous, reference_current, steps[i].t, x, dt,
            *reference_state);
        auto notebook_response = notebook_mcc->integrateStressPressureCoupled(
            notebook_previous, notebook_current, steps[i].t, x, dt,
            *notebook_state);

        ASSERT_TRUE(reference_response);
        ASSERT_TRUE(notebook_response);
        ASSERT_TRUE(reference_response->state);
        ASSERT_TRUE(notebook_response->state);

        EXPECT_TRUE((reference_response->stress - notebook_response->stress)
                        .isZero(5e-11));
        EXPECT_TRUE(
            (reference_response->dStress_dStrain -
             notebook_response->dStress_dStrain)
                .isZero(1e-8));
        EXPECT_TRUE((reference_response->dStress_dLiquidPressure -
                     notebook_response->dStress_dLiquidPressure)
                        .isZero(1e-12));
        EXPECT_TRUE((reference_response->dSaturation_dStrain -
                     notebook_response->dSaturation_dStrain)
                        .isZero(1e-12));
        EXPECT_NEAR(reference_response->saturation, notebook_response->saturation,
                    1e-15);
        EXPECT_NEAR(reference_response->dSaturation_dLiquidPressure,
                    notebook_response->dSaturation_dLiquidPressure, 1e-15);

        for (auto const& name :
             {"EquivalentPlasticStrain", "PreConsolidationPressure",
              "PlasticVolumetricStrain", "VolumeRatio"})
        {
            auto const reference_value =
                getInternalScalar(*reference, *reference_response->state, name);
            auto const notebook_value =
                getInternalScalar(*notebook_mcc, *notebook_response->state,
                                  name);
            EXPECT_NEAR(reference_value, notebook_value,
                        name == std::string("PreConsolidationPressure") ? 1e-8
                                                                        : 1e-14)
                << name;
        }

        auto const reference_elastic_strain = getInternalVector(
            *reference, *reference_response->state, "ElasticStrain");
        auto const notebook_elastic_strain = getInternalVector(
            *notebook_mcc, *notebook_response->state, "ElasticStrain");
        ASSERT_EQ(reference_elastic_strain.size(), notebook_elastic_strain.size());
        for (std::size_t k = 0; k < reference_elastic_strain.size(); ++k)
        {
            EXPECT_NEAR(reference_elastic_strain[k], notebook_elastic_strain[k],
                        1e-14);
        }

        auto const n_l =
            getInternalScalar(*notebook_mcc, *notebook_response->state, "n_l");
        auto const rho_lR =
            getInternalScalar(*notebook_mcc, *notebook_response->state,
                              "rho_lR");
        auto const epsilon_sw = getInternalScalar(*notebook_mcc,
                                                  *notebook_response->state,
                                                  "epsilon_sw");
        auto const phi_m =
            getInternalScalar(*notebook_mcc, *notebook_response->state, "phi_m");
        auto const phi_M =
            getInternalScalar(*notebook_mcc, *notebook_response->state, "phi_M");
        auto const mu_lR = getInternalScalar(*notebook_mcc,
                                             *notebook_response->state,
                                             "mu_lR");
        auto const rho_l_hat =
            getInternalScalar(*notebook_mcc, *notebook_response->state,
                              "rho_l_hat");

        EXPECT_TRUE(std::isfinite(n_l));
        EXPECT_TRUE(std::isfinite(rho_lR));
        EXPECT_TRUE(std::isfinite(epsilon_sw));
        EXPECT_TRUE(std::isfinite(phi_m));
        EXPECT_TRUE(std::isfinite(phi_M));
        EXPECT_TRUE(std::isfinite(mu_lR));
        EXPECT_TRUE(std::isfinite(rho_l_hat));
        EXPECT_GT(n_l, 0.0);
        EXPECT_GT(rho_lR, 1000.0);
        EXPECT_NEAR(phi_m, n_l, 1e-14);
        EXPECT_NEAR(phi_M, 0.432 - n_l, 1e-14);
        EXPECT_NEAR(epsilon_sw, 0.0, 1e-14);
        EXPECT_NEAR(rho_l_hat, 0.0, 1e-14);

        reference_state = std::move(reference_response->state);
        notebook_state = std::move(notebook_response->state);

        reference_previous = reference_current;
        reference_previous.stress.emplace<KV>(reference_response->stress);
        notebook_previous = notebook_current;
        notebook_previous.stress.emplace<KV>(notebook_response->stress);
    }
}

TEST(MaterialLib_RMBridgeMFront_NotebookMCC,
     SwellingFeedbackChangesStressButKeepsCarrierSaturation)
{
    auto parameters = createNotebookMCCParameters(2.0, 1e-4, 0.1, 1300.0, 0.0);
    auto reference = createReferenceMCCModel(parameters);
    auto notebook_mcc = createNotebookMCCModel(parameters);

    auto reference_state = reference->createMaterialStateVariables();
    auto notebook_state = notebook_mcc->createMaterialStateVariables();
    initializeState(*reference, *reference_state);
    initializeState(*notebook_mcc, *notebook_state);

    ParameterLib::SpatialPosition x{};

    MPL::VariableArray previous;
    previous.mechanical_strain.emplace<KV>(KV::Zero());
    previous.stress.emplace<KV>(notebookMCCIsotropicStress(-5e3));
    previous.liquid_phase_pressure = -5e3;
    previous.temperature = 293.15;

    constexpr double t = 1.0;
    constexpr double dt = 1.0;
    constexpr double eps_v = -2e-5;
    constexpr double liquid_pressure = -4e4;

    MPL::VariableArray current = previous;
    current.mechanical_strain.emplace<KV>(
        notebookMCCIsotropicStrainFromVolumetric(eps_v));
    current.liquid_phase_pressure = liquid_pressure;

    auto reference_response = reference->integrateStressPressureCoupled(
        previous, current, t, x, dt, *reference_state);
    auto notebook_response = notebook_mcc->integrateStressPressureCoupled(
        previous, current, t, x, dt, *notebook_state);

    ASSERT_TRUE(reference_response);
    ASSERT_TRUE(notebook_response);
    ASSERT_TRUE(reference_response->state);
    ASSERT_TRUE(notebook_response->state);

    auto const epsilon_sw =
        getInternalScalar(*notebook_mcc, *notebook_response->state,
                          "epsilon_sw");
    EXPECT_GT(std::abs(epsilon_sw), 1e-12);

    double const bulk_modulus = 52e6 / (3.0 * (1.0 - 2.0 * 0.3));
    auto const expected_stress =
        reference_response->stress -
        notebookMCCIsotropicStress(bulk_modulus * epsilon_sw);
    EXPECT_TRUE((expected_stress - notebook_response->stress).isZero(5e-8));

    EXPECT_NEAR(reference_response->saturation, notebook_response->saturation,
                1e-15);
    EXPECT_TRUE((reference_response->dSaturation_dStrain -
                 notebook_response->dSaturation_dStrain)
                    .isZero(1e-12));
    EXPECT_NEAR(reference_response->dSaturation_dLiquidPressure,
                notebook_response->dSaturation_dLiquidPressure, 1e-15);

    auto evaluate_notebook_response =
        [&](double const eps_v_value, double const pressure_value)
    {
        auto local_parameters =
            createNotebookMCCParameters(2.0, 1e-4, 0.1, 1300.0, 0.0);
        auto local_model = createNotebookMCCModel(local_parameters);
        auto local_state = local_model->createMaterialStateVariables();
        initializeState(*local_model, *local_state);

        MPL::VariableArray local_current = previous;
        local_current.mechanical_strain.emplace<KV>(
            notebookMCCIsotropicStrainFromVolumetric(eps_v_value));
        local_current.liquid_phase_pressure = pressure_value;

        return local_model->integrateStressPressureCoupled(
            previous, local_current, t, x, dt, *local_state);
    };

    double const dp_fd = std::max(1e-6, std::abs(liquid_pressure) * 1e-6);
    auto notebook_p_plus = evaluate_notebook_response(eps_v, liquid_pressure + dp_fd);
    auto notebook_p_minus = evaluate_notebook_response(eps_v, liquid_pressure - dp_fd);
    ASSERT_TRUE(notebook_p_plus);
    ASSERT_TRUE(notebook_p_minus);
    auto const dsigma_dp_fd =
        (notebook_p_plus->stress - notebook_p_minus->stress) / (2.0 * dp_fd);
    EXPECT_TRUE((notebook_response->dStress_dLiquidPressure - dsigma_dp_fd)
                    .isZero(5e-4));
}

TEST(MaterialLib_RMBridgeMFront_NotebookMCC,
     NotebookSaturationModeMatchesTullerLawAndKeepsStressSurface)
{
    auto parameters = createNotebookMCCParameters(0.0, 0.0, 0.1, 1300.0, 0.0,
                                                  1.0);
    auto reference = createReferenceMCCModel(parameters);
    auto notebook_mcc = createNotebookMCCModel(parameters);

    auto reference_state = reference->createMaterialStateVariables();
    auto notebook_state = notebook_mcc->createMaterialStateVariables();
    initializeState(*reference, *reference_state);
    initializeState(*notebook_mcc, *notebook_state);

    ParameterLib::SpatialPosition x{};

    MPL::VariableArray previous;
    previous.mechanical_strain.emplace<KV>(KV::Zero());
    previous.stress.emplace<KV>(notebookMCCIsotropicStress(-5e3));
    previous.liquid_phase_pressure = -5e3;
    previous.temperature = 293.15;

    constexpr double t = 1.0;
    constexpr double dt = 1.0;
    constexpr double eps_v = -2e-5;
    constexpr double liquid_pressure = -4e4;

    MPL::VariableArray current = previous;
    current.mechanical_strain.emplace<KV>(
        notebookMCCIsotropicStrainFromVolumetric(eps_v));
    current.liquid_phase_pressure = liquid_pressure;

    auto reference_response = reference->integrateStressPressureCoupled(
        previous, current, t, x, dt, *reference_state);
    auto notebook_response = notebook_mcc->integrateStressPressureCoupled(
        previous, current, t, x, dt, *notebook_state);

    ASSERT_TRUE(reference_response);
    ASSERT_TRUE(notebook_response);
    ASSERT_TRUE(reference_response->state);
    ASSERT_TRUE(notebook_response->state);

    EXPECT_TRUE((reference_response->stress - notebook_response->stress)
                    .isZero(5e-11));
    EXPECT_TRUE((reference_response->dStress_dStrain -
                 notebook_response->dStress_dStrain)
                    .isZero(1e-8));
    EXPECT_TRUE((reference_response->dStress_dLiquidPressure -
                 notebook_response->dStress_dLiquidPressure)
                    .isZero(1e-12));

    EXPECT_NEAR(notebook_response->saturation,
                notebookTullerSaturation(liquid_pressure), 1e-15);
    EXPECT_NEAR(notebook_response->dSaturation_dLiquidPressure,
                notebookTullerDSaturationDLiquidPressure(liquid_pressure),
                1e-15);
    EXPECT_TRUE(notebook_response->dSaturation_dStrain.isZero(1e-12));

    auto const epsilon_sw =
        getInternalScalar(*notebook_mcc, *notebook_response->state,
                          "epsilon_sw");
    EXPECT_NEAR(epsilon_sw, 0.0, 1e-14);

    MPL::VariableArray saturated_current = previous;
    saturated_current.mechanical_strain.emplace<KV>(
        notebookMCCIsotropicStrainFromVolumetric(eps_v));
    saturated_current.liquid_phase_pressure = 1e3;

    auto saturated_response = notebook_mcc->integrateStressPressureCoupled(
        previous, saturated_current, t, x, dt, *notebook_state);
    ASSERT_TRUE(saturated_response);
    EXPECT_NEAR(saturated_response->saturation, 1.0, 1e-15);
    EXPECT_NEAR(saturated_response->dSaturation_dLiquidPressure, 0.0, 1e-15);
    EXPECT_TRUE(saturated_response->dSaturation_dStrain.isZero(1e-12));
}

#endif  // OGS_USE_MFRONT
