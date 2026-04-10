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

#include "BaseLib/ConfigTree.h"
#include "InfoLib/TestInfo.h"
#include "MaterialLib/SolidModels/CreateConstitutiveRelation.h"
#include "MaterialLib/SolidModels/MFront/MFrontGeneric.h"
#include "ParameterLib/ConstantParameter.h"
#include "Tests/TestTools.h"

namespace MPL = MaterialPropertyLib;
namespace MSM = MaterialLib::Solids::MFront;
using KV = MathLib::KelvinVector::KelvinVectorType<3>;
using MB = MaterialLib::Solids::MechanicsBase<3>;

namespace
{
std::vector<std::unique_ptr<ParameterLib::ParameterBase>>
createDSMMicroMacroMCCParameters(double const swelling_slope = 0.0,
                            double const mass_exchange_coefficient = 0.0,
                            double const n_l0 = 0.1,
                            double const rho_lR0 = 1300.0,
                            double const epsilon_sw0 = 0.0,
                            double const macro_saturation_law_mode = 0.0,
                            double const microstate_local_solve_mode = 0.0,
                            double const micro_potential_convention = 0.0,
                            double const hamaker_constant = -6e-20,
                            double const macro_viscosity = 1000.0)
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
    add_param("MacroSaturationLawMode", macro_saturation_law_mode);
    add_param("MicroStateLocalSolveMode", microstate_local_solve_mode);
    add_param("MicroPotentialConvention", micro_potential_convention);

    // Neutral first-step dsm_micromacro coupling: state is updated and visible, but
    // the verified MCC stress/saturation surface stays unchanged.
    add_param("SwellingSlope", swelling_slope);
    add_param("MassExchangeCoefficient", mass_exchange_coefficient);
    add_param("MacroViscosity", macro_viscosity);
    add_param("ReferenceLiquidDensityMacro", 1000.0);
    add_param("ReferenceLiquidDensityMicro", 1300.0);
    add_param("ReferenceDensitySolid", 2470.0);
    add_param("MicroLiquidDensityA", 1.3);
    add_param("MicroLiquidDensityB", 1.0);
    add_param("HamakerConstant", hamaker_constant);
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

std::vector<std::unique_ptr<ParameterLib::ParameterBase>>
createDSMMicroMacroSupportParameters(double const macro_saturation_law_mode = 0.0,
                                double const microstate_local_solve_mode = 0.0,
                                double const micro_potential_convention = 0.0,
                                double const hamaker_constant = -6e-20)
{
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> parameters;

    auto add_param = [&parameters](char const* name, double const value)
    {
        parameters.push_back(
            std::make_unique<ParameterLib::ConstantParameter<double>>(name,
                                                                      value));
    };

    // Use the original dsm_micromacro support-state settings for the auxiliary
    // surface, while keeping the MCC carrier mechanically stable.
    add_param("YoungModulus", 1e10);
    add_param("PoissonRatio", 0.25);
    add_param("CriticalStateLineSlope", 1.2);
    add_param("SwellingLineSlope", 6.6e-3);
    add_param("VirginConsolidationLineSlope", 7.7e-2);
    add_param("InitialPreConsolidationPressure", 1e9);
    add_param("InitialVolumeRatio", 1.7857142857142858);
    add_param("ResidualLiquidSaturation", 0.0);
    add_param("ResidualGasSaturation", 0.0);
    add_param("BubblePressure", 1e4);
    add_param("VanGenuchtenExponent_m", 0.4);
    add_param("MacroSaturationLawMode", macro_saturation_law_mode);
    add_param("MicroStateLocalSolveMode", microstate_local_solve_mode);
    add_param("MicroPotentialConvention", micro_potential_convention);

    add_param("SwellingSlope", 0.1);
    add_param("MassExchangeCoefficient", 1.0);
    add_param("MacroViscosity", 1000.0);
    add_param("ReferenceLiquidDensityMacro", 1000.0);
    add_param("ReferenceLiquidDensityMicro", 1300.0);
    add_param("ReferenceDensitySolid", 2470.0);
    add_param("MicroLiquidDensityA", 1.3);
    add_param("MicroLiquidDensityB", 1.0);
    add_param("HamakerConstant", hamaker_constant);
    add_param("SpecificSurface", 100.0);
    add_param("AreaFactorTuller", 1.0);
    add_param("PoreAreaShapeFactorTuller", 0.8584073464102069);
    add_param("CharacteristicPoreSize", 1e-5);
    add_param("SurfaceTension", 0.0715);
    add_param("InitialPorosity", 0.2);
    add_param("n_l0", 0.1);
    add_param("rho_lR0", 1300.0);
    add_param("epsilon_sw0", 0.0);

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

std::unique_ptr<MB> createDSMMicroMacroMCCModel(
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters)
{
    char const* xml = R"XML(
        <type>MFrontRichardsMechanics</type>
        <behaviour>RichardsMechanicsDSMMicroMacroBridge_MCC</behaviour>
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
            <material_property name="MacroSaturationLawMode" parameter="MacroSaturationLawMode"/>
            <material_property name="MicroStateLocalSolveMode" parameter="MicroStateLocalSolveMode"/>
            <material_property name="MicroPotentialConvention" parameter="MicroPotentialConvention"/>
            <material_property name="SwellingSlope" parameter="SwellingSlope"/>
            <material_property name="MassExchangeCoefficient" parameter="MassExchangeCoefficient"/>
            <material_property name="MacroViscosity" parameter="MacroViscosity"/>
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

void setInternalScalar(MB const& model,
                       MB::MaterialStateVariables& state,
                       std::string const& name,
                       double const value)
{
    auto const internal_variables = model.getInternalVariables();
    auto const it = std::find_if(
        internal_variables.begin(), internal_variables.end(),
        [&name](auto const& internal_variable)
        { return internal_variable.name == name; });

    ASSERT_TRUE(it != internal_variables.end()) << name;
    auto values = it->reference(state);
    ASSERT_EQ(values.size(), 1u);
    values[0] = value;
}

void setDSMMicroMacroMCCThermodynamicForces(MB::MaterialStateVariables& state,
                                       KV const& stress,
                                       double const saturation)
{
    auto* const mfront_state =
        dynamic_cast<MSM::MaterialStateVariablesMFront<3>*>(&state);
    ASSERT_TRUE(mfront_state != nullptr);

    auto& thermodynamic_forces =
        mfront_state->_behaviour_data.s1.thermodynamic_forces;
    thermodynamic_forces[0] = stress[0];
    thermodynamic_forces[1] = stress[1];
    thermodynamic_forces[2] = stress[2];
    thermodynamic_forces[3] = stress[3];
    thermodynamic_forces[4] = stress[4];
    thermodynamic_forces[5] = stress[5];
    thermodynamic_forces[6] = saturation;
}

KV dsm_micromacroMCCIsotropicStress(double const value)
{
    KV sigma = KV::Zero();
    sigma[0] = value;
    sigma[1] = value;
    sigma[2] = value;
    return sigma;
}

KV dsm_micromacroMCCIsotropicStrainFromVolumetric(double const epsilon_v)
{
    KV eps = KV::Zero();
    eps[0] = epsilon_v / 3.0;
    eps[1] = epsilon_v / 3.0;
    eps[2] = epsilon_v / 3.0;
    return eps;
}

double dsm_micromacroTullerSaturation(double const liquid_pressure)
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

double dsm_micromacroTullerDSaturationDLiquidPressure(double const liquid_pressure)
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

std::string dsm_micromacroMCCStripQuotes(std::string value)
{
    value.erase(0, value.find_first_not_of(" \t\r\n\""));
    value.erase(value.find_last_not_of(" \t\r\n\"") + 1);
    return value;
}

std::vector<std::string> dsm_micromacroMCCSplitCommaLine(std::string const& line)
{
    std::stringstream ss(line);
    std::string field;
    std::vector<std::string> fields;
    while (std::getline(ss, field, ','))
    {
        fields.push_back(dsm_micromacroMCCStripQuotes(field));
    }
    return fields;
}

struct DSMMicroMacroSupportBaselineRow
{
    int step = 0;
    double pressure = 0.0;
    double epsilon_v_total = 0.0;
    double saturation = 0.0;
    double n_l = 0.0;
    double phi_m = 0.0;
    double phi_M = 0.0;
    double phi = 0.0;
    double n_S = 0.0;
    double n_L = 0.0;
    double rho_lR = 0.0;
    double rho_LR = 0.0;
    double omega_l = 0.0;
    double mu_lR = 0.0;
    double rho_l_hat = 0.0;
    double delta_epsilon_sw = 0.0;
    double epsilon_sw = 0.0;
    double sigma_S_xx = 0.0;
};

std::vector<DSMMicroMacroSupportBaselineRow> loadDSMMicroMacroSupportBaselineRows(
    std::string const& filename)
{
    std::ifstream in(filename);
    EXPECT_TRUE(in.good()) << filename;

    std::string header_line;
    std::getline(in, header_line);
    auto const headers = dsm_micromacroMCCSplitCommaLine(header_line);

    std::unordered_map<std::string, std::size_t> column;
    for (std::size_t i = 0; i < headers.size(); ++i)
    {
        column[headers[i]] = i;
    }

    auto const get_value = [&](std::vector<std::string> const& fields,
                               std::string const& key) -> double
    { return std::stod(fields.at(column.at(key))); };

    auto const get_optional_value = [&](std::vector<std::string> const& fields,
                                        std::string const& key,
                                        double const default_value) -> double
    {
        auto const it = column.find(key);
        if (it == column.end())
        {
            return default_value;
        }
        return std::stod(fields.at(it->second));
    };

    std::vector<DSMMicroMacroSupportBaselineRow> rows;
    std::string line;
    while (std::getline(in, line))
    {
        if (line.empty())
        {
            continue;
        }

        auto const fields = dsm_micromacroMCCSplitCommaLine(line);
        DSMMicroMacroSupportBaselineRow row;
        row.step = static_cast<int>(get_value(fields, "step"));
        row.pressure = get_value(fields, "pressure");
        row.epsilon_v_total =
            get_optional_value(fields, "epsilon_v_total", 0.0);
        row.saturation = get_value(fields, "S_L");
        row.n_l = get_value(fields, "n_l");
        row.phi_m = get_value(fields, "phi_m");
        row.phi_M = get_value(fields, "phi_M");
        row.phi = get_value(fields, "phi");
        row.n_S = get_value(fields, "n_S");
        row.n_L = get_value(fields, "n_L");
        row.rho_lR = get_value(fields, "rho_lR");
        row.rho_LR = get_value(fields, "rho_LR");
        row.omega_l = get_value(fields, "omega_l");
        row.mu_lR = get_value(fields, "mu_lR");
        row.rho_l_hat = get_value(fields, "rho_l_hat");
        row.delta_epsilon_sw = get_value(fields, "delta_epsilon_sw");
        row.epsilon_sw = get_value(fields, "epsilon_sw");
        row.sigma_S_xx = get_value(fields, "sigma_S_xx");
        rows.push_back(row);
    }

    return rows;
}

double dsm_micromacroMCCScaledTolerance(double const value,
                                  double const relative = 1e-8,
                                  double const absolute = 1e-12)
{
    return absolute + relative * std::max(1.0, std::abs(value));
}
}  // namespace

TEST(MaterialLib_RMBridgeMFront_DSMMicroMacroMCC,
     NeutralDSMMicroMacroStateMatchesVerifiedMCCBridge)
{
    auto parameters = createDSMMicroMacroMCCParameters();
    auto reference = createReferenceMCCModel(parameters);
    auto dsm_micromacro_mcc = createDSMMicroMacroMCCModel(parameters);

    auto reference_state = reference->createMaterialStateVariables();
    auto dsm_micromacro_state = dsm_micromacro_mcc->createMaterialStateVariables();
    initializeState(*reference, *reference_state);
    initializeState(*dsm_micromacro_mcc, *dsm_micromacro_state);

    ParameterLib::SpatialPosition x{};

    MPL::VariableArray previous;
    previous.mechanical_strain.emplace<KV>(KV::Zero());
    previous.stress.emplace<KV>(dsm_micromacroMCCIsotropicStress(-5e3));
    previous.liquid_phase_pressure = -5e3;
    previous.temperature = 293.15;

    auto reference_previous = previous;
    auto dsm_micromacro_previous = previous;

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
            dsm_micromacroMCCIsotropicStrainFromVolumetric(steps[i].eps_v));
        reference_current.liquid_phase_pressure = steps[i].liquid_pressure;

        MPL::VariableArray dsm_micromacro_current = dsm_micromacro_previous;
        dsm_micromacro_current.mechanical_strain.emplace<KV>(
            dsm_micromacroMCCIsotropicStrainFromVolumetric(steps[i].eps_v));
        dsm_micromacro_current.liquid_phase_pressure = steps[i].liquid_pressure;

        double const dt = i == 0 ? steps[i].t : steps[i].t - steps[i - 1].t;

        auto reference_response = reference->integrateStressPressureCoupled(
            reference_previous, reference_current, steps[i].t, x, dt,
            *reference_state);
        auto dsm_micromacro_response = dsm_micromacro_mcc->integrateStressPressureCoupled(
            dsm_micromacro_previous, dsm_micromacro_current, steps[i].t, x, dt,
            *dsm_micromacro_state);

        ASSERT_TRUE(reference_response);
        ASSERT_TRUE(dsm_micromacro_response);
        ASSERT_TRUE(reference_response->state);
        ASSERT_TRUE(dsm_micromacro_response->state);

        EXPECT_TRUE((reference_response->stress - dsm_micromacro_response->stress)
                        .isZero(5e-11));
        EXPECT_TRUE(
            (reference_response->dStress_dStrain -
             dsm_micromacro_response->dStress_dStrain)
                .isZero(1e-8));
        EXPECT_TRUE((reference_response->dStress_dLiquidPressure -
                     dsm_micromacro_response->dStress_dLiquidPressure)
                        .isZero(1e-12));
        EXPECT_TRUE((reference_response->dSaturation_dStrain -
                     dsm_micromacro_response->dSaturation_dStrain)
                        .isZero(1e-12));
        EXPECT_NEAR(reference_response->saturation, dsm_micromacro_response->saturation,
                    1e-15);
        EXPECT_NEAR(reference_response->dSaturation_dLiquidPressure,
                    dsm_micromacro_response->dSaturation_dLiquidPressure, 1e-15);

        for (auto const& name :
             {"EquivalentPlasticStrain", "PreConsolidationPressure",
              "PlasticVolumetricStrain", "VolumeRatio"})
        {
            auto const reference_value =
                getInternalScalar(*reference, *reference_response->state, name);
            auto const dsm_micromacro_value =
                getInternalScalar(*dsm_micromacro_mcc, *dsm_micromacro_response->state,
                                  name);
            EXPECT_NEAR(reference_value, dsm_micromacro_value,
                        name == std::string("PreConsolidationPressure") ? 1e-8
                                                                        : 1e-14)
                << name;
        }

        auto const reference_elastic_strain = getInternalVector(
            *reference, *reference_response->state, "ElasticStrain");
        auto const dsm_micromacro_elastic_strain = getInternalVector(
            *dsm_micromacro_mcc, *dsm_micromacro_response->state, "ElasticStrain");
        ASSERT_EQ(reference_elastic_strain.size(), dsm_micromacro_elastic_strain.size());
        for (std::size_t k = 0; k < reference_elastic_strain.size(); ++k)
        {
            EXPECT_NEAR(reference_elastic_strain[k], dsm_micromacro_elastic_strain[k],
                        1e-14);
        }

        auto const n_l =
            getInternalScalar(*dsm_micromacro_mcc, *dsm_micromacro_response->state, "n_l");
        auto const rho_lR =
            getInternalScalar(*dsm_micromacro_mcc, *dsm_micromacro_response->state,
                              "rho_lR");
        auto const epsilon_sw = getInternalScalar(*dsm_micromacro_mcc,
                                                  *dsm_micromacro_response->state,
                                                  "epsilon_sw");
        auto const phi_m =
            getInternalScalar(*dsm_micromacro_mcc, *dsm_micromacro_response->state, "phi_m");
        auto const phi_M =
            getInternalScalar(*dsm_micromacro_mcc, *dsm_micromacro_response->state, "phi_M");
        auto const mu_lR = getInternalScalar(*dsm_micromacro_mcc,
                                             *dsm_micromacro_response->state,
                                             "mu_lR");
        auto const rho_l_hat =
            getInternalScalar(*dsm_micromacro_mcc, *dsm_micromacro_response->state,
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
        dsm_micromacro_state = std::move(dsm_micromacro_response->state);

        reference_previous = reference_current;
        reference_previous.stress.emplace<KV>(reference_response->stress);
        dsm_micromacro_previous = dsm_micromacro_current;
        dsm_micromacro_previous.stress.emplace<KV>(dsm_micromacro_response->stress);
    }
}

TEST(MaterialLib_RMBridgeMFront_DSMMicroMacroMCC,
     SwellingFeedbackChangesStressButKeepsCarrierSaturation)
{
    auto parameters = createDSMMicroMacroMCCParameters(2.0, 1e-4, 0.1, 1300.0, 0.0);
    auto reference = createReferenceMCCModel(parameters);
    auto dsm_micromacro_mcc = createDSMMicroMacroMCCModel(parameters);

    auto reference_state = reference->createMaterialStateVariables();
    auto dsm_micromacro_state = dsm_micromacro_mcc->createMaterialStateVariables();
    initializeState(*reference, *reference_state);
    initializeState(*dsm_micromacro_mcc, *dsm_micromacro_state);

    ParameterLib::SpatialPosition x{};

    MPL::VariableArray previous;
    previous.mechanical_strain.emplace<KV>(KV::Zero());
    previous.stress.emplace<KV>(dsm_micromacroMCCIsotropicStress(-5e3));
    previous.liquid_phase_pressure = -5e3;
    previous.temperature = 293.15;

    constexpr double t = 1.0;
    constexpr double dt = 1.0;
    constexpr double eps_v = -2e-5;
    constexpr double liquid_pressure = -4e4;

    MPL::VariableArray current = previous;
    current.mechanical_strain.emplace<KV>(
        dsm_micromacroMCCIsotropicStrainFromVolumetric(eps_v));
    current.liquid_phase_pressure = liquid_pressure;

    auto reference_response = reference->integrateStressPressureCoupled(
        previous, current, t, x, dt, *reference_state);
    auto dsm_micromacro_response = dsm_micromacro_mcc->integrateStressPressureCoupled(
        previous, current, t, x, dt, *dsm_micromacro_state);

    ASSERT_TRUE(reference_response);
    ASSERT_TRUE(dsm_micromacro_response);
    ASSERT_TRUE(reference_response->state);
    ASSERT_TRUE(dsm_micromacro_response->state);

    auto const epsilon_sw =
        getInternalScalar(*dsm_micromacro_mcc, *dsm_micromacro_response->state,
                          "epsilon_sw");
    EXPECT_GT(std::abs(epsilon_sw), 1e-12);

    double const bulk_modulus = 52e6 / (3.0 * (1.0 - 2.0 * 0.3));
    auto const expected_stress =
        reference_response->stress -
        dsm_micromacroMCCIsotropicStress(bulk_modulus * epsilon_sw);
    EXPECT_TRUE((expected_stress - dsm_micromacro_response->stress).isZero(5e-8));

    auto const sigma_S_values =
        getInternalVector(*dsm_micromacro_mcc, *dsm_micromacro_response->state, "sigma_S");
    ASSERT_EQ(sigma_S_values.size(), static_cast<std::size_t>(KV::RowsAtCompileTime));
    KV sigma_S = KV::Zero();
    for (Eigen::Index i = 0; i < KV::RowsAtCompileTime; ++i)
    {
        sigma_S[i] = sigma_S_values[static_cast<std::size_t>(i)];
    }
    auto const added_swelling_stress =
        dsm_micromacro_response->stress - reference_response->stress;
    EXPECT_TRUE((added_swelling_stress - sigma_S).isZero(5e-8));

    EXPECT_NEAR(reference_response->saturation, dsm_micromacro_response->saturation,
                1e-15);
    EXPECT_TRUE((reference_response->dSaturation_dStrain -
                 dsm_micromacro_response->dSaturation_dStrain)
                    .isZero(1e-12));
    EXPECT_NEAR(reference_response->dSaturation_dLiquidPressure,
                dsm_micromacro_response->dSaturation_dLiquidPressure, 1e-15);

    auto evaluate_dsm_micromacro_response =
        [&](double const eps_v_value, double const pressure_value)
    {
        auto local_parameters =
            createDSMMicroMacroMCCParameters(2.0, 1e-4, 0.1, 1300.0, 0.0);
        auto local_model = createDSMMicroMacroMCCModel(local_parameters);
        auto local_state = local_model->createMaterialStateVariables();
        initializeState(*local_model, *local_state);

        MPL::VariableArray local_current = previous;
        local_current.mechanical_strain.emplace<KV>(
            dsm_micromacroMCCIsotropicStrainFromVolumetric(eps_v_value));
        local_current.liquid_phase_pressure = pressure_value;

        return local_model->integrateStressPressureCoupled(
            previous, local_current, t, x, dt, *local_state);
    };

    double const dp_fd = std::max(1e-6, std::abs(liquid_pressure) * 1e-6);
    auto dsm_micromacro_p_plus = evaluate_dsm_micromacro_response(eps_v, liquid_pressure + dp_fd);
    auto dsm_micromacro_p_minus = evaluate_dsm_micromacro_response(eps_v, liquid_pressure - dp_fd);
    ASSERT_TRUE(dsm_micromacro_p_plus);
    ASSERT_TRUE(dsm_micromacro_p_minus);
    auto const dsigma_dp_fd =
        (dsm_micromacro_p_plus->stress - dsm_micromacro_p_minus->stress) / (2.0 * dp_fd);
    EXPECT_TRUE((dsm_micromacro_response->dStress_dLiquidPressure - dsigma_dp_fd)
                    .isZero(5e-4));
}

TEST(MaterialLib_RMBridgeMFront_DSMMicroMacroMCC,
     MicroOnlySwellingKeepsMacroExchangeOffAndPreservesCarrierSaturation)
{
    auto parameters = createDSMMicroMacroMCCParameters(2.0, 0.0, 0.1, 1300.0, 0.0);
    auto reference = createReferenceMCCModel(parameters);
    auto dsm_micromacro_mcc = createDSMMicroMacroMCCModel(parameters);

    auto reference_state = reference->createMaterialStateVariables();
    auto dsm_micromacro_state = dsm_micromacro_mcc->createMaterialStateVariables();
    initializeState(*reference, *reference_state);
    initializeState(*dsm_micromacro_mcc, *dsm_micromacro_state);

    ParameterLib::SpatialPosition x{};

    MPL::VariableArray previous;
    previous.mechanical_strain.emplace<KV>(KV::Zero());
    previous.stress.emplace<KV>(dsm_micromacroMCCIsotropicStress(-5e3));
    previous.liquid_phase_pressure = -5e3;
    previous.temperature = 293.15;

    constexpr double t = 1.0;
    constexpr double dt = 1.0;
    constexpr double eps_v = -2e-5;
    constexpr double liquid_pressure = -4e4;

    MPL::VariableArray current = previous;
    current.mechanical_strain.emplace<KV>(
        dsm_micromacroMCCIsotropicStrainFromVolumetric(eps_v));
    current.liquid_phase_pressure = liquid_pressure;

    auto reference_response = reference->integrateStressPressureCoupled(
        previous, current, t, x, dt, *reference_state);
    auto dsm_micromacro_response = dsm_micromacro_mcc->integrateStressPressureCoupled(
        previous, current, t, x, dt, *dsm_micromacro_state);

    ASSERT_TRUE(reference_response);
    ASSERT_TRUE(dsm_micromacro_response);
    ASSERT_TRUE(dsm_micromacro_response->state);

    auto const epsilon_sw =
        getInternalScalar(*dsm_micromacro_mcc, *dsm_micromacro_response->state,
                          "epsilon_sw");
    EXPECT_GT(std::abs(epsilon_sw), 1e-12);

    EXPECT_NEAR(reference_response->liquid_mass_exchange_source, 0.0, 1e-18);
    EXPECT_NEAR(dsm_micromacro_response->liquid_mass_exchange_source, 0.0, 1e-18);

    EXPECT_NEAR(reference_response->saturation, dsm_micromacro_response->saturation,
                1e-15);
    EXPECT_TRUE((reference_response->dSaturation_dStrain -
                 dsm_micromacro_response->dSaturation_dStrain)
                    .isZero(1e-12));
    EXPECT_NEAR(reference_response->dSaturation_dLiquidPressure,
                dsm_micromacro_response->dSaturation_dLiquidPressure, 1e-15);

    auto const sigma_S_values =
        getInternalVector(*dsm_micromacro_mcc, *dsm_micromacro_response->state, "sigma_S");
    ASSERT_EQ(sigma_S_values.size(), static_cast<std::size_t>(KV::RowsAtCompileTime));
    KV sigma_S = KV::Zero();
    for (Eigen::Index i = 0; i < KV::RowsAtCompileTime; ++i)
    {
        sigma_S[i] = sigma_S_values[static_cast<std::size_t>(i)];
    }

    auto const added_swelling_stress =
        dsm_micromacro_response->stress - reference_response->stress;
    EXPECT_TRUE((added_swelling_stress - sigma_S).isZero(5e-8));
}

TEST(MaterialLib_RMBridgeMFront_DSMMicroMacroMCC,
     NativeAlignedStageOneStressGapIsMicroSupportStress)
{
    auto parameters = createDSMMicroMacroMCCParameters(
        0.1, 1e-13, 0.1, 2095.3222465784393, 0.0, 0.0, 0.0, 1.0, 6e-20, 1e-3);
    auto reference = createReferenceMCCModel(parameters);
    auto dsm_micromacro_mcc = createDSMMicroMacroMCCModel(parameters);

    auto reference_state = reference->createMaterialStateVariables();
    auto dsm_micromacro_state = dsm_micromacro_mcc->createMaterialStateVariables();
    initializeState(*reference, *reference_state);
    initializeState(*dsm_micromacro_mcc, *dsm_micromacro_state);

    ParameterLib::SpatialPosition x{};

    MPL::VariableArray previous;
    previous.mechanical_strain.emplace<KV>(KV::Zero());
    previous.stress.emplace<KV>(dsm_micromacroMCCIsotropicStress(0.0));
    previous.liquid_phase_pressure = -1e6;
    previous.temperature = 293.15;

    MPL::VariableArray current = previous;
    current.liquid_phase_pressure = 2e3;

    constexpr double t = 1e3;
    constexpr double dt = 1e3;

    auto reference_response = reference->integrateStressPressureCoupled(
        previous, current, t, x, dt, *reference_state);
    auto dsm_micromacro_response = dsm_micromacro_mcc->integrateStressPressureCoupled(
        previous, current, t, x, dt, *dsm_micromacro_state);

    ASSERT_TRUE(reference_response);
    ASSERT_TRUE(dsm_micromacro_response);
    ASSERT_TRUE(reference_response->state);
    ASSERT_TRUE(dsm_micromacro_response->state);

    auto const sigma_S_values =
        getInternalVector(*dsm_micromacro_mcc, *dsm_micromacro_response->state, "sigma_S");
    ASSERT_EQ(sigma_S_values.size(), static_cast<std::size_t>(KV::RowsAtCompileTime));
    KV sigma_S = KV::Zero();
    for (Eigen::Index i = 0; i < KV::RowsAtCompileTime; ++i)
    {
        sigma_S[i] = sigma_S_values[static_cast<std::size_t>(i)];
    }

    auto const added_swelling_stress =
        dsm_micromacro_response->stress - reference_response->stress;
    EXPECT_GT(sigma_S.norm(), 1e-3);
    EXPECT_TRUE((added_swelling_stress - sigma_S).isZero(1e-6));
    EXPECT_NEAR(reference_response->saturation, dsm_micromacro_response->saturation,
                1e-15);
}

TEST(MaterialLib_RMBridgeMFront_DSMMicroMacroMCC,
     PressureCoupledResponseExposesMacroExchangeSource)
{
    auto parameters = createDSMMicroMacroMCCParameters(
        0.1, 1e-13, 0.01, 2276.031917690513, 0.0, 0.0, 1.0, 1.0, 5.1e-21, 1e-3);
    auto reference = createReferenceMCCModel(parameters);
    auto dsm_micromacro_mcc = createDSMMicroMacroMCCModel(parameters);

    auto reference_state = reference->createMaterialStateVariables();
    auto dsm_micromacro_state = dsm_micromacro_mcc->createMaterialStateVariables();
    initializeState(*reference, *reference_state);
    initializeState(*dsm_micromacro_mcc, *dsm_micromacro_state);

    ParameterLib::SpatialPosition x{};

    MPL::VariableArray previous;
    previous.mechanical_strain.emplace<KV>(KV::Zero());
    previous.stress.emplace<KV>(KV::Zero());
    previous.liquid_phase_pressure = 2e3;
    previous.temperature = 293.15;

    MPL::VariableArray current = previous;

    constexpr double t = 1e5;
    constexpr double dt = 1e5;

    auto reference_response = reference->integrateStressPressureCoupled(
        previous, current, t, x, dt, *reference_state);
    auto dsm_micromacro_response = dsm_micromacro_mcc->integrateStressPressureCoupled(
        previous, current, t, x, dt, *dsm_micromacro_state);

    ASSERT_TRUE(reference_response);
    ASSERT_TRUE(dsm_micromacro_response);
    ASSERT_TRUE(dsm_micromacro_response->state);

    auto const rho_l_hat =
        getInternalScalar(*dsm_micromacro_mcc, *dsm_micromacro_response->state,
                          "rho_l_hat");
    auto const sigma_S_values =
        getInternalVector(*dsm_micromacro_mcc, *dsm_micromacro_response->state, "sigma_S");
    KV sigma_S = KV::Zero();
    for (Eigen::Index i = 0; i < KV::RowsAtCompileTime; ++i)
    {
        sigma_S[i] = sigma_S_values[static_cast<std::size_t>(i)];
    }

    EXPECT_NEAR(reference_response->liquid_mass_exchange_source, 0.0, 1e-18);
    EXPECT_GT(rho_l_hat, 0.0);
    EXPECT_NEAR(dsm_micromacro_response->liquid_mass_exchange_source, -rho_l_hat,
                std::max(1e-18, 1e-12 * std::abs(rho_l_hat)));
    EXPECT_TRUE(
        (dsm_micromacro_response->swelling_stress - sigma_S).isZero(1e-12));
}

TEST(MaterialLib_RMBridgeMFront_DSMMicroMacroMCC,
     NegativeAttractiveConventionMatchesLegacyNegativeHamakerPath)
{
    auto legacy_parameters = createDSMMicroMacroMCCParameters(
        8.0, 1e-13, 0.1, 1300.0, 0.0, 0.0, 0.0, 0.0, -6e-20);
    auto native_aligned_parameters = createDSMMicroMacroMCCParameters(
        8.0, 1e-13, 0.1, 1300.0, 0.0, 0.0, 0.0, 1.0, 6e-20);

    auto legacy_model = createDSMMicroMacroMCCModel(legacy_parameters);
    auto native_aligned_model =
        createDSMMicroMacroMCCModel(native_aligned_parameters);

    auto legacy_state = legacy_model->createMaterialStateVariables();
    auto native_aligned_state =
        native_aligned_model->createMaterialStateVariables();
    initializeState(*legacy_model, *legacy_state);
    initializeState(*native_aligned_model, *native_aligned_state);

    ParameterLib::SpatialPosition x{};

    MPL::VariableArray previous;
    previous.mechanical_strain.emplace<KV>(KV::Zero());
    previous.stress.emplace<KV>(dsm_micromacroMCCIsotropicStress(0.0));
    previous.liquid_phase_pressure = -1e6;
    previous.temperature = 293.15;

    MPL::VariableArray current = previous;
    current.liquid_phase_pressure = 2e3;

    constexpr double t = 1e5;
    constexpr double dt = 1e5;

    auto legacy_response = legacy_model->integrateStressPressureCoupled(
        previous, current, t, x, dt, *legacy_state);
    auto native_aligned_response =
        native_aligned_model->integrateStressPressureCoupled(
            previous, current, t, x, dt, *native_aligned_state);

    ASSERT_TRUE(legacy_response);
    ASSERT_TRUE(native_aligned_response);
    ASSERT_TRUE(legacy_response->state);
    ASSERT_TRUE(native_aligned_response->state);

    EXPECT_TRUE((legacy_response->stress - native_aligned_response->stress)
                    .isZero(1e-8));
    EXPECT_TRUE(
        (legacy_response->dStress_dStrain -
         native_aligned_response->dStress_dStrain)
            .isZero(1e-8));
    EXPECT_TRUE((legacy_response->dStress_dLiquidPressure -
                 native_aligned_response->dStress_dLiquidPressure)
                    .isZero(1e-8));
    EXPECT_NEAR(legacy_response->saturation, native_aligned_response->saturation,
                1e-15);
    EXPECT_NEAR(legacy_response->dSaturation_dLiquidPressure,
                native_aligned_response->dSaturation_dLiquidPressure, 1e-15);

    auto const legacy_epsilon_sw =
        getInternalScalar(*legacy_model, *legacy_response->state, "epsilon_sw");
    auto const aligned_epsilon_sw = getInternalScalar(
        *native_aligned_model, *native_aligned_response->state, "epsilon_sw");
    EXPECT_GT(std::abs(legacy_epsilon_sw), 1e-8);
    EXPECT_NEAR(legacy_epsilon_sw, aligned_epsilon_sw, 1e-12);
}

TEST(MaterialLib_RMBridgeMFront_DSMMicroMacroMCC,
     DSMMicroMacroSaturationModeMatchesTullerLawAndKeepsStressSurface)
{
    auto parameters = createDSMMicroMacroMCCParameters(0.0, 0.0, 0.1, 1300.0, 0.0,
                                                  1.0);
    auto reference = createReferenceMCCModel(parameters);
    auto dsm_micromacro_mcc = createDSMMicroMacroMCCModel(parameters);

    auto reference_state = reference->createMaterialStateVariables();
    auto dsm_micromacro_state = dsm_micromacro_mcc->createMaterialStateVariables();
    initializeState(*reference, *reference_state);
    initializeState(*dsm_micromacro_mcc, *dsm_micromacro_state);

    ParameterLib::SpatialPosition x{};

    MPL::VariableArray previous;
    previous.mechanical_strain.emplace<KV>(KV::Zero());
    previous.stress.emplace<KV>(dsm_micromacroMCCIsotropicStress(-5e3));
    previous.liquid_phase_pressure = -5e3;
    previous.temperature = 293.15;

    constexpr double t = 1.0;
    constexpr double dt = 1.0;
    constexpr double eps_v = -2e-5;
    constexpr double liquid_pressure = -4e4;

    MPL::VariableArray current = previous;
    current.mechanical_strain.emplace<KV>(
        dsm_micromacroMCCIsotropicStrainFromVolumetric(eps_v));
    current.liquid_phase_pressure = liquid_pressure;

    auto reference_response = reference->integrateStressPressureCoupled(
        previous, current, t, x, dt, *reference_state);
    auto dsm_micromacro_response = dsm_micromacro_mcc->integrateStressPressureCoupled(
        previous, current, t, x, dt, *dsm_micromacro_state);

    ASSERT_TRUE(reference_response);
    ASSERT_TRUE(dsm_micromacro_response);
    ASSERT_TRUE(reference_response->state);
    ASSERT_TRUE(dsm_micromacro_response->state);

    EXPECT_TRUE((reference_response->stress - dsm_micromacro_response->stress)
                    .isZero(5e-11));
    EXPECT_TRUE((reference_response->dStress_dStrain -
                 dsm_micromacro_response->dStress_dStrain)
                    .isZero(1e-8));
    EXPECT_TRUE((reference_response->dStress_dLiquidPressure -
                 dsm_micromacro_response->dStress_dLiquidPressure)
                    .isZero(1e-12));

    EXPECT_NEAR(dsm_micromacro_response->saturation,
                dsm_micromacroTullerSaturation(liquid_pressure), 1e-15);
    EXPECT_NEAR(dsm_micromacro_response->dSaturation_dLiquidPressure,
                dsm_micromacroTullerDSaturationDLiquidPressure(liquid_pressure),
                1e-15);
    EXPECT_TRUE(dsm_micromacro_response->dSaturation_dStrain.isZero(1e-12));

    auto const epsilon_sw =
        getInternalScalar(*dsm_micromacro_mcc, *dsm_micromacro_response->state,
                          "epsilon_sw");
    EXPECT_NEAR(epsilon_sw, 0.0, 1e-14);

    MPL::VariableArray saturated_current = previous;
    saturated_current.mechanical_strain.emplace<KV>(
        dsm_micromacroMCCIsotropicStrainFromVolumetric(eps_v));
    saturated_current.liquid_phase_pressure = 1e3;

    auto saturated_response = dsm_micromacro_mcc->integrateStressPressureCoupled(
        previous, saturated_current, t, x, dt, *dsm_micromacro_state);
    ASSERT_TRUE(saturated_response);
    EXPECT_NEAR(saturated_response->saturation, 1.0, 1e-15);
    EXPECT_NEAR(saturated_response->dSaturation_dLiquidPressure, 0.0, 1e-15);
    EXPECT_TRUE(saturated_response->dSaturation_dStrain.isZero(1e-12));
}

TEST(MaterialLib_RMBridgeMFront_DSMMicroMacroMCC,
     DSMMicroMacroSupportStateMatchesOverlapTransferBaseline)
{
    auto const baseline_rows = loadDSMMicroMacroSupportBaselineRows(
        TestInfoLib::TestInfo::data_path +
        "/MaterialLib/MFront/RichardsMechanicsDSMMicroMacroBridge_overlap_transfer_baseline.csv");
    ASSERT_EQ(baseline_rows.size(), 5);

    auto parameters = createDSMMicroMacroSupportParameters();
    auto dsm_micromacro_mcc = createDSMMicroMacroMCCModel(parameters);

    auto dsm_micromacro_state = dsm_micromacro_mcc->createMaterialStateVariables();
    initializeState(*dsm_micromacro_mcc, *dsm_micromacro_state);

    ParameterLib::SpatialPosition x{};

    MPL::VariableArray previous;
    previous.mechanical_strain.emplace<KV>(KV::Zero());
    previous.stress.emplace<KV>(KV::Zero());
    previous.liquid_phase_pressure = 0.0;
    previous.liquid_saturation = 1.0;
    previous.temperature = 293.15;

    for (auto const& row : baseline_rows)
    {
        MPL::VariableArray current = previous;
        current.liquid_phase_pressure = row.pressure;

        auto dsm_micromacro_response = dsm_micromacro_mcc->integrateStressPressureCoupled(
            previous, current, static_cast<double>(row.step + 1), x, 1.0,
            *dsm_micromacro_state);

        ASSERT_TRUE(dsm_micromacro_response);
        ASSERT_TRUE(dsm_micromacro_response->state);

        auto const n_l =
            getInternalScalar(*dsm_micromacro_mcc, *dsm_micromacro_response->state, "n_l");
        auto const phi_m =
            getInternalScalar(*dsm_micromacro_mcc, *dsm_micromacro_response->state, "phi_m");
        auto const phi_M =
            getInternalScalar(*dsm_micromacro_mcc, *dsm_micromacro_response->state, "phi_M");
        auto const phi =
            getInternalScalar(*dsm_micromacro_mcc, *dsm_micromacro_response->state, "phi");
        auto const n_S =
            getInternalScalar(*dsm_micromacro_mcc, *dsm_micromacro_response->state, "n_S");
        auto const n_L =
            getInternalScalar(*dsm_micromacro_mcc, *dsm_micromacro_response->state, "n_L");
        auto const rho_lR =
            getInternalScalar(*dsm_micromacro_mcc, *dsm_micromacro_response->state,
                              "rho_lR");
        auto const rho_LR =
            getInternalScalar(*dsm_micromacro_mcc, *dsm_micromacro_response->state,
                              "rho_LR");
        auto const omega_l =
            getInternalScalar(*dsm_micromacro_mcc, *dsm_micromacro_response->state,
                              "omega_l");
        auto const mu_lR =
            getInternalScalar(*dsm_micromacro_mcc, *dsm_micromacro_response->state,
                              "mu_lR");
        auto const rho_l_hat =
            getInternalScalar(*dsm_micromacro_mcc, *dsm_micromacro_response->state,
                              "rho_l_hat");
        auto const delta_epsilon_sw = getInternalScalar(
            *dsm_micromacro_mcc, *dsm_micromacro_response->state, "delta_epsilon_sw");
        auto const epsilon_sw = getInternalScalar(
            *dsm_micromacro_mcc, *dsm_micromacro_response->state, "epsilon_sw");
        auto const sigma_S =
            getInternalVector(*dsm_micromacro_mcc, *dsm_micromacro_response->state,
                              "sigma_S");

        EXPECT_NEAR(dsm_micromacro_response->saturation, row.saturation,
                    dsm_micromacroMCCScaledTolerance(row.saturation));
        EXPECT_NEAR(n_l, row.n_l,
                    dsm_micromacroMCCScaledTolerance(row.n_l, 5e-3, 1e-8));
        EXPECT_NEAR(phi_m, row.phi_m,
                    dsm_micromacroMCCScaledTolerance(row.phi_m, 5e-3, 1e-8));
        EXPECT_NEAR(phi_M, row.phi_M,
                    dsm_micromacroMCCScaledTolerance(row.phi_M, 5e-3, 1e-8));
        EXPECT_NEAR(phi, row.phi,
                    dsm_micromacroMCCScaledTolerance(row.phi, 1e-12, 1e-12));
        EXPECT_NEAR(n_S, row.n_S,
                    dsm_micromacroMCCScaledTolerance(row.n_S, 1e-12, 1e-12));
        EXPECT_NEAR(n_L, row.n_L,
                    dsm_micromacroMCCScaledTolerance(row.n_L, 5e-3, 1e-8));
        EXPECT_NEAR(rho_lR, row.rho_lR,
                    dsm_micromacroMCCScaledTolerance(row.rho_lR, 1e-3, 1e-6));
        EXPECT_NEAR(rho_LR, row.rho_LR,
                    dsm_micromacroMCCScaledTolerance(row.rho_LR, 1e-12, 1e-12));
        EXPECT_NEAR(omega_l, row.omega_l,
                    dsm_micromacroMCCScaledTolerance(row.omega_l, 5e-3, 1e-8));
        EXPECT_NEAR(mu_lR, row.mu_lR,
                    dsm_micromacroMCCScaledTolerance(row.mu_lR, 1e-2, 1e-8));
        EXPECT_NEAR(rho_l_hat, row.rho_l_hat,
                    dsm_micromacroMCCScaledTolerance(row.rho_l_hat, 1e-2, 1e-8));
        EXPECT_NEAR(delta_epsilon_sw, row.delta_epsilon_sw,
                    dsm_micromacroMCCScaledTolerance(row.delta_epsilon_sw, 5e-3,
                                               1e-8));
        EXPECT_NEAR(epsilon_sw, row.epsilon_sw,
                    dsm_micromacroMCCScaledTolerance(row.epsilon_sw, 5e-3, 1e-8));
        ASSERT_EQ(sigma_S.size(), 6u);
        EXPECT_NEAR(sigma_S[0], row.sigma_S_xx,
                    dsm_micromacroMCCScaledTolerance(row.sigma_S_xx, 5e-3, 1e-6));
        EXPECT_NEAR(sigma_S[1], row.sigma_S_xx,
                    dsm_micromacroMCCScaledTolerance(row.sigma_S_xx, 5e-3, 1e-6));
        EXPECT_NEAR(sigma_S[2], row.sigma_S_xx,
                    dsm_micromacroMCCScaledTolerance(row.sigma_S_xx, 5e-3, 1e-6));
        EXPECT_NEAR(sigma_S[3], 0.0, 1e-12);
        EXPECT_NEAR(sigma_S[4], 0.0, 1e-12);
        EXPECT_NEAR(sigma_S[5], 0.0, 1e-12);

        dsm_micromacro_state = std::move(dsm_micromacro_response->state);
        dsm_micromacro_state->pushBackState();

        previous = current;
        previous.stress.emplace<KV>(dsm_micromacro_response->stress);
        previous.liquid_saturation = dsm_micromacro_response->saturation;
    }
}

TEST(MaterialLib_RMBridgeMFront_DSMMicroMacroMCC,
     DSMMicroMacroSupportStateMatchesStrainCoupledBaseline)
{
    auto const overlap_rows = loadDSMMicroMacroSupportBaselineRows(
        TestInfoLib::TestInfo::data_path +
        "/MaterialLib/MFront/RichardsMechanicsDSMMicroMacroBridge_overlap_transfer_baseline.csv");
    auto const strain_rows = loadDSMMicroMacroSupportBaselineRows(
        TestInfoLib::TestInfo::data_path +
        "/MaterialLib/MFront/RichardsMechanicsDSMMicroMacroBridge_strain_coupled_overlap_baseline.csv");
    ASSERT_EQ(overlap_rows.size(), 5);
    ASSERT_EQ(strain_rows.size(), 5);
    auto const& anchor = overlap_rows.back();

    auto parameters = createDSMMicroMacroSupportParameters();
    auto dsm_micromacro_mcc = createDSMMicroMacroMCCModel(parameters);

    auto dsm_micromacro_state = dsm_micromacro_mcc->createMaterialStateVariables();
    initializeState(*dsm_micromacro_mcc, *dsm_micromacro_state);
    setInternalScalar(*dsm_micromacro_mcc, *dsm_micromacro_state, "n_l", anchor.n_l);
    setInternalScalar(*dsm_micromacro_mcc, *dsm_micromacro_state, "rho_lR", anchor.rho_lR);
    setInternalScalar(*dsm_micromacro_mcc, *dsm_micromacro_state, "epsilon_sw",
                      anchor.epsilon_sw);

    KV anchor_stress = KV::Zero();
    anchor_stress[0] = anchor.sigma_S_xx;
    anchor_stress[1] = anchor.sigma_S_xx;
    anchor_stress[2] = anchor.sigma_S_xx;
    setDSMMicroMacroMCCThermodynamicForces(*dsm_micromacro_state, anchor_stress,
                                      anchor.saturation);
    dsm_micromacro_state->pushBackState();

    ParameterLib::SpatialPosition x{};

    MPL::VariableArray previous;
    previous.mechanical_strain.emplace<KV>(KV::Zero());
    previous.stress.emplace<KV>(anchor_stress);
    previous.liquid_phase_pressure = anchor.pressure;
    previous.liquid_saturation = anchor.saturation;
    previous.temperature = 293.15;

    for (auto const& row : strain_rows)
    {
        MPL::VariableArray current = previous;
        current.liquid_phase_pressure = row.pressure;
        current.mechanical_strain.emplace<KV>(
            dsm_micromacroMCCIsotropicStrainFromVolumetric(row.epsilon_v_total));

        auto dsm_micromacro_response = dsm_micromacro_mcc->integrateStressPressureCoupled(
            previous, current, static_cast<double>(row.step + 1), x, 1.0,
            *dsm_micromacro_state);

        ASSERT_TRUE(dsm_micromacro_response);
        ASSERT_TRUE(dsm_micromacro_response->state);

        auto const n_l =
            getInternalScalar(*dsm_micromacro_mcc, *dsm_micromacro_response->state, "n_l");
        auto const phi_m =
            getInternalScalar(*dsm_micromacro_mcc, *dsm_micromacro_response->state, "phi_m");
        auto const phi_M =
            getInternalScalar(*dsm_micromacro_mcc, *dsm_micromacro_response->state, "phi_M");
        auto const phi =
            getInternalScalar(*dsm_micromacro_mcc, *dsm_micromacro_response->state, "phi");
        auto const n_S =
            getInternalScalar(*dsm_micromacro_mcc, *dsm_micromacro_response->state, "n_S");
        auto const n_L =
            getInternalScalar(*dsm_micromacro_mcc, *dsm_micromacro_response->state, "n_L");
        auto const rho_lR =
            getInternalScalar(*dsm_micromacro_mcc, *dsm_micromacro_response->state,
                              "rho_lR");
        auto const rho_LR =
            getInternalScalar(*dsm_micromacro_mcc, *dsm_micromacro_response->state,
                              "rho_LR");
        auto const omega_l =
            getInternalScalar(*dsm_micromacro_mcc, *dsm_micromacro_response->state,
                              "omega_l");
        auto const mu_lR =
            getInternalScalar(*dsm_micromacro_mcc, *dsm_micromacro_response->state,
                              "mu_lR");
        auto const rho_l_hat =
            getInternalScalar(*dsm_micromacro_mcc, *dsm_micromacro_response->state,
                              "rho_l_hat");
        auto const delta_epsilon_sw = getInternalScalar(
            *dsm_micromacro_mcc, *dsm_micromacro_response->state, "delta_epsilon_sw");
        auto const epsilon_sw = getInternalScalar(
            *dsm_micromacro_mcc, *dsm_micromacro_response->state, "epsilon_sw");
        auto const sigma_S =
            getInternalVector(*dsm_micromacro_mcc, *dsm_micromacro_response->state,
                              "sigma_S");

        EXPECT_NEAR(dsm_micromacro_response->saturation, row.saturation,
                    dsm_micromacroMCCScaledTolerance(row.saturation));
        EXPECT_NEAR(n_l, row.n_l,
                    dsm_micromacroMCCScaledTolerance(row.n_l, 5e-3, 1e-8));
        EXPECT_NEAR(phi_m, row.phi_m,
                    dsm_micromacroMCCScaledTolerance(row.phi_m, 5e-3, 1e-8));
        EXPECT_NEAR(phi_M, row.phi_M,
                    dsm_micromacroMCCScaledTolerance(row.phi_M, 5e-3, 1e-8));
        EXPECT_NEAR(phi, row.phi,
                    dsm_micromacroMCCScaledTolerance(row.phi, 1e-12, 1e-12));
        EXPECT_NEAR(n_S, row.n_S,
                    dsm_micromacroMCCScaledTolerance(row.n_S, 1e-12, 1e-12));
        EXPECT_NEAR(n_L, row.n_L,
                    dsm_micromacroMCCScaledTolerance(row.n_L, 5e-3, 1e-8));
        EXPECT_NEAR(rho_lR, row.rho_lR,
                    dsm_micromacroMCCScaledTolerance(row.rho_lR, 1e-3, 1e-6));
        EXPECT_NEAR(rho_LR, row.rho_LR,
                    dsm_micromacroMCCScaledTolerance(row.rho_LR, 1e-12, 1e-12));
        EXPECT_NEAR(omega_l, row.omega_l,
                    dsm_micromacroMCCScaledTolerance(row.omega_l, 5e-3, 1e-8));
        EXPECT_NEAR(mu_lR, row.mu_lR,
                    dsm_micromacroMCCScaledTolerance(row.mu_lR, 1e-2, 1e-8));
        EXPECT_NEAR(rho_l_hat, row.rho_l_hat,
                    dsm_micromacroMCCScaledTolerance(row.rho_l_hat, 1e-2, 1e-8));
        EXPECT_NEAR(delta_epsilon_sw, row.delta_epsilon_sw,
                    dsm_micromacroMCCScaledTolerance(row.delta_epsilon_sw, 5e-3,
                                               1e-8));
        EXPECT_NEAR(epsilon_sw, row.epsilon_sw,
                    dsm_micromacroMCCScaledTolerance(row.epsilon_sw, 5e-3, 1e-8));
        double const bulk_modulus = 1e10 / (3.0 * (1.0 - 2.0 * 0.25));
        double const expected_sigma_S_xx = -bulk_modulus * row.epsilon_sw;
        ASSERT_EQ(sigma_S.size(), 6u);
        EXPECT_NEAR(sigma_S[0], expected_sigma_S_xx,
                    dsm_micromacroMCCScaledTolerance(expected_sigma_S_xx, 5e-3,
                                               1e-6));
        EXPECT_NEAR(sigma_S[1], expected_sigma_S_xx,
                    dsm_micromacroMCCScaledTolerance(expected_sigma_S_xx, 5e-3,
                                               1e-6));
        EXPECT_NEAR(sigma_S[2], expected_sigma_S_xx,
                    dsm_micromacroMCCScaledTolerance(expected_sigma_S_xx, 5e-3,
                                               1e-6));
        EXPECT_NEAR(sigma_S[3], 0.0, 1e-12);
        EXPECT_NEAR(sigma_S[4], 0.0, 1e-12);
        EXPECT_NEAR(sigma_S[5], 0.0, 1e-12);

        dsm_micromacro_state = std::move(dsm_micromacro_response->state);
        dsm_micromacro_state->pushBackState();

        previous = current;
        previous.stress.emplace<KV>(dsm_micromacro_response->stress);
        previous.liquid_saturation = dsm_micromacro_response->saturation;
    }
}

#endif  // OGS_USE_MFRONT
