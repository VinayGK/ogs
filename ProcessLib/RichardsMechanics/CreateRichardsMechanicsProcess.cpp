// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "CreateRichardsMechanicsProcess.h"

#include <algorithm>
#include <cassert>

#include "BaseLib/DemangleTypeInfo.h"
#include "InfoLib/GitInfo.h"
#include "MaterialLib/MPL/CreateMaterialSpatialDistributionMap.h"
#include "MaterialLib/MPL/MaterialSpatialDistributionMap.h"
#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/SolidModels/CreateConstitutiveRelation.h"
#include "MaterialLib/SolidModels/MechanicsBase.h"
#include "NumLib/CreateNewtonRaphsonSolverParameters.h"
#include "ParameterLib/Utils.h"
#include "ProcessLib/Common/HydroMechanics/CreateInitialStress.h"
#include "ProcessLib/Output/CreateSecondaryVariables.h"
#include "ProcessLib/Utils/ProcessUtils.h"
#include "RichardsMechanicsProcess.h"
#include "RichardsMechanicsProcessData.h"

namespace ProcessLib
{
namespace RichardsMechanics
{
namespace
{
char const* toString(MaterialLib::Solids::ConstitutiveModel const model)
{
    using MaterialLib::Solids::ConstitutiveModel;
    switch (model)
    {
        case ConstitutiveModel::Ehlers:
            return "Ehlers";
        case ConstitutiveModel::LinearElasticIsotropic:
            return "LinearElasticIsotropic";
        case ConstitutiveModel::Lubby2:
            return "Lubby2";
        case ConstitutiveModel::CreepBGRa:
            return "CreepBGRa";
        case ConstitutiveModel::Invalid:
            return "Invalid";
    }
    return "Unknown";
}

char const* toString(VKPotentialExchangeMode const mode)
{
    switch (mode)
    {
        case VKPotentialExchangeMode::FullPotential:
            return "full_potential";
    }
    return "unknown";
}

VKMicroPotentialConvention parseVKMicroPotentialConvention(
    std::string const& convention)
{
    if (convention == "positive_reduced")
    {
        return VKMicroPotentialConvention::PositiveReduced;
    }
    if (convention == "negative_attractive")
    {
        return VKMicroPotentialConvention::NegativeAttractive;
    }

    OGS_FATAL(
        "RichardsMechanics: unsupported vk_potential_exchange "
        "micro_potential_convention '{}'. Currently supported: "
        "'positive_reduced', 'negative_attractive'.",
        convention);
}

VKLocalNonlinearSolveMode parseVKLocalNonlinearSolveMode(
    std::string const& mode)
{
    if (mode == "scalar_exchange")
    {
        return VKLocalNonlinearSolveMode::ScalarExchange;
    }
    if (mode == "scalar_notebook_storage")
    {
        return VKLocalNonlinearSolveMode::ScalarNotebookStorage;
    }
    if (mode == "scalar_notebook_mass_storage")
    {
        return VKLocalNonlinearSolveMode::ScalarNotebookMassStorage;
    }

    OGS_FATAL(
        "RichardsMechanics: unsupported vk_potential_exchange "
        "local_nonlinear_solve_mode '{}'. Currently supported: "
        "'scalar_exchange', 'scalar_notebook_storage', "
        "'scalar_notebook_mass_storage'.",
        mode);
}

VKMacroPorosityUpdateMode parseVKMacroPorosityUpdateMode(
    std::string const& mode)
{
    if (mode == "algebraic_split")
    {
        return VKMacroPorosityUpdateMode::AlgebraicSplit;
    }
    if (mode == "notebook_additive_rate")
    {
        return VKMacroPorosityUpdateMode::NotebookAdditiveRate;
    }

    OGS_FATAL(
        "RichardsMechanics: unsupported vk_potential_exchange "
        "macro_porosity_update_mode '{}'. Currently supported: "
        "'algebraic_split', 'notebook_additive_rate'.",
        mode);
}

VKMicroSolidVolumeFractionMode parseVKMicroSolidVolumeFractionMode(
    std::string const& mode)
{
    if (mode == "reference")
    {
        return VKMicroSolidVolumeFractionMode::Reference;
    }
    if (mode == "current_porosity_split")
    {
        return VKMicroSolidVolumeFractionMode::CurrentPorositySplit;
    }

    OGS_FATAL(
        "RichardsMechanics: unsupported vk_potential_exchange "
        "micro_solid_volume_fraction_mode '{}'. Currently supported: "
        "'reference', 'current_porosity_split'.",
        mode);
}

VKPotentialExchangeRoleMapping parseVKPotentialExchangeRoleMapping(
    std::string const& mapping)
{
    if (mapping == "current_ogs")
    {
        return VKPotentialExchangeRoleMapping::CurrentOgs;
    }
    if (mapping == "notebook_roles")
    {
        return VKPotentialExchangeRoleMapping::NotebookRoles;
    }

    OGS_FATAL(
        "RichardsMechanics: unsupported vk_potential_exchange "
        "potential_role_mapping '{}'. Currently supported: 'current_ogs', "
        "'notebook_roles'.",
        mapping);
}

VKPotentialExchangeMode parseVKPotentialExchangeMode(std::string const& mode)
{
    if (mode == "full_potential")
    {
        return VKPotentialExchangeMode::FullPotential;
    }

    OGS_FATAL(
        "RichardsMechanics: unsupported vk_potential_exchange mode '{}'. "
        "Currently supported: 'full_potential'.",
        mode);
}

template <int DisplacementDim>
void logPhase0TransitionAudit(
    std::string const& process_name,
    bool const use_monolithic_scheme,
    ProcessVariable const& pressure_process_variable,
    ProcessVariable const& displacement_process_variable,
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media,
    std::map<int,
             std::shared_ptr<MaterialLib::Solids::MechanicsBase<DisplacementDim>>> const&
        solid_constitutive_relations,
    std::optional<MicroPorosityParameters> const& micro_porosity_parameters,
    std::optional<VKPotentialExchangeParameters> const&
        vk_potential_exchange_parameters,
    std::map<int, VKPotentialExchangeParameters> const&
        vk_potential_exchange_parameters_by_material)
{
    namespace MPL = MaterialPropertyLib;

    INFO(
        "[RM Phase0 audit] process='{}', OGS='{}', scheme='{}', pressure PV='{}' (1 comp), displacement PV='{}' ({} comp).",
        process_name, GitInfoLib::GitInfo::ogs_version,
        use_monolithic_scheme ? "monolithic" : "staggered",
        pressure_process_variable.getName(), displacement_process_variable.getName(),
        displacement_process_variable.getNumberOfGlobalComponents());

    INFO(
        "[RM Phase0 audit] Hydraulic convention in DS-RM local assembler: capillary pressure is derived from the FE pressure variable as p_c = -p_L (e.g. interpolate(-p_L, ...)); the process config has no explicit gas-pressure parameter.");

    INFO(
        "[RM Phase0 audit] Current DS-RM implementation also sets MPL variable gas_phase_pressure = 1.0e5 Pa in RichardsMechanicsFEM-impl.h (hardcoded, with TODO comment).");

    if (micro_porosity_parameters)
    {
        INFO(
            "[RM Phase0 audit] Micro-porosity constitutive hook: ENABLED (mass_exchange_coefficient = {}).",
            micro_porosity_parameters->mass_exchange_coefficient);
    }
    else
    {
        INFO("[RM Phase0 audit] Micro-porosity constitutive hook: DISABLED.");
    }

    if (vk_potential_exchange_parameters)
    {
        auto const& vkp = *vk_potential_exchange_parameters;
        INFO(
            "[RM Phase0 audit] VK potential-exchange config block: PRESENT (enabled={}, mode='{}', pressure_tolerance={} Pa, hamaker_constant={}, specific_surface={}, rho_SR_ref={}, n_S_ref={}, rho_l0={}, a_rho={}, b_rho={}, micro_potential_convention='{}', local_nonlinear_solve_mode='{}', macro_porosity_update_mode='{}', micro_solid_volume_fraction_mode='{}', potential_role_mapping='{}', initial_n_l={}, fd_jacobian_for_exchange={}, fd_jacobian_perturbation={}, check_local_jacobian={}, local_jacobian_perturbation={}, local_jacobian_relative_tolerance={}, vdw_relaxation_stress_gain={}, micro_water_content_stress_gain={}, micro_water_content_swelling_slope={} ).",
            vkp.enabled ? "true" : "false", toString(vkp.mode),
            vkp.pressure_tolerance, vkp.hamaker_constant, vkp.specific_surface,
            vkp.micro_solid_density_reference,
            vkp.micro_solid_volume_fraction_reference,
            vkp.micro_liquid_density_reference, vkp.micro_liquid_density_a,
            vkp.micro_liquid_density_b,
            toString(vkp.micro_potential_convention),
            toString(vkp.local_nonlinear_solve_mode),
            toString(vkp.macro_porosity_update_mode),
            toString(vkp.micro_solid_volume_fraction_mode),
            toString(vkp.potential_role_mapping),
            vkp.initial_micro_water_content
                ? std::to_string(*vkp.initial_micro_water_content)
                : std::string{"<unset>"},
            vkp.use_fd_jacobian_for_exchange ? "true" : "false",
            vkp.fd_jacobian_perturbation,
            vkp.check_local_jacobian ? "true" : "false",
            vkp.local_jacobian_perturbation,
            vkp.local_jacobian_relative_tolerance,
            vkp.vdw_relaxation_stress_gain,
            vkp.micro_water_content_stress_gain,
            vkp.micro_water_content_swelling_slope);
    }
    else
    {
        INFO("[RM Phase0 audit] VK potential-exchange config block: ABSENT.");
    }

    if (!vk_potential_exchange_parameters_by_material.empty())
    {
        INFO(
            "[RM Phase0 audit] VK potential-exchange medium-specific overrides: {} material ids configured.",
            vk_potential_exchange_parameters_by_material.size());
    }
    else
    {
        INFO(
            "[RM Phase0 audit] VK potential-exchange medium-specific overrides: none.");
    }

    INFO(
        "[RM Phase0 audit] Output/postprocessing note: the RM primary process variable named '{}' is the FE pressure variable (liquid-pressure convention in the current DS-RM code path); capillary pressure is an internal derived quantity.",
        pressure_process_variable.getName());

    for (auto const& [material_id, medium] : media)
    {
        auto const& solid_phase = medium->phase("Solid");
        bool const has_swelling =
            solid_phase.hasProperty(MPL::PropertyType::swelling_stress_rate);
        bool const has_saturation_micro =
            medium->hasProperty(MPL::PropertyType::saturation_micro);

        auto const cr_it = solid_constitutive_relations.find(material_id);
        char const* solid_model_name = "Missing";
        std::string solid_model_type = "Missing";
        if (cr_it != solid_constitutive_relations.end() && cr_it->second)
        {
            solid_model_name = toString(cr_it->second->getConstitutiveModel());
            solid_model_type = BaseLib::demangle(typeid(*cr_it->second).name());
        }

        INFO(
            "[RM Phase0 audit] material_id={} solid_model_enum={} solid_model_type='{}' swelling_stress_rate={} saturation_micro={} bishops_effective_stress={} saturation={} porosity={} relative_permeability={}.",
            material_id, solid_model_name, solid_model_type,
            has_swelling ? "yes" : "no",
            has_saturation_micro ? "yes" : "no",
            medium->hasProperty(MPL::PropertyType::bishops_effective_stress)
                ? "yes"
                : "no",
            medium->hasProperty(MPL::PropertyType::saturation) ? "yes" : "no",
            medium->hasProperty(MPL::PropertyType::porosity) ? "yes" : "no",
            medium->hasProperty(MPL::PropertyType::relative_permeability)
                ? "yes"
                : "no");
    }
}
}  // namespace

void checkMPLProperties(
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media)
{
    std::array const required_medium_properties = {
        MaterialPropertyLib::reference_temperature,
        MaterialPropertyLib::bishops_effective_stress,
        MaterialPropertyLib::relative_permeability,
        MaterialPropertyLib::saturation,
        MaterialPropertyLib::porosity,
        MaterialPropertyLib::biot_coefficient};
    std::array const required_liquid_properties = {
        MaterialPropertyLib::viscosity, MaterialPropertyLib::density};
    std::array const required_solid_properties = {MaterialPropertyLib::density};

    for (auto const& m : media)
    {
        checkRequiredProperties(*m.second, required_medium_properties);
        checkRequiredProperties(m.second->phase("AqueousLiquid"),
                                required_liquid_properties);
        checkRequiredProperties(m.second->phase("Solid"),
                                required_solid_properties);
    }
}

void validateMicroPorosityAndVKConfiguration(
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media,
    std::optional<MicroPorosityParameters> const& micro_porosity_parameters,
    std::optional<VKPotentialExchangeParameters> const&
        vk_potential_exchange_parameters,
    std::map<int, VKPotentialExchangeParameters> const&
        vk_potential_exchange_parameters_by_material)
{
    namespace MPL = MaterialPropertyLib;

    bool const micro_porosity_enabled = micro_porosity_parameters.has_value();
    bool any_saturation_micro = false;
    bool const any_vk_enabled =
        (vk_potential_exchange_parameters &&
         vk_potential_exchange_parameters->enabled) ||
        std::any_of(vk_potential_exchange_parameters_by_material.begin(),
                    vk_potential_exchange_parameters_by_material.end(),
                    [](auto const& item) { return item.second.enabled; });

    for (auto const& [material_id, medium] : media)
    {
        bool const has_saturation_micro =
            medium->hasProperty(MPL::PropertyType::saturation_micro);
        any_saturation_micro = any_saturation_micro || has_saturation_micro;

        if (has_saturation_micro && !micro_porosity_enabled)
        {
            OGS_FATAL(
                "RichardsMechanics: medium {} defines 'saturation_micro' but "
                "the process has no <micro_porosity> block. Define "
                "<micro_porosity> or remove 'saturation_micro'.",
                material_id);
        }
    }

    if (micro_porosity_enabled && !any_saturation_micro && !any_vk_enabled)
    {
        OGS_FATAL(
            "RichardsMechanics: <micro_porosity> is configured, but no medium "
            "defines 'saturation_micro'. Define 'saturation_micro' in at least "
            "one medium or remove <micro_porosity>.");
    }

    if (any_vk_enabled && !micro_porosity_enabled)
    {
        OGS_FATAL(
            "RichardsMechanics: vk_potential_exchange.enabled=true requires "
            "a <micro_porosity> process block.");
    }

    for (auto const& [material_id, vkp] :
         vk_potential_exchange_parameters_by_material)
    {
        if (media.find(material_id) == media.end())
        {
            OGS_FATAL(
                "RichardsMechanics: vk_potential_exchange medium override "
                "references unknown material id {}.",
                material_id);
        }

    }
}

VKPotentialExchangeParameters parseVKPotentialExchangeParameters(
    BaseLib::ConfigTree const& config,
    std::optional<VKPotentialExchangeParameters> const& defaults,
    std::string const& context)
{
    auto const enabled =
        config.getConfigParameter<bool>("enabled",
                                        defaults ? defaults->enabled : false);

    auto const mode = parseVKPotentialExchangeMode(
        config.getConfigParameter<std::string>(
            "mode", defaults ? toString(defaults->mode) : "full_potential"));

    auto const pressure_tolerance = config.getConfigParameter<double>(
        "pressure_tolerance",
        defaults ? defaults->pressure_tolerance : 0.0);
    if (pressure_tolerance < 0.0)
    {
        OGS_FATAL(
            "RichardsMechanics: {} pressure_tolerance must be >= 0, got {:g}.",
            context, pressure_tolerance);
    }

    auto const micro_potential_convention = parseVKMicroPotentialConvention(
        config.getConfigParameter<std::string>(
            "micro_potential_convention",
            defaults ? toString(defaults->micro_potential_convention)
                     : "positive_reduced"));
    auto const local_nonlinear_solve_mode = parseVKLocalNonlinearSolveMode(
        config.getConfigParameter<std::string>(
            "local_nonlinear_solve_mode",
            defaults ? toString(defaults->local_nonlinear_solve_mode)
                     : "scalar_exchange"));
    auto const macro_porosity_update_mode = parseVKMacroPorosityUpdateMode(
        config.getConfigParameter<std::string>(
            "macro_porosity_update_mode",
            defaults ? toString(defaults->macro_porosity_update_mode)
                     : "algebraic_split"));
    auto const micro_solid_volume_fraction_mode =
        parseVKMicroSolidVolumeFractionMode(
            config.getConfigParameter<std::string>(
                "micro_solid_volume_fraction_mode",
                defaults
                    ? toString(defaults->micro_solid_volume_fraction_mode)
                    : "reference"));
    auto const potential_role_mapping =
        parseVKPotentialExchangeRoleMapping(
            config.getConfigParameter<std::string>(
                "potential_role_mapping",
                defaults ? toString(defaults->potential_role_mapping)
                         : "current_ogs"));

    auto get_positive_required_or_default =
        [&](char const* const key, double const fallback)
    {
        auto const value = config.getConfigParameterOptional<double>(key);
        double const selected = value ? *value : fallback;
        if (!(selected > 0.0))
        {
            OGS_FATAL(
                "RichardsMechanics: {} {} must be > 0, got {:g}.", context,
                key, selected);
        }
        return selected;
    };

    auto get_positive_optional_or_default =
        [&](char const* const key, std::optional<double> const fallback)
            -> std::optional<double>
    {
        auto const value = config.getConfigParameterOptional<double>(key);
        std::optional<double> selected = value ? std::optional<double>{*value}
                                               : fallback;
        if (selected && !(*selected > 0.0))
        {
            OGS_FATAL(
                "RichardsMechanics: {} {} must be > 0 if provided, got {:g}.",
                context, key, *selected);
        }
        return selected;
    };

    double const default_hamaker =
        defaults ? defaults->hamaker_constant : 0.0;
    double const default_surface =
        defaults ? defaults->specific_surface : 0.0;
    double const default_rho_sr =
        defaults ? defaults->micro_solid_density_reference : 0.0;
    double const default_ns =
        defaults ? defaults->micro_solid_volume_fraction_reference : 0.0;
    double const default_rho_l0 =
        defaults ? defaults->micro_liquid_density_reference : 0.0;
    double const default_a_rho =
        defaults ? defaults->micro_liquid_density_a : 0.0;
    double const default_b_rho =
        defaults ? defaults->micro_liquid_density_b : 0.0;

    double hamaker_constant = 0.0;
    double specific_surface = 0.0;
    double micro_solid_density_reference = 0.0;
    double micro_solid_volume_fraction_reference = 0.0;
    double micro_liquid_density_reference = 0.0;
    double micro_liquid_density_a = 0.0;
    double micro_liquid_density_b = 0.0;
    bool const uses_micro_liquid_density_eos =
        local_nonlinear_solve_mode ==
        VKLocalNonlinearSolveMode::ScalarNotebookMassStorage;

    if (enabled)
    {
        hamaker_constant =
            get_positive_required_or_default("hamaker_constant",
                                             default_hamaker);
        specific_surface =
            get_positive_required_or_default("specific_surface",
                                             default_surface);
        micro_solid_density_reference = get_positive_required_or_default(
            "micro_solid_density_reference", default_rho_sr);
        micro_solid_volume_fraction_reference =
            get_positive_required_or_default(
                "micro_solid_volume_fraction_reference", default_ns);
        if (uses_micro_liquid_density_eos)
        {
            micro_liquid_density_reference =
                get_positive_required_or_default(
                    "micro_liquid_density_reference", default_rho_l0);
            micro_liquid_density_a = get_positive_required_or_default(
                "micro_liquid_density_a", default_a_rho);
            micro_liquid_density_b = get_positive_required_or_default(
                "micro_liquid_density_b", default_b_rho);
        }
        else
        {
            micro_liquid_density_reference =
                get_positive_optional_or_default(
                    "micro_liquid_density_reference",
                    defaults ? std::optional<double>{
                                   defaults->micro_liquid_density_reference}
                             : std::nullopt)
                    .value_or(0.0);
            micro_liquid_density_a = get_positive_optional_or_default(
                                         "micro_liquid_density_a",
                                         defaults
                                             ? std::optional<double>{
                                                   defaults->micro_liquid_density_a}
                                             : std::nullopt)
                                         .value_or(0.0);
            micro_liquid_density_b = get_positive_optional_or_default(
                                         "micro_liquid_density_b",
                                         defaults
                                             ? std::optional<double>{
                                                   defaults->micro_liquid_density_b}
                                             : std::nullopt)
                                         .value_or(0.0);
        }
    }
    else
    {
        hamaker_constant = get_positive_optional_or_default(
                               "hamaker_constant",
                               defaults ? std::optional<double>{
                                              defaults->hamaker_constant}
                                        : std::nullopt)
                               .value_or(0.0);
        specific_surface = get_positive_optional_or_default(
                               "specific_surface",
                               defaults ? std::optional<double>{
                                              defaults->specific_surface}
                                        : std::nullopt)
                               .value_or(0.0);
        micro_solid_density_reference =
            get_positive_optional_or_default(
                "micro_solid_density_reference",
                defaults ? std::optional<double>{
                               defaults->micro_solid_density_reference}
                         : std::nullopt)
                .value_or(0.0);
        micro_solid_volume_fraction_reference =
            get_positive_optional_or_default(
                "micro_solid_volume_fraction_reference",
                defaults ? std::optional<double>{
                               defaults->micro_solid_volume_fraction_reference}
                         : std::nullopt)
                .value_or(0.0);
        micro_liquid_density_reference =
            get_positive_optional_or_default(
                "micro_liquid_density_reference",
                defaults ? std::optional<double>{
                               defaults->micro_liquid_density_reference}
                         : std::nullopt)
                .value_or(0.0);
        micro_liquid_density_a =
            get_positive_optional_or_default(
                "micro_liquid_density_a",
                defaults
                    ? std::optional<double>{defaults->micro_liquid_density_a}
                    : std::nullopt)
                .value_or(0.0);
        micro_liquid_density_b =
            get_positive_optional_or_default(
                "micro_liquid_density_b",
                defaults
                    ? std::optional<double>{defaults->micro_liquid_density_b}
                    : std::nullopt)
                .value_or(0.0);
    }

    auto const initial_micro_water_content = get_positive_optional_or_default(
        "initial_micro_water_content",
        defaults ? defaults->initial_micro_water_content : std::nullopt);

    auto const use_fd_jacobian_for_exchange = config.getConfigParameter<bool>(
        "fd_jacobian_for_exchange",
        defaults ? defaults->use_fd_jacobian_for_exchange : false);

    auto const fd_jacobian_perturbation = config.getConfigParameter<double>(
        "fd_jacobian_perturbation",
        defaults ? defaults->fd_jacobian_perturbation : 1e-8);
    if (!(fd_jacobian_perturbation > 0.0))
    {
        OGS_FATAL(
            "RichardsMechanics: {} fd_jacobian_perturbation must be > 0, got {:g}.",
            context, fd_jacobian_perturbation);
    }

    auto const check_local_jacobian = config.getConfigParameter<bool>(
        "check_local_jacobian",
        defaults ? defaults->check_local_jacobian : false);

    auto const local_jacobian_perturbation = config.getConfigParameter<double>(
        "local_jacobian_perturbation",
        defaults ? defaults->local_jacobian_perturbation : 1e-8);
    if (!(local_jacobian_perturbation > 0.0))
    {
        OGS_FATAL(
            "RichardsMechanics: {} local_jacobian_perturbation must be > 0, got {:g}.",
            context, local_jacobian_perturbation);
    }

    auto const local_jacobian_relative_tolerance =
        config.getConfigParameter<double>(
            "local_jacobian_relative_tolerance",
            defaults ? defaults->local_jacobian_relative_tolerance : 1e-3);
    if (!(local_jacobian_relative_tolerance >= 0.0))
    {
        OGS_FATAL(
            "RichardsMechanics: {} local_jacobian_relative_tolerance must be >= 0, got {:g}.",
            context, local_jacobian_relative_tolerance);
    }

    auto const vdw_relaxation_stress_gain = config.getConfigParameter<double>(
        "vdw_relaxation_stress_gain",
        defaults ? defaults->vdw_relaxation_stress_gain : 0.0);
    if (!(vdw_relaxation_stress_gain >= 0.0))
    {
        OGS_FATAL(
            "RichardsMechanics: {} vdw_relaxation_stress_gain must be >= 0, got {:g}.",
            context, vdw_relaxation_stress_gain);
    }

    auto const micro_water_content_stress_gain =
        config.getConfigParameter<double>(
            "micro_water_content_stress_gain",
            defaults ? defaults->micro_water_content_stress_gain : 0.0);
    if (!(micro_water_content_stress_gain >= 0.0))
    {
        OGS_FATAL(
            "RichardsMechanics: {} micro_water_content_stress_gain must be >= 0, got {:g}.",
            context, micro_water_content_stress_gain);
    }

    auto const micro_water_content_swelling_slope =
        config.getConfigParameter<double>(
            "micro_water_content_swelling_slope",
            defaults ? defaults->micro_water_content_swelling_slope : 0.0);
    if (!(micro_water_content_swelling_slope >= 0.0))
    {
        OGS_FATAL(
            "RichardsMechanics: {} micro_water_content_swelling_slope must be >= 0, got {:g}.",
            context, micro_water_content_swelling_slope);
    }

    return VKPotentialExchangeParameters{
        enabled,
        mode,
        pressure_tolerance,
        hamaker_constant,
        specific_surface,
        micro_solid_density_reference,
        micro_solid_volume_fraction_reference,
        micro_liquid_density_reference,
        micro_liquid_density_a,
        micro_liquid_density_b,
        micro_potential_convention,
        local_nonlinear_solve_mode,
        macro_porosity_update_mode,
        micro_solid_volume_fraction_mode,
        potential_role_mapping,
        initial_micro_water_content,
        use_fd_jacobian_for_exchange,
        fd_jacobian_perturbation,
        check_local_jacobian,
        local_jacobian_perturbation,
        local_jacobian_relative_tolerance,
        vdw_relaxation_stress_gain,
        micro_water_content_stress_gain,
        micro_water_content_swelling_slope};
}

template <int DisplacementDim>
std::unique_ptr<Process> createRichardsMechanicsProcess(
    std::string const& name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config,
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media)
{
    //! \ogs_file_param{prj__processes__process__type}
    config.checkConfigParameter("type", "RICHARDS_MECHANICS");
    DBUG("Create RichardsMechanicsProcess.");

    auto const coupling_scheme =
        //! \ogs_file_param{prj__processes__process__RICHARDS_MECHANICS__coupling_scheme}
        config.getConfigParameterOptional<std::string>("coupling_scheme");
    const bool use_monolithic_scheme =
        !(coupling_scheme && (*coupling_scheme == "staggered"));

    /// \section processvariablesrm Process Variables

    //! \ogs_file_param{prj__processes__process__RICHARDS_MECHANICS__process_variables}
    auto const pv_config = config.getConfigSubtree("process_variables");

    ProcessVariable* variable_p;
    ProcessVariable* variable_u;
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>
        process_variables;
    if (use_monolithic_scheme)  // monolithic scheme.
    {
        /// Primary process variables as they appear in the global component
        /// vector:
        auto per_process_variables = findProcessVariables(
            variables, pv_config,
            {//! \ogs_file_param_special{prj__processes__process__RICHARDS_MECHANICS__process_variables__pressure}
             "pressure",
             //! \ogs_file_param_special{prj__processes__process__RICHARDS_MECHANICS__process_variables__displacement}
             "displacement"});
        variable_p = &per_process_variables[0].get();
        variable_u = &per_process_variables[1].get();
        process_variables.push_back(std::move(per_process_variables));
    }
    else  // staggered scheme.
    {
        using namespace std::string_literals;
        for (auto const& variable_name : {"pressure"s, "displacement"s})
        {
            auto per_process_variables =
                findProcessVariables(variables, pv_config, {variable_name});
            process_variables.push_back(std::move(per_process_variables));
        }
        variable_p = &process_variables[0][0].get();
        variable_u = &process_variables[1][0].get();
    }

    DBUG("Associate displacement with process variable '{:s}'.",
         variable_u->getName());

    if (variable_u->getNumberOfGlobalComponents() != DisplacementDim)
    {
        OGS_FATAL(
            "Number of components of the process variable '{:s}' is different "
            "from the displacement dimension: got {:d}, expected {:d}",
            variable_u->getName(),
            variable_u->getNumberOfGlobalComponents(),
            DisplacementDim);
    }

    DBUG("Associate pressure with process variable '{:s}'.",
         variable_p->getName());
    if (variable_p->getNumberOfGlobalComponents() != 1)
    {
        OGS_FATAL(
            "Pressure process variable '{:s}' is not a scalar variable but has "
            "{:d} components.",
            variable_p->getName(),
            variable_p->getNumberOfGlobalComponents());
    }

    auto solid_constitutive_relations =
        MaterialLib::Solids::createConstitutiveRelations<DisplacementDim>(
            parameters, local_coordinate_system, materialIDs(mesh), config);

    /// \section parametersrm Process Parameters
    // Specific body force
    Eigen::Matrix<double, DisplacementDim, 1> specific_body_force;
    {
        std::vector<double> const b =
            //! \ogs_file_param{prj__processes__process__RICHARDS_MECHANICS__specific_body_force}
            config.getConfigParameter<std::vector<double>>(
                "specific_body_force");
        if (b.size() != DisplacementDim)
        {
            OGS_FATAL(
                "The size of the specific body force vector does not match the "
                "displacement dimension. Vector size is {:d}, displacement "
                "dimension is {:d}",
                b.size(), DisplacementDim);
        }

        std::copy_n(b.data(), b.size(), specific_body_force.data());
    }

    auto media_map =
        MaterialPropertyLib::createMaterialSpatialDistributionMap(media, mesh);
    DBUG("Check the media properties of RichardsMechanics process ...");
    checkMPLProperties(media);
    DBUG("Media properties verified.");

    // Initial stress conditions
    auto const initial_stress =
        ProcessLib::createInitialStress<DisplacementDim>(config, parameters,
                                                         mesh);

    std::optional<MicroPorosityParameters> micro_porosity_parameters;
    if (auto const micro_porosity_config =
            //! \ogs_file_param{prj__processes__process__RICHARDS_MECHANICS__micro_porosity}
        config.getConfigSubtreeOptional("micro_porosity"))
    {
        micro_porosity_parameters = MicroPorosityParameters{
            NumLib::createNewtonRaphsonSolverParameters(
                //! \ogs_file_param{prj__processes__process__RICHARDS_MECHANICS__micro_porosity__nonlinear_solver}
                micro_porosity_config->getConfigSubtree("nonlinear_solver")),
            //! \ogs_file_param{prj__processes__process__RICHARDS_MECHANICS__micro_porosity__mass_exchange_coefficient}
            micro_porosity_config->getConfigParameter<double>(
                "mass_exchange_coefficient")};
    }

    std::optional<VKPotentialExchangeParameters> vk_potential_exchange_parameters;
    std::map<int, VKPotentialExchangeParameters>
        vk_potential_exchange_parameters_by_material;
    if (auto const vk_potential_exchange_config =
            //! \ogs_file_param{prj__processes__process__RICHARDS_MECHANICS__vk_potential_exchange}
            config.getConfigSubtreeOptional("vk_potential_exchange"))
    {
        vk_potential_exchange_parameters = parseVKPotentialExchangeParameters(
            *vk_potential_exchange_config, std::nullopt,
            "vk_potential_exchange");

        for (auto medium_config :
             vk_potential_exchange_config->getConfigSubtreeList("medium"))
        {
            int const material_id = medium_config.getConfigAttribute<int>("id");
            if (!vk_potential_exchange_parameters_by_material
                     .emplace(material_id,
                              parseVKPotentialExchangeParameters(
                                  medium_config,
                                  vk_potential_exchange_parameters,
                                  fmt::format(
                                      "vk_potential_exchange medium id {}",
                                      material_id)))
                     .second)
            {
                OGS_FATAL(
                    "RichardsMechanics: duplicate vk_potential_exchange medium override for material id {}.",
                    material_id);
            }
        }
    }

    validateMicroPorosityAndVKConfiguration(
        media, micro_porosity_parameters, vk_potential_exchange_parameters,
        vk_potential_exchange_parameters_by_material);

    auto const mass_lumping =
        //! \ogs_file_param{prj__processes__process__RICHARDS_MECHANICS__mass_lumping}
        config.getConfigParameter<bool>("mass_lumping", false);

    auto const explicit_hm_coupling_in_unsaturated_zone =
        //! \ogs_file_param{prj__processes__process__RICHARDS_MECHANICS__explicit_hm_coupling_in_unsaturated_zone}
        config.getConfigParameter<bool>(
            "explicit_hm_coupling_in_unsaturated_zone", false);

    auto const is_linear =
        //! \ogs_file_param{prj__processes__process__linear}
        config.getConfigParameter("linear", false);

    bool const use_numerical_jacobian =
        jacobian_assembler->isPerturbationEnabled();

    logPhase0TransitionAudit<DisplacementDim>(
        name, use_monolithic_scheme, *variable_p, *variable_u, media,
        solid_constitutive_relations, micro_porosity_parameters,
        vk_potential_exchange_parameters,
        vk_potential_exchange_parameters_by_material);

    RichardsMechanicsProcessData<DisplacementDim> process_data{
        materialIDs(mesh),
        std::move(media_map),
        std::move(solid_constitutive_relations),
        initial_stress,
        specific_body_force,
        micro_porosity_parameters,
        vk_potential_exchange_parameters,
        vk_potential_exchange_parameters_by_material,
        mass_lumping,
        explicit_hm_coupling_in_unsaturated_zone,
        use_numerical_jacobian};

    SecondaryVariableCollection secondary_variables;

    ProcessLib::createSecondaryVariables(config, secondary_variables);

    return std::make_unique<RichardsMechanicsProcess<DisplacementDim>>(
        std::move(name), mesh, std::move(jacobian_assembler), parameters,
        integration_order, std::move(process_variables),
        std::move(process_data), std::move(secondary_variables),
        use_monolithic_scheme, is_linear);
}

template std::unique_ptr<Process> createRichardsMechanicsProcess<2>(
    std::string const& name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config,
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media);

template std::unique_ptr<Process> createRichardsMechanicsProcess<3>(
    std::string const& name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config,
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media);

}  // namespace RichardsMechanics
}  // namespace ProcessLib
