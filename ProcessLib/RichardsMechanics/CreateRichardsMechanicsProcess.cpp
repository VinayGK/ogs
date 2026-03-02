// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "CreateRichardsMechanicsProcess.h"

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
        vk_potential_exchange_parameters)
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
            "[RM Phase0 audit] VK potential-exchange config block: PRESENT (enabled={}, mode='{}', pressure_tolerance={} Pa, hamaker_constant={}, specific_surface={}, rho_SR_ref={}, n_S_ref={}, initial_n_l={}, fd_jacobian_for_exchange={}, fd_jacobian_perturbation={} ).",
            vkp.enabled ? "true" : "false", toString(vkp.mode),
            vkp.pressure_tolerance, vkp.hamaker_constant, vkp.specific_surface,
            vkp.micro_solid_density_reference,
            vkp.micro_solid_volume_fraction_reference,
            vkp.initial_micro_water_content
                ? std::to_string(*vkp.initial_micro_water_content)
                : std::string{"<unset>"},
            vkp.use_fd_jacobian_for_exchange ? "true" : "false",
            vkp.fd_jacobian_perturbation);
    }
    else
    {
        INFO("[RM Phase0 audit] VK potential-exchange config block: ABSENT.");
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
        vk_potential_exchange_parameters)
{
    namespace MPL = MaterialPropertyLib;

    bool const micro_porosity_enabled = micro_porosity_parameters.has_value();
    bool any_saturation_micro = false;

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

    if (micro_porosity_enabled && !any_saturation_micro)
    {
        OGS_FATAL(
            "RichardsMechanics: <micro_porosity> is configured, but no medium "
            "defines 'saturation_micro'. Define 'saturation_micro' in at least "
            "one medium or remove <micro_porosity>.");
    }

    bool const vk_enabled = vk_potential_exchange_parameters &&
                            vk_potential_exchange_parameters->enabled;
    if (vk_enabled && (!micro_porosity_enabled || !any_saturation_micro))
    {
        OGS_FATAL(
            "RichardsMechanics: vk_potential_exchange.enabled=true requires "
            "both a <micro_porosity> process block and medium property "
            "'saturation_micro'.");
    }
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
    if (auto const vk_potential_exchange_config =
            //! \ogs_file_param{prj__processes__process__RICHARDS_MECHANICS__vk_potential_exchange}
            config.getConfigSubtreeOptional("vk_potential_exchange"))
    {
        auto const enabled =
            //! \ogs_file_param{prj__processes__process__RICHARDS_MECHANICS__vk_potential_exchange__enabled}
            vk_potential_exchange_config->getConfigParameter<bool>("enabled",
                                                                   false);

        auto const mode = parseVKPotentialExchangeMode(
            //! \ogs_file_param{prj__processes__process__RICHARDS_MECHANICS__vk_potential_exchange__mode}
            vk_potential_exchange_config->getConfigParameter<std::string>(
                "mode", "full_potential"));

        auto const pressure_tolerance =
            //! \ogs_file_param{prj__processes__process__RICHARDS_MECHANICS__vk_potential_exchange__pressure_tolerance}
            vk_potential_exchange_config->getConfigParameter<double>(
                "pressure_tolerance", 0.0);

        if (pressure_tolerance < 0.0)
        {
            OGS_FATAL(
                "RichardsMechanics: vk_potential_exchange.pressure_tolerance "
                "must be >= 0, got {:g}.",
                pressure_tolerance);
        }

        auto get_positive_required = [&](char const* const key)
        {
            double const value =
                vk_potential_exchange_config->getConfigParameter<double>(key);
            if (!(value > 0.0))
            {
                OGS_FATAL(
                    "RichardsMechanics: vk_potential_exchange.{} must be > 0, "
                    "got {:g}.",
                    key, value);
            }
            return value;
        };

        auto get_positive_optional = [&](char const* const key)
            -> std::optional<double>
        {
            auto const value =
                vk_potential_exchange_config->getConfigParameterOptional<double>(
                    key);
            if (value && !(*value > 0.0))
            {
                OGS_FATAL(
                    "RichardsMechanics: vk_potential_exchange.{} must be > 0 "
                    "if provided, got {:g}.",
                    key, *value);
            }
            return value;
        };

        // Parse required values when the opt-in mode is enabled. If the block
        // is present but disabled, keep values optional to allow staged
        // configuration in project files without forcing full parameter input.
        double hamaker_constant = 0.0;
        double specific_surface = 0.0;
        double micro_solid_density_reference = 0.0;
        double micro_solid_volume_fraction_reference = 0.0;

        if (enabled)
        {
            hamaker_constant = get_positive_required("hamaker_constant");
            specific_surface = get_positive_required("specific_surface");
            micro_solid_density_reference =
                get_positive_required("micro_solid_density_reference");
            micro_solid_volume_fraction_reference =
                get_positive_required("micro_solid_volume_fraction_reference");
        }
        else
        {
            if (auto const v = get_positive_optional("hamaker_constant"))
            {
                hamaker_constant = *v;
            }
            if (auto const v = get_positive_optional("specific_surface"))
            {
                specific_surface = *v;
            }
            if (auto const v =
                    get_positive_optional("micro_solid_density_reference"))
            {
                micro_solid_density_reference = *v;
            }
            if (auto const v = get_positive_optional(
                    "micro_solid_volume_fraction_reference"))
            {
                micro_solid_volume_fraction_reference = *v;
            }
        }

        auto const initial_micro_water_content =
            //! \ogs_file_param{prj__processes__process__RICHARDS_MECHANICS__vk_potential_exchange__initial_micro_water_content}
            get_positive_optional("initial_micro_water_content");

        auto const use_fd_jacobian_for_exchange =
            //! \ogs_file_param{prj__processes__process__RICHARDS_MECHANICS__vk_potential_exchange__fd_jacobian_for_exchange}
            vk_potential_exchange_config->getConfigParameter<bool>(
                "fd_jacobian_for_exchange", false);

        auto const fd_jacobian_perturbation =
            //! \ogs_file_param{prj__processes__process__RICHARDS_MECHANICS__vk_potential_exchange__fd_jacobian_perturbation}
            vk_potential_exchange_config->getConfigParameter<double>(
                "fd_jacobian_perturbation", 1e-8);
        if (!(fd_jacobian_perturbation > 0.0))
        {
            OGS_FATAL(
                "RichardsMechanics: vk_potential_exchange.fd_jacobian_perturbation "
                "must be > 0, got {:g}.",
                fd_jacobian_perturbation);
        }

        vk_potential_exchange_parameters = VKPotentialExchangeParameters{
            enabled,
            mode,
            pressure_tolerance,
            hamaker_constant,
            specific_surface,
            micro_solid_density_reference,
            micro_solid_volume_fraction_reference,
            initial_micro_water_content,
            use_fd_jacobian_for_exchange,
            fd_jacobian_perturbation};
    }

    validateMicroPorosityAndVKConfiguration(
        media, micro_porosity_parameters, vk_potential_exchange_parameters);

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
        vk_potential_exchange_parameters);

    RichardsMechanicsProcessData<DisplacementDim> process_data{
        materialIDs(mesh),
        std::move(media_map),
        std::move(solid_constitutive_relations),
        initial_stress,
        specific_body_force,
        micro_porosity_parameters,
        vk_potential_exchange_parameters,
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
