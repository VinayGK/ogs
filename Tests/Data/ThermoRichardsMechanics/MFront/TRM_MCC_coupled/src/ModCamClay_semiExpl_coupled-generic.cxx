/*!
* \file   ModCamClay_semiExpl_coupled-generic.cxx
* \brief  This file implements the umat interface for the ModCamClay_semiExpl_coupled behaviour law
* \author Aqeel Afzal Chaudhry , Christian Silbermann , Eric Simo , Miguel Mánica , Thomas Helfer , Thomas Nagel
* \date   01 / 09 / 2025
*/

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif /* NOMINMAX */
#include <windows.h>
#ifdef small
#undef small
#endif /* small */
#endif /* _WIN32 */

#ifndef MFRONT_SHAREDOBJ
#define MFRONT_SHAREDOBJ TFEL_VISIBILITY_EXPORT
#endif /* MFRONT_SHAREDOBJ */

#ifndef MFRONT_EXPORT_SYMBOL
#define MFRONT_EXPORT_SYMBOL(TYPE, NAME, VALUE) \
  MFRONT_SHAREDOBJ extern TYPE NAME;            \
  MFRONT_SHAREDOBJ TYPE NAME = VALUE
#endif /* MFRONT_EXPORT_SYMBOL*/

#ifndef MFRONT_EXPORT_ARRAY_ARGUMENTS
#define MFRONT_EXPORT_ARRAY_ARGUMENTS(...) __VA_ARGS__
#endif /* MFRONT_EXPORT_ARRAY_ARGUMENTS */

#ifndef MFRONT_EXPORT_ARRAY_OF_SYMBOLS
#define MFRONT_EXPORT_ARRAY_OF_SYMBOLS(TYPE, NAME, SIZE, VALUE) \
  MFRONT_SHAREDOBJ extern TYPE NAME[SIZE];                      \
  MFRONT_SHAREDOBJ TYPE NAME[SIZE] = {VALUE}
#endif /* MFRONT_EXPORT_ARRAY_OF_SYMBOLS*/

#include<iostream>
#include<cstdlib>
#include"TFEL/Material/OutOfBoundsPolicy.hxx"
#include"TFEL/Math/t2tot2.hxx"
#include"TFEL/Math/t2tost2.hxx"
#include"TFEL/Material/ModCamClay_semiExpl_coupled.hxx"
#include"MFront/GenericBehaviour/GenericBehaviourTraits.hxx"
#include"MFront/GenericBehaviour/Integrate.hxx"
#include"MFront/GenericBehaviour/ModCamClay_semiExpl_coupled-generic.hxx"

static tfel::material::OutOfBoundsPolicy&
ModCamClay_semiExpl_coupled_getOutOfBoundsPolicy(){
static auto policy = tfel::material::None;
return policy;
}

namespace mfront::gb{

template<>
struct GenericBehaviourTraits<tfel::material::ModCamClay_semiExpl_coupled<tfel::material::ModellingHypothesis::AXISYMMETRICALGENERALISEDPLANESTRAIN, real, false>>{
static constexpr auto hypothesis = tfel::material::ModellingHypothesis::AXISYMMETRICALGENERALISEDPLANESTRAIN;
static constexpr auto N = tfel::material::ModellingHypothesisToSpaceDimension<hypothesis>::value;
static constexpr auto StensorSize = tfel::material::ModellingHypothesisToStensorSize<hypothesis>::value;
static constexpr auto TensorSize = tfel::material::ModellingHypothesisToTensorSize<hypothesis>::value;
};

template<>
struct GenericBehaviourTraits<tfel::material::ModCamClay_semiExpl_coupled<tfel::material::ModellingHypothesis::AXISYMMETRICALGENERALISEDPLANESTRESS, real, false>>{
static constexpr auto hypothesis = tfel::material::ModellingHypothesis::AXISYMMETRICALGENERALISEDPLANESTRESS;
static constexpr auto N = tfel::material::ModellingHypothesisToSpaceDimension<hypothesis>::value;
static constexpr auto StensorSize = tfel::material::ModellingHypothesisToStensorSize<hypothesis>::value;
static constexpr auto TensorSize = tfel::material::ModellingHypothesisToTensorSize<hypothesis>::value;
};

template<>
struct GenericBehaviourTraits<tfel::material::ModCamClay_semiExpl_coupled<tfel::material::ModellingHypothesis::AXISYMMETRICAL, real, false>>{
static constexpr auto hypothesis = tfel::material::ModellingHypothesis::AXISYMMETRICAL;
static constexpr auto N = tfel::material::ModellingHypothesisToSpaceDimension<hypothesis>::value;
static constexpr auto StensorSize = tfel::material::ModellingHypothesisToStensorSize<hypothesis>::value;
static constexpr auto TensorSize = tfel::material::ModellingHypothesisToTensorSize<hypothesis>::value;
};

template<>
struct GenericBehaviourTraits<tfel::material::ModCamClay_semiExpl_coupled<tfel::material::ModellingHypothesis::PLANESTRESS, real, false>>{
static constexpr auto hypothesis = tfel::material::ModellingHypothesis::PLANESTRESS;
static constexpr auto N = tfel::material::ModellingHypothesisToSpaceDimension<hypothesis>::value;
static constexpr auto StensorSize = tfel::material::ModellingHypothesisToStensorSize<hypothesis>::value;
static constexpr auto TensorSize = tfel::material::ModellingHypothesisToTensorSize<hypothesis>::value;
};

template<>
struct GenericBehaviourTraits<tfel::material::ModCamClay_semiExpl_coupled<tfel::material::ModellingHypothesis::PLANESTRAIN, real, false>>{
static constexpr auto hypothesis = tfel::material::ModellingHypothesis::PLANESTRAIN;
static constexpr auto N = tfel::material::ModellingHypothesisToSpaceDimension<hypothesis>::value;
static constexpr auto StensorSize = tfel::material::ModellingHypothesisToStensorSize<hypothesis>::value;
static constexpr auto TensorSize = tfel::material::ModellingHypothesisToTensorSize<hypothesis>::value;
};

template<>
struct GenericBehaviourTraits<tfel::material::ModCamClay_semiExpl_coupled<tfel::material::ModellingHypothesis::GENERALISEDPLANESTRAIN, real, false>>{
static constexpr auto hypothesis = tfel::material::ModellingHypothesis::GENERALISEDPLANESTRAIN;
static constexpr auto N = tfel::material::ModellingHypothesisToSpaceDimension<hypothesis>::value;
static constexpr auto StensorSize = tfel::material::ModellingHypothesisToStensorSize<hypothesis>::value;
static constexpr auto TensorSize = tfel::material::ModellingHypothesisToTensorSize<hypothesis>::value;
};

template<>
struct GenericBehaviourTraits<tfel::material::ModCamClay_semiExpl_coupled<tfel::material::ModellingHypothesis::TRIDIMENSIONAL, real, false>>{
static constexpr auto hypothesis = tfel::material::ModellingHypothesis::TRIDIMENSIONAL;
static constexpr auto N = tfel::material::ModellingHypothesisToSpaceDimension<hypothesis>::value;
static constexpr auto StensorSize = tfel::material::ModellingHypothesisToStensorSize<hypothesis>::value;
static constexpr auto TensorSize = tfel::material::ModellingHypothesisToTensorSize<hypothesis>::value;
};

} // end of namespace mfront::gb

#ifdef __cplusplus
extern "C"{
#endif /* __cplusplus */

MFRONT_EXPORT_SYMBOL(const char*, ModCamClay_semiExpl_coupled_author, "Aqeel Afzal Chaudhry , Christian Silbermann , Eric Simo , Miguel Mánica , Thomas Helfer , Thomas Nagel");

MFRONT_EXPORT_SYMBOL(const char*, ModCamClay_semiExpl_coupled_date, "01 / 09 / 2025");

MFRONT_EXPORT_SYMBOL(const char*, ModCamClay_semiExpl_coupled_description, "The modified cam - clay model according to Callari ( 1998 ) :\nA finite-strain cam-clay model in the framework of multiplicative elasto-plasticity\nbut here in a consistent geometrically linear form ( linearized volume ratio evolution )\nsemi - explicit due to explicit volume ratio update at the end of time step ,\nnonlinear hypoelastic behavior : pressure - dependent bulk modulus , constant Poisson ratio ,\nincremental formulation assuming constant elastic parameters over the time step ,\nnormalized plastic flow direction , lower limit for a minimal pre - consolidation pressure ,\nwith hydraulic and thermal coupling , saturation dependent swelling ,");

MFRONT_EXPORT_SYMBOL(const char*, ModCamClay_semiExpl_coupled_validator, "");

MFRONT_EXPORT_SYMBOL(const char*, ModCamClay_semiExpl_coupled_build_id, "");

MFRONT_EXPORT_SYMBOL(const char*, ModCamClay_semiExpl_coupled_mfront_ept, "ModCamClay_semiExpl_coupled");

MFRONT_EXPORT_SYMBOL(const char*, ModCamClay_semiExpl_coupled_tfel_version, "4.2.4-dev");

MFRONT_EXPORT_SYMBOL(const char*, ModCamClay_semiExpl_coupled_unit_system, "");

MFRONT_EXPORT_SYMBOL(unsigned short, ModCamClay_semiExpl_coupled_mfront_mkt, 1u);

MFRONT_EXPORT_SYMBOL(const char*, ModCamClay_semiExpl_coupled_mfront_interface, "Generic");

MFRONT_EXPORT_SYMBOL(const char*, ModCamClay_semiExpl_coupled_src, "ModCamClay_semiExpl_coupled.mfront");

MFRONT_EXPORT_SYMBOL(unsigned short, ModCamClay_semiExpl_coupled_nModellingHypotheses, 7u);

MFRONT_EXPORT_ARRAY_OF_SYMBOLS(const char *, ModCamClay_semiExpl_coupled_ModellingHypotheses, 7, MFRONT_EXPORT_ARRAY_ARGUMENTS("AxisymmetricalGeneralisedPlaneStrain",
"AxisymmetricalGeneralisedPlaneStress","Axisymmetrical","PlaneStress","PlaneStrain","GeneralisedPlaneStrain",
"Tridimensional"));

MFRONT_EXPORT_SYMBOL(unsigned short, ModCamClay_semiExpl_coupled_nMainVariables, 2u);

MFRONT_EXPORT_SYMBOL(unsigned short, ModCamClay_semiExpl_coupled_nGradients, 2u);

MFRONT_EXPORT_ARRAY_OF_SYMBOLS(int, ModCamClay_semiExpl_coupled_GradientsTypes, 2, MFRONT_EXPORT_ARRAY_ARGUMENTS(1,
0));

MFRONT_EXPORT_ARRAY_OF_SYMBOLS(const char *, ModCamClay_semiExpl_coupled_Gradients, 2, MFRONT_EXPORT_ARRAY_ARGUMENTS("Strain",
"LiquidPressure"));

MFRONT_EXPORT_SYMBOL(unsigned short, ModCamClay_semiExpl_coupled_nThermodynamicForces, 2u);

MFRONT_EXPORT_ARRAY_OF_SYMBOLS(int, ModCamClay_semiExpl_coupled_ThermodynamicForcesTypes, 2, MFRONT_EXPORT_ARRAY_ARGUMENTS(1,
0));

MFRONT_EXPORT_ARRAY_OF_SYMBOLS(const char *, ModCamClay_semiExpl_coupled_ThermodynamicForces, 2, MFRONT_EXPORT_ARRAY_ARGUMENTS("Stress",
"Saturation"));

MFRONT_EXPORT_SYMBOL(unsigned short, ModCamClay_semiExpl_coupled_nTangentOperatorBlocks, 8u);

MFRONT_EXPORT_ARRAY_OF_SYMBOLS(const char *, ModCamClay_semiExpl_coupled_TangentOperatorBlocks, 8, MFRONT_EXPORT_ARRAY_ARGUMENTS("Stress",
"Strain","Stress","Temperature","Stress","LiquidPressure",
"Saturation","LiquidPressure"));

MFRONT_EXPORT_SYMBOL(unsigned short, ModCamClay_semiExpl_coupled_BehaviourType, 0u);

MFRONT_EXPORT_SYMBOL(unsigned short, ModCamClay_semiExpl_coupled_BehaviourKinematic, 0u);

MFRONT_EXPORT_SYMBOL(unsigned short, ModCamClay_semiExpl_coupled_SymmetryType, 0u);

MFRONT_EXPORT_SYMBOL(unsigned short, ModCamClay_semiExpl_coupled_ElasticSymmetryType, 0u);

MFRONT_EXPORT_SYMBOL(unsigned short, ModCamClay_semiExpl_coupled_api_version, 1u);

MFRONT_EXPORT_SYMBOL(unsigned short, ModCamClay_semiExpl_coupled_TemperatureRemovedFromExternalStateVariables, 1u);

MFRONT_EXPORT_SYMBOL(unsigned short, ModCamClay_semiExpl_coupled_UsableInPurelyImplicitResolution, 0u);

MFRONT_EXPORT_SYMBOL(unsigned short, ModCamClay_semiExpl_coupled_nMaterialProperties, 18u);

MFRONT_EXPORT_ARRAY_OF_SYMBOLS(const char *, ModCamClay_semiExpl_coupled_MaterialProperties, 18, MFRONT_EXPORT_ARRAY_ARGUMENTS("PoissonRatio",
"CriticalStateLineSlope","SwellingLineSlope","VirginConsolidationLineSlope","CharacteristicPreConsolidationPressure","InitialVolumeRatio",
"ThermalExpansion","BiotCoefficient","ResidualLiquidSaturation","ResidualGasSaturation","BubblePressure",
"VanGenuchtenExponent_m","BishopsModelType","BishopsModelValue","SwellingPressures","Exponents",
"LowerSaturationLimit","UpperSaturationLimit"));

MFRONT_EXPORT_SYMBOL(unsigned short, ModCamClay_semiExpl_coupled_nInternalStateVariables, 6u);
MFRONT_EXPORT_ARRAY_OF_SYMBOLS(const char *, ModCamClay_semiExpl_coupled_InternalStateVariables, 6, MFRONT_EXPORT_ARRAY_ARGUMENTS("swelling_stress",
"ElasticStrain","EquivalentPlasticStrain","PreConsolidationPressure","PlasticVolumetricStrain","VolumeRatio"));

MFRONT_EXPORT_ARRAY_OF_SYMBOLS(int, ModCamClay_semiExpl_coupled_InternalStateVariablesTypes, 6, MFRONT_EXPORT_ARRAY_ARGUMENTS(1,1,0,0,0,0));

MFRONT_EXPORT_SYMBOL(unsigned short, ModCamClay_semiExpl_coupled_nExternalStateVariables, 0u);
MFRONT_EXPORT_SYMBOL(const char * const *, ModCamClay_semiExpl_coupled_ExternalStateVariables, nullptr);

MFRONT_EXPORT_SYMBOL(const int *, ModCamClay_semiExpl_coupled_ExternalStateVariablesTypes, nullptr);

MFRONT_EXPORT_SYMBOL(unsigned short, ModCamClay_semiExpl_coupled_nParameters, 6u);

MFRONT_EXPORT_ARRAY_OF_SYMBOLS(const char *, ModCamClay_semiExpl_coupled_Parameters, 6, MFRONT_EXPORT_ARRAY_ARGUMENTS("theta",
"epsilon","iterMax","minimal_time_step_scaling_factor","maximal_time_step_scaling_factor","numerical_jacobian_epsilon"));

MFRONT_EXPORT_ARRAY_OF_SYMBOLS(int, ModCamClay_semiExpl_coupled_ParametersTypes, 6, MFRONT_EXPORT_ARRAY_ARGUMENTS(0,0,2,0,0,0));

MFRONT_EXPORT_SYMBOL(double, ModCamClay_semiExpl_coupled_theta_ParameterDefaultValue, 1);

MFRONT_EXPORT_SYMBOL(double, ModCamClay_semiExpl_coupled_epsilon_ParameterDefaultValue, 1e-14);

MFRONT_EXPORT_SYMBOL(unsigned short, ModCamClay_semiExpl_coupled_iterMax_ParameterDefaultValue, 100);

MFRONT_EXPORT_SYMBOL(double, ModCamClay_semiExpl_coupled_minimal_time_step_scaling_factor_ParameterDefaultValue, 0.1);

MFRONT_EXPORT_SYMBOL(double, ModCamClay_semiExpl_coupled_maximal_time_step_scaling_factor_ParameterDefaultValue, 1.7976931348623e+308);

MFRONT_EXPORT_SYMBOL(double, ModCamClay_semiExpl_coupled_numerical_jacobian_epsilon_ParameterDefaultValue, 1e-15);

MFRONT_EXPORT_SYMBOL(unsigned short, ModCamClay_semiExpl_coupled_requiresStiffnessTensor, 0u);

MFRONT_EXPORT_SYMBOL(unsigned short, ModCamClay_semiExpl_coupled_requiresThermalExpansionCoefficientTensor, 0u);

MFRONT_EXPORT_SYMBOL(unsigned short, ModCamClay_semiExpl_coupled_nInitializeFunctions, 0u);

MFRONT_EXPORT_SYMBOL(const char * const *, ModCamClay_semiExpl_coupled_InitializeFunctions, nullptr);


MFRONT_EXPORT_SYMBOL(unsigned short, ModCamClay_semiExpl_coupled_nPostProcessings, 0u);

MFRONT_EXPORT_SYMBOL(const char * const *, ModCamClay_semiExpl_coupled_PostProcessings, nullptr);


MFRONT_EXPORT_SYMBOL(unsigned short, ModCamClay_semiExpl_coupled_ComputesInternalEnergy, 0u);

MFRONT_EXPORT_SYMBOL(unsigned short, ModCamClay_semiExpl_coupled_ComputesDissipatedEnergy, 0u);

MFRONT_SHAREDOBJ void
ModCamClay_semiExpl_coupled_setOutOfBoundsPolicy(const int p){
if(p==0){
ModCamClay_semiExpl_coupled_getOutOfBoundsPolicy() = tfel::material::None;
} else if(p==1){
ModCamClay_semiExpl_coupled_getOutOfBoundsPolicy() = tfel::material::Warning;
} else if(p==2){
ModCamClay_semiExpl_coupled_getOutOfBoundsPolicy() = tfel::material::Strict;
} else {
std::cerr << "ModCamClay_semiExpl_coupled_setOutOfBoundsPolicy: invalid argument\n";
}
}

MFRONT_SHAREDOBJ int
ModCamClay_semiExpl_coupled_setParameter(const char *const key,const double value){
using tfel::material::ModCamClay_semiExpl_coupledParametersInitializer;
auto& i = ModCamClay_semiExpl_coupledParametersInitializer::get();
try{
i.set(key,value);
} catch(std::runtime_error& e){
std::cerr << e.what() << std::endl;
return 0;
}
return 1;
}

MFRONT_SHAREDOBJ int
ModCamClay_semiExpl_coupled_setUnsignedShortParameter(const char *const key,const unsigned short value){
using tfel::material::ModCamClay_semiExpl_coupledParametersInitializer;
auto& i = ModCamClay_semiExpl_coupledParametersInitializer::get();
try{
i.set(key,value);
} catch(std::runtime_error& e){
std::cerr << e.what() << std::endl;
return 0;
}
return 1;
}

MFRONT_SHAREDOBJ int ModCamClay_semiExpl_coupled_AxisymmetricalGeneralisedPlaneStrain(mfront_gb_BehaviourData* const d){
using namespace tfel::material;
using real = mfront::gb::real;
constexpr auto h = ModellingHypothesis::AXISYMMETRICALGENERALISEDPLANESTRAIN;
using Behaviour = ModCamClay_semiExpl_coupled<h,real,false>;
const auto r = mfront::gb::integrate<Behaviour>(*d, Behaviour::STANDARDTANGENTOPERATOR, ModCamClay_semiExpl_coupled_getOutOfBoundsPolicy());
return r;
} // end of ModCamClay_semiExpl_coupled_AxisymmetricalGeneralisedPlaneStrain

MFRONT_SHAREDOBJ int ModCamClay_semiExpl_coupled_AxisymmetricalGeneralisedPlaneStress(mfront_gb_BehaviourData* const d){
using namespace tfel::material;
using real = mfront::gb::real;
constexpr auto h = ModellingHypothesis::AXISYMMETRICALGENERALISEDPLANESTRESS;
using Behaviour = ModCamClay_semiExpl_coupled<h,real,false>;
const auto r = mfront::gb::integrate<Behaviour>(*d, Behaviour::STANDARDTANGENTOPERATOR, ModCamClay_semiExpl_coupled_getOutOfBoundsPolicy());
return r;
} // end of ModCamClay_semiExpl_coupled_AxisymmetricalGeneralisedPlaneStress

MFRONT_SHAREDOBJ int ModCamClay_semiExpl_coupled_Axisymmetrical(mfront_gb_BehaviourData* const d){
using namespace tfel::material;
using real = mfront::gb::real;
constexpr auto h = ModellingHypothesis::AXISYMMETRICAL;
using Behaviour = ModCamClay_semiExpl_coupled<h,real,false>;
const auto r = mfront::gb::integrate<Behaviour>(*d, Behaviour::STANDARDTANGENTOPERATOR, ModCamClay_semiExpl_coupled_getOutOfBoundsPolicy());
return r;
} // end of ModCamClay_semiExpl_coupled_Axisymmetrical

MFRONT_SHAREDOBJ int ModCamClay_semiExpl_coupled_PlaneStress(mfront_gb_BehaviourData* const d){
using namespace tfel::material;
using real = mfront::gb::real;
constexpr auto h = ModellingHypothesis::PLANESTRESS;
using Behaviour = ModCamClay_semiExpl_coupled<h,real,false>;
const auto r = mfront::gb::integrate<Behaviour>(*d, Behaviour::STANDARDTANGENTOPERATOR, ModCamClay_semiExpl_coupled_getOutOfBoundsPolicy());
return r;
} // end of ModCamClay_semiExpl_coupled_PlaneStress

MFRONT_SHAREDOBJ int ModCamClay_semiExpl_coupled_PlaneStrain(mfront_gb_BehaviourData* const d){
using namespace tfel::material;
using real = mfront::gb::real;
constexpr auto h = ModellingHypothesis::PLANESTRAIN;
using Behaviour = ModCamClay_semiExpl_coupled<h,real,false>;
const auto r = mfront::gb::integrate<Behaviour>(*d, Behaviour::STANDARDTANGENTOPERATOR, ModCamClay_semiExpl_coupled_getOutOfBoundsPolicy());
return r;
} // end of ModCamClay_semiExpl_coupled_PlaneStrain

MFRONT_SHAREDOBJ int ModCamClay_semiExpl_coupled_GeneralisedPlaneStrain(mfront_gb_BehaviourData* const d){
using namespace tfel::material;
using real = mfront::gb::real;
constexpr auto h = ModellingHypothesis::GENERALISEDPLANESTRAIN;
using Behaviour = ModCamClay_semiExpl_coupled<h,real,false>;
const auto r = mfront::gb::integrate<Behaviour>(*d, Behaviour::STANDARDTANGENTOPERATOR, ModCamClay_semiExpl_coupled_getOutOfBoundsPolicy());
return r;
} // end of ModCamClay_semiExpl_coupled_GeneralisedPlaneStrain

MFRONT_SHAREDOBJ int ModCamClay_semiExpl_coupled_Tridimensional(mfront_gb_BehaviourData* const d){
using namespace tfel::material;
using real = mfront::gb::real;
constexpr auto h = ModellingHypothesis::TRIDIMENSIONAL;
using Behaviour = ModCamClay_semiExpl_coupled<h,real,false>;
const auto r = mfront::gb::integrate<Behaviour>(*d, Behaviour::STANDARDTANGENTOPERATOR, ModCamClay_semiExpl_coupled_getOutOfBoundsPolicy());
return r;
} // end of ModCamClay_semiExpl_coupled_Tridimensional

#ifdef __cplusplus
}
#endif /* __cplusplus */

