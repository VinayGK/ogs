/*!
* \file   ModCamClay_semiExpl_coupled-generic.hxx
* \brief  This file declares the umat interface for the ModCamClay_semiExpl_coupled behaviour law
* \author Aqeel Afzal Chaudhry , Christian Silbermann , Eric Simo , Miguel Mánica , Thomas Helfer , Thomas Nagel
* \date   01 / 09 / 2025
*/

#ifndef LIB_GENERIC_MODCAMCLAY_SEMIEXPL_COUPLED_HXX
#define LIB_GENERIC_MODCAMCLAY_SEMIEXPL_COUPLED_HXX

#include"TFEL/Config/TFELConfig.hxx"
#include"MFront/GenericBehaviour/BehaviourData.h"

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

#ifdef __cplusplus
extern "C"{
#endif /* __cplusplus */

MFRONT_SHAREDOBJ void
ModCamClay_semiExpl_coupled_setOutOfBoundsPolicy(const int);

MFRONT_SHAREDOBJ int
ModCamClay_semiExpl_coupled_setParameter(const char *const,const double);

MFRONT_SHAREDOBJ int
ModCamClay_semiExpl_coupled_setUnsignedShortParameter(const char *const,const unsigned short);

/*!
 * \param[in,out] d: material data
 */
MFRONT_SHAREDOBJ int ModCamClay_semiExpl_coupled_AxisymmetricalGeneralisedPlaneStrain(mfront_gb_BehaviourData* const);

/*!
 * \param[in,out] d: material data
 */
MFRONT_SHAREDOBJ int ModCamClay_semiExpl_coupled_AxisymmetricalGeneralisedPlaneStress(mfront_gb_BehaviourData* const);

/*!
 * \param[in,out] d: material data
 */
MFRONT_SHAREDOBJ int ModCamClay_semiExpl_coupled_Axisymmetrical(mfront_gb_BehaviourData* const);

/*!
 * \param[in,out] d: material data
 */
MFRONT_SHAREDOBJ int ModCamClay_semiExpl_coupled_PlaneStress(mfront_gb_BehaviourData* const);

/*!
 * \param[in,out] d: material data
 */
MFRONT_SHAREDOBJ int ModCamClay_semiExpl_coupled_PlaneStrain(mfront_gb_BehaviourData* const);

/*!
 * \param[in,out] d: material data
 */
MFRONT_SHAREDOBJ int ModCamClay_semiExpl_coupled_GeneralisedPlaneStrain(mfront_gb_BehaviourData* const);

/*!
 * \param[in,out] d: material data
 */
MFRONT_SHAREDOBJ int ModCamClay_semiExpl_coupled_Tridimensional(mfront_gb_BehaviourData* const);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* LIB_GENERIC_MODCAMCLAY_SEMIEXPL_COUPLED_HXX */
