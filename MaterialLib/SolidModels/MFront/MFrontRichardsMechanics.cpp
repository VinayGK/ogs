/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "MFrontRichardsMechanics.h"

namespace MaterialLib::Solids::MFront
{
template class MFrontRichardsMechanics<2>;
template class MFrontRichardsMechanics<3>;
}  // namespace MaterialLib::Solids::MFront
