// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <string_view>

#include "BaseLib/StrongType.h"

namespace ProcessLib::RichardsMechanics
{
using VKMicroLiquidDensity =
    BaseLib::StrongType<double, struct VKMicroLiquidDensityTag>;

constexpr std::string_view ioName(struct VKMicroLiquidDensityTag*)
{
    return "micro_liquid_density";
}
}  // namespace ProcessLib::RichardsMechanics
