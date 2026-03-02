// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <string_view>

#include "BaseLib/StrongType.h"

namespace ProcessLib::RichardsMechanics
{
using VKMicroWaterContent =
    BaseLib::StrongType<double, struct VKMicroWaterContentTag>;

constexpr std::string_view ioName(struct VKMicroWaterContentTag*)
{
    return "vk_micro_water_content";
}
}  // namespace ProcessLib::RichardsMechanics

