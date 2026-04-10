// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <string_view>

#include "BaseLib/StrongType.h"

namespace ProcessLib::RichardsMechanics
{
using VKMicroPorosity = BaseLib::StrongType<double, struct VKMicroPorosityTag>;

constexpr std::string_view ioName(struct VKMicroPorosityTag*)
{
    return "micro_porosity";
}
}  // namespace ProcessLib::RichardsMechanics
