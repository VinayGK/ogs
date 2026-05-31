// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <limits>

#include "MaterialLib/MPL/Properties/CapillaryPressureSaturation/CreateSaturationTuller.h"
#include "MaterialLib/MPL/Properties/CapillaryPressureSaturation/SaturationTuller.h"
#include "Tests/MaterialLib/TestMPL.h"

namespace MPL = MaterialPropertyLib;

TEST(MaterialPropertyLib, SaturationTuller)
{
    double const S_L_res = 0.05;
    double const S_L_max = 1.0;
    double const An = 1.0;
    double const F_gamma = 0.858407;
    double const L = 1e-5;
    double const sigma = 0.072;
    double const p_tol = 1.0;
    double const p_cav = std::numeric_limits<double>::infinity();

    MPL::Property const& saturation = MPL::SaturationTuller{
        "saturation", S_L_res, S_L_max, An, F_gamma, L, sigma, p_tol, p_cav};

    MPL::VariableArray vars;
    ParameterLib::SpatialPosition const pos;
    double const t = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();

    // Saturated branch: p_c <= p_tol
    vars.capillary_pressure = 0.0;
    ASSERT_DOUBLE_EQ(saturation.template value<double>(vars, pos, t, dt), S_L_max);
    ASSERT_DOUBLE_EQ(
        saturation.template dValue<double>(
            vars, MPL::Variable::capillary_pressure, pos, t, dt),
        0.0);

    // Unsaturated branch monotonicity and derivative checks.
    double const p1 = 5e4;
    double const p2 = 1e5;
    vars.capillary_pressure = p1;
    double const S1 = saturation.template value<double>(vars, pos, t, dt);
    vars.capillary_pressure = p2;
    double const S2 = saturation.template value<double>(vars, pos, t, dt);
    ASSERT_GT(S1, S2);
    ASSERT_GE(S2, S_L_res);
    ASSERT_LE(S1, S_L_max);

    // Compare analytical derivatives against finite differences away from the branch.
    for (double const p_cap : {5e3, 1e4, 5e4, 1e5, 2e5})
    {
        vars.capillary_pressure = p_cap;
        double const S = saturation.template value<double>(vars, pos, t, dt);
        double const dS = saturation.template dValue<double>(
            vars, MPL::Variable::capillary_pressure, pos, t, dt);
        double const d2S = saturation.template d2Value<double>(
            vars, MPL::Variable::capillary_pressure,
            MPL::Variable::capillary_pressure, pos, t, dt);

        double const eps = 1e-3 * p_cap;
        vars.capillary_pressure = p_cap - eps;
        double const S_minus = saturation.template value<double>(vars, pos, t, dt);
        vars.capillary_pressure = p_cap + eps;
        double const S_plus = saturation.template value<double>(vars, pos, t, dt);

        double const dS_fd = (S_plus - S_minus) / (2.0 * eps);
        double const d2S_fd = (S_plus - 2.0 * S + S_minus) / (eps * eps);

        EXPECT_NEAR(dS, dS_fd, 1e-10);
        EXPECT_NEAR(d2S, d2S_fd, 1e-12);
    }
}

TEST(MaterialPropertyLib, SaturationTullerCavitationCutoff)
{
    double const S_L_res = 0.05;
    double const S_L_max = 1.0;
    double const An = 1.0;
    double const F_gamma = 0.858407;
    double const L = 1e-5;
    double const sigma = 0.072;
    double const p_tol = 1.0;
    double const p_cav = 1e5;

    // Reference curve without cutoff.
    MPL::Property const& uncut = MPL::SaturationTuller{
        "uncut",  S_L_res, S_L_max, An, F_gamma, L, sigma, p_tol,
        std::numeric_limits<double>::infinity()};
    // Same curve, cut at the cavitation pressure.
    MPL::Property const& cut = MPL::SaturationTuller{
        "cut", S_L_res, S_L_max, An, F_gamma, L, sigma, p_tol, p_cav};

    MPL::VariableArray vars;
    ParameterLib::SpatialPosition const pos;
    double const t = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();

    auto value = [&](MPL::Property const& p, double p_cap)
    {
        vars.capillary_pressure = p_cap;
        return p.template value<double>(vars, pos, t, dt);
    };
    auto dvalue = [&](MPL::Property const& p, double p_cap)
    {
        vars.capillary_pressure = p_cap;
        return p.template dValue<double>(
            vars, MPL::Variable::capillary_pressure, pos, t, dt);
    };

    // Below cavitation: identical to the uncut curve.
    EXPECT_DOUBLE_EQ(value(cut, 5e4), value(uncut, 5e4));
    EXPECT_DOUBLE_EQ(dvalue(cut, 5e4), dvalue(uncut, 5e4));

    // Continuity at the cavitation pressure (C^0).
    double const S_at_cav = value(uncut, p_cav);
    EXPECT_DOUBLE_EQ(value(cut, p_cav), S_at_cav);

    // Beyond cavitation: frozen saturation, zero storage response.
    EXPECT_DOUBLE_EQ(value(cut, 2.0 * p_cav), S_at_cav);
    EXPECT_DOUBLE_EQ(value(cut, 1e7), S_at_cav);
    EXPECT_DOUBLE_EQ(dvalue(cut, 2.0 * p_cav), 0.0);
    vars.capillary_pressure = 2.0 * p_cav;
    EXPECT_DOUBLE_EQ(
        cut.template d2Value<double>(vars, MPL::Variable::capillary_pressure,
                                     MPL::Variable::capillary_pressure, pos, t,
                                     dt),
        0.0);
}

TEST(MaterialPropertyLib, CreateSaturationTuller)
{
    char const xml[] =
        "<property>"
        "  <name>saturation</name>"
        "  <type>TullerRetention</type>"
        "  <residual_liquid_saturation>0.05</residual_liquid_saturation>"
        "  <maximum_liquid_saturation>1.0</maximum_liquid_saturation>"
        "  <area_factor_tuller>1.0</area_factor_tuller>"
        "  <pore_area_shapefactor_tuller>0.858407</pore_area_shapefactor_tuller>"
        "  <characteristic_pore_size>1e-5</characteristic_pore_size>"
        "  <surface_tension>0.072</surface_tension>"
        "  <pressure_tolerance>1.0</pressure_tolerance>"
        "</property>";

    std::unique_ptr<MPL::Property> const saturation =
        Tests::createTestProperty(xml, MPL::createSaturationTuller);

    MPL::VariableArray vars;
    ParameterLib::SpatialPosition const pos;
    vars.capillary_pressure = 0.0;
    auto const S = saturation->template value<double>(
        vars, pos, std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN());
    EXPECT_DOUBLE_EQ(S, 1.0);
}

