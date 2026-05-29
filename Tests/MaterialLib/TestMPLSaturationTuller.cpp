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

    MPL::Property const& saturation = MPL::SaturationTuller{
        "saturation", S_L_res, S_L_max, An, F_gamma, L, sigma, p_tol};

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

// Option B film branch OFF (macro_specific_surface = 0) must reproduce the
// pure-Tuller relation exactly, point-for-point, in value and both
// derivatives. Physics anchor: regression equivalence (the film term is the
// only difference and is gated to zero). No new expected-value literals.
TEST(MaterialPropertyLib, SaturationTullerFilmOffRecovery)
{
    double const S_L_res = 0.05;
    double const S_L_max = 1.0;
    double const An = 1.0;
    double const F_gamma = 0.858407;
    double const L = 1e-5;
    double const sigma = 0.072;
    double const p_tol = 1.0;

    MPL::SaturationTuller const pure_obj{
        "saturation", S_L_res, S_L_max, An, F_gamma, L, sigma, p_tol};
    // Explicit film-off: macro_specific_surface = 0 disables the branch.
    MPL::SaturationTuller const film_off_obj{
        "saturation", S_L_res, S_L_max, An, F_gamma, L, sigma, p_tol,
        0.0, 0.0, 0.0};

    MPL::Property const& pure = pure_obj;
    MPL::Property const& film_off = film_off_obj;

    MPL::VariableArray vars;
    ParameterLib::SpatialPosition const pos;
    double const t = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();

    for (double const p_cap : {0.0, 5e3, 1e4, 5e4, 1e5, 2e5})
    {
        vars.capillary_pressure = p_cap;
        EXPECT_DOUBLE_EQ(pure.template value<double>(vars, pos, t, dt),
                         film_off.template value<double>(vars, pos, t, dt));
        EXPECT_DOUBLE_EQ(
            pure.template dValue<double>(
                vars, MPL::Variable::capillary_pressure, pos, t, dt),
            film_off.template dValue<double>(
                vars, MPL::Variable::capillary_pressure, pos, t, dt));
        EXPECT_DOUBLE_EQ(
            pure.template d2Value<double>(
                vars, MPL::Variable::capillary_pressure,
                MPL::Variable::capillary_pressure, pos, t, dt),
            film_off.template d2Value<double>(
                vars, MPL::Variable::capillary_pressure,
                MPL::Variable::capillary_pressure, pos, t, dt));
    }
}

// Option B film branch ON: the summed WRC stays monotone decreasing and its
// analytic derivatives match central finite differences. Physics anchors:
// (d) monotonicity and analytic-derivative consistency. No physical magnitude
// is asserted, so the film inputs below are fixtures, not material claims.
TEST(MaterialPropertyLib, SaturationTullerFilm)
{
    double const S_L_res = 0.05;
    double const S_L_max = 1.0;
    double const An = 1.0;
    double const F_gamma = 0.858407;
    double const L = 1e-5;
    double const sigma = 0.072;
    double const p_tol = 1.0;

    // A_Hamaker = 2.2e-20 J: Israelachvili & Adams 1978, J. Chem. Soc. Faraday
    // Trans. I, vol. 74, p. 975, Table 2 (mica-water-mica SFA) -- the
    // established DSM Hamaker anchor (cf. PotentialExchange.h). Not a knob.
    double const A_H = 2.2e-20;
    // a_v, phi_M: SYNTHETIC test fixtures, not physical values. Chosen so the
    // film coefficient C_film = (a_v/phi_M)*(A_H/(6 pi))^(1/3) is O(1e-2),
    // keeping corner+film unclamped across the sampled high-suction window so
    // the consistency/monotonicity assertions are well posed. A physical
    // (FEBEX-sourced) a_v is required before any PRJ uses the film branch
    // (CLAUDE.md §1.1/§12).
    double const a_v = 3.0e4;   // [1/m]
    double const phi_M = 0.3;   // [-]

    MPL::SaturationTuller const film_obj{
        "saturation", S_L_res, S_L_max, An, F_gamma, L, sigma, p_tol,
        a_v, A_H, phi_M};
    MPL::Property const& film = film_obj;

    MPL::VariableArray vars;
    ParameterLib::SpatialPosition const pos;
    double const t = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();

    // Saturated branch is unaffected by the film term.
    vars.capillary_pressure = 0.0;
    EXPECT_DOUBLE_EQ(film.template value<double>(vars, pos, t, dt), S_L_max);

    // High-suction window where corner+film stays below S_L_max (unclamped).
    double S_prev = std::numeric_limits<double>::infinity();
    for (double const p_cap : {1e5, 2e5, 5e5, 1e6, 2e6})
    {
        vars.capillary_pressure = p_cap;
        double const S = film.template value<double>(vars, pos, t, dt);
        double const dS = film.template dValue<double>(
            vars, MPL::Variable::capillary_pressure, pos, t, dt);
        double const d2S = film.template d2Value<double>(
            vars, MPL::Variable::capillary_pressure,
            MPL::Variable::capillary_pressure, pos, t, dt);

        EXPECT_LT(S, S_prev);   // monotone decreasing
        EXPECT_LT(dS, 0.0);     // dS_Macro/dp_cap < 0 on both branches
        S_prev = S;

        double const eps = 1e-3 * p_cap;
        vars.capillary_pressure = p_cap - eps;
        double const S_minus = film.template value<double>(vars, pos, t, dt);
        vars.capillary_pressure = p_cap + eps;
        double const S_plus = film.template value<double>(vars, pos, t, dt);

        double const dS_fd = (S_plus - S_minus) / (2.0 * eps);
        double const d2S_fd = (S_plus - 2.0 * S + S_minus) / (eps * eps);

        // Same FD-comparison tolerances as the pure-Tuller test in this file
        // (the analytic film term is smooth; the residual is FD truncation).
        EXPECT_NEAR(dS, dS_fd, 1e-10);
        EXPECT_NEAR(d2S, d2S_fd, 1e-12);
    }
}

// Parsing the optional film keys must activate the branch: at fixed p_cap the
// film property reports a strictly larger saturation than the same base
// property with the film off (the film adds adsorbed water). Physics anchor:
// (d) qualitative sign of the added contribution. No expected magnitude.
TEST(MaterialPropertyLib, CreateSaturationTullerFilm)
{
    // A_Hamaker cited as above; a_v, phi_M synthetic fixtures (see note above).
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
        "  <macro_specific_surface>3.0e4</macro_specific_surface>"
        "  <hamaker_constant>2.2e-20</hamaker_constant>"
        "  <macro_porosity>0.3</macro_porosity>"
        "</property>";

    std::unique_ptr<MPL::Property> const film =
        Tests::createTestProperty(xml, MPL::createSaturationTuller);

    // Same base parameters, film off, as the reference.
    MPL::SaturationTuller const pure_obj{
        "saturation", 0.05, 1.0, 1.0, 0.858407, 1e-5, 0.072, 1.0};
    MPL::Property const& pure = pure_obj;

    MPL::VariableArray vars;
    ParameterLib::SpatialPosition const pos;
    double const t = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();

    // Saturated branch unchanged.
    vars.capillary_pressure = 0.0;
    EXPECT_DOUBLE_EQ(film->template value<double>(vars, pos, t, dt), 1.0);

    // In the unclamped window the film adds saturation.
    vars.capillary_pressure = 1e6;
    EXPECT_GT(film->template value<double>(vars, pos, t, dt),
              pure.template value<double>(vars, pos, t, dt));
}

