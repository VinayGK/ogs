// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause
//
// Strained-film disjoining law h(w_m, eps_v) — unit tests.
// Design: ProcessLib/RichardsMechanics/DSM/STRAINED_FILM_IMPLEMENTATION.md.
// Physics anchors (CLAUDE.md §3): derived identities (FD-vs-analytic chains,
// force-balance inversion residual), analytical limits (zero strain + zero
// load reduction), sign-only physical limits (load raises the potential).
// No fitted expected values; tolerances derive from the FD step / solver
// residual scales.
//
// Sample-state parameter values mirror the prior approved unit tests in
// DSMMicroMacroSingleIntegrationPoint.cpp (hamaker 6.0e-20 J, Sa 1000 m^2/kg,
// rho_SR 2650 kg/m^3, nS 0.6) — citation source: prior user-approved test
// code in this repository.

#include <gtest/gtest.h>

#include <cmath>
#include <limits>

#include "ProcessLib/RichardsMechanics/RichardsMechanicsFEM-impl.h"

using namespace ProcessLib::RichardsMechanics;

namespace
{
struct StrainedFilmSampleState
{
    double n_l = 0.3;
    double rho_lR = 1100.0;  // confined micro-liquid density scale, mirrors
                             // the memory note "micro EOS ~1100 kg/m^3"
    double active_nS = 0.6;
    double rho_SR = 2650.0;
    double hamaker = 6.0e-20;
    double Sa = 1000.0;
    // NegativeAttractive convention (the MS33 PRJ family): mu_lR < 0,
    // Pi = -rho*mu_lR > 0 (repulsive operational disjoining pressure).
    double sign = -1.0;
    double K_aug = 0.0;
    double lambda_aug = 0.0;
    double floor = 0.0;
};

double bareMu(StrainedFilmSampleState const& st, double const w)
{
    return computeVanDerWaalsMicroPotential(
               w, st.rho_lR, st.active_nS, st.rho_SR, st.hamaker, st.Sa,
               st.sign, st.K_aug, st.lambda_aug, 0.0, st.floor)
        .mu_lR;
}

double barePi(StrainedFilmSampleState const& st, double const w)
{
    return -st.rho_lR * bareMu(st, w);
}

StrainedFilmStateData state(StrainedFilmSampleState const& st,
                            FilmStrainCouplingMode const mode,
                            FilmStrainKappaMode const kappa,
                            double const eps_v, double const p_conf)
{
    return computeStrainedFilmState(mode, kappa, st.n_l, st.active_nS, eps_v,
                                    p_conf, st.rho_lR, st.rho_SR, st.hamaker,
                                    st.Sa, st.sign, st.K_aug, st.lambda_aug,
                                    st.floor, st.rho_lR);
}
}  // namespace

// Anchor: approved baseline — the defaults must stay Off/Aggregate so every
// existing PRJ is bit-for-bit unaffected.
TEST(RichardsMechanicsStrainedFilm, DefaultsAreOffAggregate)
{
    PotentialExchangeParameters params;
    EXPECT_EQ(params.film_strain_coupling, FilmStrainCouplingMode::Off);
    EXPECT_EQ(params.film_strain_kappa, FilmStrainKappaMode::Aggregate);
}

// Anchor: derived identity — FD-vs-analytic chains of the kinematic state.
TEST(RichardsMechanicsStrainedFilm, KinematicChainsFDConsistent)
{
    StrainedFilmSampleState st;
    double const eps_v = -0.02;  // compression
    double const p_conf = 1.0e6;
    double const d = 1e-7;

    for (auto const kappa :
         {FilmStrainKappaMode::Aggregate, FilmStrainKappaMode::Unity})
    {
        auto const s0 =
            state(st, FilmStrainCouplingMode::Kinematic, kappa, eps_v, p_conf);
        double const kappa_value =
            kappa == FilmStrainKappaMode::Aggregate ? st.active_nS : 1.0;
        EXPECT_NEAR(s0.w_eff, st.n_l * (1.0 + kappa_value * eps_v),
                    1e-14 * st.n_l);

        // d w_eff / d eps_v by central FD on the helper itself.
        auto const sp = state(st, FilmStrainCouplingMode::Kinematic, kappa,
                              eps_v + d, p_conf);
        auto const sm = state(st, FilmStrainCouplingMode::Kinematic, kappa,
                              eps_v - d, p_conf);
        double const fd_deps = (sp.w_eff - sm.w_eff) / (2.0 * d);
        EXPECT_NEAR(s0.dw_eff_deps_v, fd_deps, 1e-6 * std::abs(fd_deps));

        // d w_eff / d n_l by central FD in n_l.
        StrainedFilmSampleState stp = st;
        stp.n_l += d;
        StrainedFilmSampleState stm = st;
        stm.n_l -= d;
        double const fd_dnl =
            (state(stp, FilmStrainCouplingMode::Kinematic, kappa, eps_v,
                   p_conf)
                 .w_eff -
             state(stm, FilmStrainCouplingMode::Kinematic, kappa, eps_v,
                   p_conf)
                 .w_eff) /
            (2.0 * d);
        EXPECT_NEAR(s0.dw_eff_dnl, fd_dnl, 1e-6 * std::abs(fd_dnl));
    }
}

// Anchor: derived identity — on the loaded branch the inverted state must
// satisfy the film force balance Pi(w_eff) = p_conf to the solver residual
// scale (1e-9 relative; the Newton terminates at 1e-12 relative).
TEST(RichardsMechanicsStrainedFilm, EquilibriumInversionSolvesForceBalance)
{
    StrainedFilmSampleState st;
    double const Pi_unloaded = barePi(st, st.n_l);
    ASSERT_GT(Pi_unloaded, 0.0);

    // Loaded branch: p_conf above the unloaded disjoining pressure.
    double const p_loaded = 2.0 * Pi_unloaded;
    auto const s_loaded =
        state(st, FilmStrainCouplingMode::Equilibrium,
              FilmStrainKappaMode::Aggregate, 0.0, p_loaded);
    EXPECT_TRUE(s_loaded.loaded_branch);
    EXPECT_LT(s_loaded.w_eff, st.n_l);  // squeezed film
    EXPECT_NEAR(barePi(st, s_loaded.w_eff), p_loaded, 1e-9 * p_loaded);

    // Unloaded branch: p_conf below the branch point -> identity state.
    auto const s_unloaded =
        state(st, FilmStrainCouplingMode::Equilibrium,
              FilmStrainKappaMode::Aggregate, 0.0, 0.5 * Pi_unloaded);
    EXPECT_FALSE(s_unloaded.loaded_branch);
    EXPECT_DOUBLE_EQ(s_unloaded.w_eff, st.n_l);
    EXPECT_DOUBLE_EQ(s_unloaded.dw_eff_dnl, 1.0);
}

// Anchor: derived identity — the equilibrium inversion with the exponential
// augmentation active (exercises the Newton off the pure cubic seed).
TEST(RichardsMechanicsStrainedFilm, EquilibriumInversionWithAugmentation)
{
    StrainedFilmSampleState st;
    // Augmentation amplitude/decay from the prior approved dd1600 PRJ family
    // (potential_augmentation_prefactor 103879 J/kg, exponent 7.5e-7 m) —
    // citation source: Tests/Data/.../ANCHORS_MS33_ModelI/ms33_modelI_dd1600.prj.
    st.K_aug = 103879.0;
    st.lambda_aug = 7.5e-7;

    double const Pi_unloaded = barePi(st, st.n_l);
    ASSERT_GT(Pi_unloaded, 0.0);
    double const p_loaded = 3.0 * Pi_unloaded;
    auto const s =
        state(st, FilmStrainCouplingMode::Equilibrium,
              FilmStrainKappaMode::Aggregate, 0.0, p_loaded);
    ASSERT_TRUE(s.loaded_branch);
    EXPECT_NEAR(barePi(st, s.w_eff), p_loaded, 1e-9 * p_loaded);
}

// Anchor: analytical limit — at zero strain and zero confining pressure both
// strained modes reduce EXACTLY to the frozen-geometry evaluation point.
TEST(RichardsMechanicsStrainedFilm, ZeroStrainZeroLoadReducesToBareState)
{
    StrainedFilmSampleState st;
    for (auto const mode : {FilmStrainCouplingMode::Kinematic,
                            FilmStrainCouplingMode::Equilibrium})
    {
        auto const s =
            state(st, mode, FilmStrainKappaMode::Aggregate, 0.0, 0.0);
        EXPECT_DOUBLE_EQ(s.w_eff, st.n_l);
        EXPECT_DOUBLE_EQ(s.dw_eff_dnl, 1.0);
        EXPECT_DOUBLE_EQ(s.dw_eff_deps_v, mode ==
                             FilmStrainCouplingMode::Kinematic
                             ? st.n_l * st.active_nS
                             : 0.0);
        EXPECT_FALSE(s.loaded_branch);
    }
}

// Anchor: physical limit (sign only) — through the fold point, raising the
// confining pressure at fixed water content must RAISE mu_lR (the Derjaguin
// load term: squeezing confined liquid raises its chemical potential; the
// expulsion channel). No magnitude asserted.
TEST(RichardsMechanicsStrainedFilm, LoadRaisesPotentialAtFixedWaterContent)
{
    StrainedFilmSampleState st;
    PotentialExchangeParameters params;
    params.enabled = true;
    params.hamaker_constant = st.hamaker;
    params.specific_surface = st.Sa;
    params.micro_solid_density_reference = st.rho_SR;
    params.micro_solid_volume_fraction_reference = st.active_nS;
    params.micro_potential_convention =
        MicroPotentialConvention::NegativeAttractive;
    params.film_pressure_coupling = true;
    params.film_strain_coupling = FilmStrainCouplingMode::Kinematic;
    params.film_strain_kappa = FilmStrainKappaMode::Aggregate;

    auto const mu_at = [&](double const p_conf)
    {
        auto out = computeVanDerWaalsMicroPotential(
            st.n_l, st.rho_lR, st.active_nS, st.rho_SR, st.hamaker, st.Sa,
            st.sign, 0.0, 0.0, 0.0, 0.0);
        PotentialExchangeLocalSolveContext ctx;
        ctx.phi = 0.4;
        ctx.volumetric_strain = -0.01;  // compression, fixed
        ctx.volumetric_strain_prev = 0.0;
        ctx.confining_pressure_p_conf = p_conf;
        ctx.biot_coefficient = 1.0;
        applyFilmPressureMicroPotential(out, st.n_l, st.rho_lR, st.active_nS,
                                        ctx, params);
        return out.mu_lR;
    };

    double const mu_low = mu_at(1.0e5);
    double const mu_high = mu_at(1.0e6);
    EXPECT_GT(mu_high, mu_low);
}

// Anchor: approved baseline — with the strain coupling Off, the fold point
// must follow the existing (shipped) path: the integrable partner is active
// and the result differs from the bare law only by that partner; with the
// coupling ON the shipped partner must NOT also be applied (no double
// counting; replacement is exclusive). Verified structurally: Off and
// Kinematic disagree at finite strain (different mechanisms), while at zero
// strain and zero load Kinematic equals the bare law exactly but Off equals
// bare law + (zero) partner = bare law as well.
TEST(RichardsMechanicsStrainedFilm, ReplacementIsExclusiveAtZeroStrain)
{
    StrainedFilmSampleState st;
    PotentialExchangeParameters params;
    params.enabled = true;
    params.hamaker_constant = st.hamaker;
    params.specific_surface = st.Sa;
    params.micro_solid_density_reference = st.rho_SR;
    params.micro_solid_volume_fraction_reference = st.active_nS;
    params.micro_potential_convention =
        MicroPotentialConvention::NegativeAttractive;
    params.film_pressure_coupling = true;

    auto const mu_for = [&](FilmStrainCouplingMode const mode,
                            double const eps_v, double const p_conf)
    {
        auto out = computeVanDerWaalsMicroPotential(
            st.n_l, st.rho_lR, st.active_nS, st.rho_SR, st.hamaker, st.Sa,
            st.sign, 0.0, 0.0, 0.0, 0.0);
        PotentialExchangeLocalSolveContext ctx;
        ctx.phi = 0.4;
        ctx.volumetric_strain = eps_v;
        ctx.volumetric_strain_prev = 0.0;
        ctx.confining_pressure_p_conf = p_conf;
        ctx.biot_coefficient = 1.0;
        ctx.drained_bulk_modulus = 1.0e9;  // structural sample stiffness
        params.film_strain_coupling = mode;
        applyFilmPressureMicroPotential(out, st.n_l, st.rho_lR, st.active_nS,
                                        ctx, params);
        return out.mu_lR;
    };

    double const mu_bare = bareMu(st, st.n_l);

    // Zero strain, zero load: both reduce to the bare law (the Off path's
    // partner vanishes at eps_v = 0 and p_conf = 0; the strained path's
    // evaluation point and load term are identities there).
    EXPECT_DOUBLE_EQ(mu_for(FilmStrainCouplingMode::Off, 0.0, 0.0), mu_bare);
    EXPECT_DOUBLE_EQ(mu_for(FilmStrainCouplingMode::Kinematic, 0.0, 0.0),
                     mu_bare);

    // Finite strain + load: the two mechanisms must DIFFER (if the shipped
    // partner were still added on top of the strained law, the strained value
    // would carry both and this distinction would collapse).
    double const eps_v = -0.02;
    double const p_conf = 1.0e6;
    EXPECT_NE(mu_for(FilmStrainCouplingMode::Off, eps_v, p_conf),
              mu_for(FilmStrainCouplingMode::Kinematic, eps_v, p_conf));
}

// ── Live K(rho_d) helper (K_OF_RHO_D_LIVE.md; Vinay 2026-06-10) ────────────
// Physics anchor (CLAUDE.md §3a): analytical limit / derived identity — the
// helper must reproduce the table exactly at its knots (both schemes are
// node-preserving) and hold the endpoint values outside the range; off-mode
// must return the parse-time scalar bit-for-bit. Since 2026-08-26 (Vinay's
// interpolation-scheme decision) the LIVE path interpolates log-linearly
// (getValueLogLinear) while the parse-time frozen-K path keeps the linear
// getValue(); the tests below anchor BOTH schemes accordingly. The knot
// values below are STRUCTURAL in-test constants (not physical material
// parameters); expected values are derived in-file from the respective
// interpolation identity.
namespace
{
PotentialExchangeParameters liveKSampleParams()
{
    PotentialExchangeParameters params;
    params.micro_solid_density_reference = 2650.0;  // rho_SR, mirrors the
        // prior approved sample state above (kg/m^3).
    params.potential_augmentation_prefactor = 7.0;  // structural scalar K
    // Structural knots: K(1000) = 10, K(2000) = 30 (J/kg vs kg/m^3).
    params.potential_augmentation_prefactor_vs_dry_density =
        std::make_shared<AugmentationPrefactorTable const>(
            std::vector<double>{1000.0, 2000.0},
            std::vector<double>{10.0, 30.0});
    return params;
}

// phi such that rho_d = rho_SR*(1-phi) equals the requested dry density.
double phiForDryDensity(PotentialExchangeParameters const& params,
                        double const rho_d)
{
    return 1.0 - rho_d / params.micro_solid_density_reference;
}
}  // namespace

TEST(RichardsMechanicsLiveKOfRhoD, OffModeReturnsScalar)
{
    auto params = liveKSampleParams();
    params.potential_augmentation_prefactor_live_dry_density = false;
    // Off mode: scalar, regardless of table presence and finite phi.
    EXPECT_DOUBLE_EQ(7.0, effectiveAugmentationPrefactor(
                              params, phiForDryDensity(params, 1500.0)));

    // No table: scalar even when the mode flag is on.
    auto params_no_table = liveKSampleParams();
    params_no_table.potential_augmentation_prefactor_live_dry_density = true;
    params_no_table.potential_augmentation_prefactor_vs_dry_density = nullptr;
    EXPECT_DOUBLE_EQ(7.0,
                     effectiveAugmentationPrefactor(
                         params_no_table, phiForDryDensity(params, 1500.0)));
}

// 2026-08-26 RESTRUCTURE (Vinay's decision, 2026-08-26: "keep linear
// tests, add log-linear tests"): the LIVE free-function path switched to
// log-linear interpolation (getValueLogLinear). This test — formerly
// LiveModeEvaluatesTable, which asserted the same linear literals on
// effectiveAugmentationPrefactor() — was RETARGETED to the table's own
// LINEAR getValue(), preserving its linear expected literals unchanged.
// It anchors the still-present K-linear scheme, which the parse-time
// frozen-K path (CreateRichardsMechanicsProcess.cpp, ->getValue(
// *dry_density)) still uses. The live path is now covered by
// LiveModeEvaluatesTableLogLinear below. Physics anchor (CLAUDE.md §3a):
// linear-interpolation identity, derived in-file.
TEST(RichardsMechanicsLiveKOfRhoD, TableLinearGetValueAnchorsFrozenKParsePath)
{
    auto const params = liveKSampleParams();
    auto const& table = *params.potential_augmentation_prefactor_vs_dry_density;

    // At the knots: exact knot values.
    EXPECT_DOUBLE_EQ(10.0, table.getValue(1000.0));
    EXPECT_DOUBLE_EQ(30.0, table.getValue(2000.0));

    // Interior points: linear-interpolation identity
    // K(rho_d) = 10 + 20*(rho_d - 1000)/1000, derived in-file.
    double const rho_d_1 = 1500.0;
    double const rho_d_2 = 1750.0;
    double const expected_1 = 10.0 + 20.0 * (rho_d_1 - 1000.0) / 1000.0;
    double const expected_2 = 10.0 + 20.0 * (rho_d_2 - 1000.0) / 1000.0;
    EXPECT_DOUBLE_EQ(expected_1, table.getValue(rho_d_1));
    EXPECT_DOUBLE_EQ(expected_2, table.getValue(rho_d_2));
}

// 2026-08-26 NEW (companion to the retarget above; Vinay's decision,
// 2026-08-26): the LIVE path now evaluates the table log-linearly — ln(K)
// linear in rho_d between knots (getValueLogLinear). Physics anchor
// (CLAUDE.md §3a, analytical limit): expected values derived in-file from
// the log-linear identity K(x) = K_l * exp(t * ln(K_r/K_l)),
// t = (x - 1000)/1000, on the structural knots K(1000)=10, K(2000)=30;
// approved Vinay 2026-08-26.
TEST(RichardsMechanicsLiveKOfRhoD, LiveModeEvaluatesTableLogLinear)
{
    auto params = liveKSampleParams();
    params.potential_augmentation_prefactor_live_dry_density = true;

    // At the knots: node-preserving — exact knot values (the boundary-knot
    // clamp branches are byte-identical to getValue's).
    EXPECT_DOUBLE_EQ(10.0, effectiveAugmentationPrefactor(
                               params, phiForDryDensity(params, 1000.0)));
    EXPECT_DOUBLE_EQ(30.0, effectiveAugmentationPrefactor(
                               params, phiForDryDensity(params, 2000.0)));

    // Interior points, derived in-file (approved Vinay 2026-08-26):
    //   K(1500) = 10*exp(0.5*ln(30/10)) = 10*sqrt(3)
    //           = 17.320508075688775
    //   K(1750) = 10*exp(0.75*ln(3))   = 10*3^0.75
    //           = 22.795070569547775
    // Tolerance derived in-file: the phi round-trip rho_SR*(1 - phi) is
    // not ulp-exact (perturbs rho_d by O(rho_SR*eps) ~ 5e-13 kg/m^3;
    // |dK/drho_d| ~ 0.02 -> |dK| <~ 1e-14), so assert to 1e-12 relative,
    // two orders above that round-trip noise.
    double const expected_1500 = 17.320508075688775;
    double const expected_1750 = 22.795070569547775;
    EXPECT_NEAR(expected_1500,
                effectiveAugmentationPrefactor(
                    params, phiForDryDensity(params, 1500.0)),
                1e-12 * expected_1500);
    EXPECT_NEAR(expected_1750,
                effectiveAugmentationPrefactor(
                    params, phiForDryDensity(params, 1750.0)),
                1e-12 * expected_1750);

    // Non-finite phi (the context sentinel): fall back to the scalar —
    // dispatch is untouched by the interpolation scheme (moved here from
    // the retargeted linear test, 2026-08-26).
    EXPECT_DOUBLE_EQ(7.0,
                     effectiveAugmentationPrefactor(
                         params, std::numeric_limits<double>::infinity()));
    EXPECT_DOUBLE_EQ(7.0,
                     effectiveAugmentationPrefactor(
                         params, std::numeric_limits<double>::quiet_NaN()));
}

TEST(RichardsMechanicsLiveKOfRhoD, ClampsAtTableRangeEnds)
{
    auto params = liveKSampleParams();
    params.potential_augmentation_prefactor_live_dry_density = true;

    // Below rho_d_min and above rho_d_max the endpoint values are held
    // (PiecewiseLinearInterpolation::getValue endpoint hold).
    EXPECT_DOUBLE_EQ(10.0, effectiveAugmentationPrefactor(
                               params, phiForDryDensity(params, 500.0)));
    EXPECT_DOUBLE_EQ(30.0, effectiveAugmentationPrefactor(
                               params, phiForDryDensity(params, 2600.0)));
}

// ── Live K(rho_d) analytic tangent (Vinay 2026-06-12 Jacobian completion) ──
// Physics anchor (per Vinay's task spec, 2026-06-12): FD-vs-analytic
// agreement (derived identity) of the live-K tangent against a central
// finite difference of the VALUE actually used in the residual. Inside a
// table segment the value is exactly linear in phi, so the central FD is
// exact to roundoff; tolerances are derived in-file from the tangent scale.
// Knot values are STRUCTURAL in-test constants (not physical parameters),
// mirroring the approved live-K tests above.
// 2026-08-26 RESTRUCTURE (Vinay's decision, 2026-08-26: "keep linear
// tests, add log-linear tests"): the LIVE tangent
// effectiveAugmentationPrefactorPhiDerivative() switched to the log-linear
// companion slope (getSegmentSlopeLogLinear). This test — formerly
// AnalyticPhiTangentMatchesFDInsideSegment, which asserted the linear
// segment-slope chain on the free function — was RETARGETED to the table's
// own LINEAR getSegmentSlope(), preserving its linear expected literals
// unchanged. It anchors the K-linear scheme still used by the parse-time
// frozen-K path (CreateRichardsMechanicsProcess.cpp). The live-path
// FD-vs-analytic, zero-case, and mu-level Jacobian-chain checks moved to
// AnalyticPhiTangentMatchesFDLogLinear below. Physics anchor (CLAUDE.md
// §3a): linear segment-slope identity, derived in-file.
TEST(RichardsMechanicsLiveKOfRhoD, TableLinearSegmentSlopeAnchorsFrozenKParsePath)
{
    auto const params = liveKSampleParams();
    auto const& table = *params.potential_augmentation_prefactor_vs_dry_density;

    // Derived in-file: segment slope dK/drho_d = (30-10)/(2000-1000)
    // = 0.02 (J/kg)/(kg/m^3); chain dK/dphi = -rho_SR * slope (= -53, the
    // pre-2026-08-26 live-path expectation, preserved here as the
    // LINEAR-scheme identity).
    double const slope = (30.0 - 10.0) / (2000.0 - 1000.0);
    EXPECT_NEAR(slope, table.getSegmentSlope(1500.0), 1e-12 * slope);
    double const expected_dK_dphi =
        -params.micro_solid_density_reference * slope;
    EXPECT_NEAR(expected_dK_dphi,
                -params.micro_solid_density_reference *
                    table.getSegmentSlope(1500.0),
                1e-12 * std::abs(expected_dK_dphi));
}

// 2026-08-26 NEW (successor of the retargeted linear tangent test above,
// for the LIVE path; Vinay's decision, 2026-08-26): analytic-vs-FD
// self-consistency of effectiveAugmentationPrefactorPhiDerivative() under
// the log-linear tangent. Physics anchor (CLAUDE.md §3a, analytical
// derivation in-file; approved Vinay 2026-08-26): on the structural knots
// K(1000)=10, K(2000)=30,
//   dK/drho_d = K(rho_d) * ln(30/10)/(2000-1000) = K(rho_d)*ln(3)/1000
//   dK/dphi   = -rho_SR * dK/drho_d
//             = -2650 * 17.320508075688775 * ln(3)/1000
//             = -50.425585997506346 at rho_d = 1500.
TEST(RichardsMechanicsLiveKOfRhoD, AnalyticPhiTangentMatchesFDLogLinear)
{
    auto params = liveKSampleParams();
    params.potential_augmentation_prefactor_live_dry_density = true;

    double const phi_mid = phiForDryDensity(params, 1500.0);  // mid-segment
    double const analytic =
        effectiveAugmentationPrefactorPhiDerivative(params, phi_mid);
    // Closed form derived in the comment above (approved Vinay
    // 2026-08-26); tolerance 1e-12 relative, derived in-file from the phi
    // round-trip noise exactly as in LiveModeEvaluatesTableLogLinear.
    double const expected_dK_dphi = -50.425585997506346;
    EXPECT_NEAR(expected_dK_dphi, analytic,
                1e-12 * std::abs(expected_dK_dphi));

    // Central FD of the residual-side value (stays within the segment).
    double const d_phi = 1e-6;  // step on the O(1) phi scale
    double const fd = (effectiveAugmentationPrefactor(params, phi_mid + d_phi) -
                       effectiveAugmentationPrefactor(params, phi_mid - d_phi)) /
                      (2.0 * d_phi);
    // The log-linear value is smooth (exponential) in phi inside a segment;
    // central-FD relative truncation error ~ (ln(3)/1000 * rho_SR*d_phi)^2/6
    // ~ 1.4e-12 (derived in-file) — far inside the 1e-9 tolerance.
    EXPECT_NEAR(fd, analytic, 1e-9 * std::abs(analytic));

    // Off mode / sentinel phi / no table: tangent identically zero
    // (residual uses the parse-time scalar there).
    auto params_off = liveKSampleParams();
    params_off.potential_augmentation_prefactor_live_dry_density = false;
    EXPECT_DOUBLE_EQ(
        0.0, effectiveAugmentationPrefactorPhiDerivative(params_off, phi_mid));
    EXPECT_DOUBLE_EQ(0.0,
                     effectiveAugmentationPrefactorPhiDerivative(
                         params, std::numeric_limits<double>::quiet_NaN()));

    // The mu-level aug K-partials feeding the Jacobian chain: mu_aug is
    // LINEAR in K, so central FD in K of the vdW helper is exact to
    // roundoff. Sample state mirrors the approved StrainedFilmSampleState
    // (file header citation); lambda chosen so xi = h/lambda is O(1)
    // (structural, not physical).
    StrainedFilmSampleState st;
    st.K_aug = 10.0;     // structural K (J/kg), matches the table knot
    st.lambda_aug = 2e-7;  // m, structural; h = n_l/(nS*rho_SR*Sa) ~ 1.9e-7 m
    auto const vdw = computeVanDerWaalsMicroPotential(
        st.n_l, st.rho_lR, st.active_nS, st.rho_SR, st.hamaker, st.Sa, st.sign,
        st.K_aug, st.lambda_aug, 0.0, st.floor);
    double const dK = 1e-3 * st.K_aug;  // step derived from K scale
    auto const mu_at_K = [&](double const K)
    {
        return computeVanDerWaalsMicroPotential(
            st.n_l, st.rho_lR, st.active_nS, st.rho_SR, st.hamaker, st.Sa,
            st.sign, K, st.lambda_aug, 0.0, st.floor);
    };
    double const fd_dmu_dK =
        (mu_at_K(st.K_aug + dK).mu_lR - mu_at_K(st.K_aug - dK).mu_lR) /
        (2.0 * dK);
    EXPECT_NEAR(fd_dmu_dK, vdw.dmu_lR_dK, 1e-9 * std::abs(vdw.dmu_lR_dK));
    double const fd_ddmudnl_dK = (mu_at_K(st.K_aug + dK).dmu_lR_dnl -
                                  mu_at_K(st.K_aug - dK).dmu_lR_dnl) /
                                 (2.0 * dK);
    EXPECT_NEAR(fd_ddmudnl_dK, vdw.ddmu_lR_dnl_dK,
                1e-9 * std::abs(vdw.ddmu_lR_dnl_dK));
}

TEST(RichardsMechanicsLiveKOfRhoD, AnalyticPhiTangentClampedEdgesAndKnots)
{
    // Three structural knots so an INTERIOR knot exists:
    // K(1000)=10, K(1500)=16, K(2000)=30 (J/kg vs kg/m^3). The K-LINEAR
    // chord slopes 0.012 and 0.028 (derived in-file) belong to
    // getValue()/getSegmentSlope(), i.e. to the parse-time frozen-K path.
    // The LIVE path (effectiveAugmentationPrefactor and its Phi
    // derivative) has used the LOG-LINEAR pair since 2026-08-26; its
    // tangents are derived in-file below from the same knots.
    PotentialExchangeParameters params;
    params.micro_solid_density_reference = 2650.0;  // mirrors sample state
    params.potential_augmentation_prefactor_live_dry_density = true;
    params.potential_augmentation_prefactor_vs_dry_density =
        std::make_shared<AugmentationPrefactorTable const>(
            std::vector<double>{1000.0, 1500.0, 2000.0},
            std::vector<double>{10.0, 16.0, 30.0});
    auto const phi_of = [&](double const rho_d)
    { return 1.0 - rho_d / params.micro_solid_density_reference; };

    // Outside the range: the clamped evaluation is flat -> tangent 0.
    // (Strictly-outside points only here: the phi round-trip
    // rho_SR*(1 - phi_of(rho_d)) is not exact to the last ulp, so the
    // AT-knot convention is tested on the table directly below.)
    for (double const rho_d : {500.0, 2600.0})
    {
        EXPECT_DOUBLE_EQ(0.0, effectiveAugmentationPrefactorPhiDerivative(
                                  params, phi_of(rho_d)));
    }
    // AT the edge knots (exact arguments): slope 0, the documented
    // one-sided/zero convention of AugmentationPrefactorTable, mirroring
    // getValue's <=/>= clamp branches. Asserted for BOTH tangent methods:
    // getSegmentSlope() pins the K-linear frozen-K parse path (kept, it is
    // that path's coverage), getSegmentSlopeLogLinear() the method the
    // LIVE path actually calls since the 2026-08-26 scheme switch. Both
    // boundary knots and both strictly-outside points are asserted for the
    // log-linear tangent, so a regression that narrowed the <=/>= clamp to
    // </> (which would index below the first segment at rho_d = 1000)
    // fails here.
    auto const& table = *params.potential_augmentation_prefactor_vs_dry_density;
    EXPECT_DOUBLE_EQ(0.0, table.getSegmentSlope(1000.0));
    EXPECT_DOUBLE_EQ(0.0, table.getSegmentSlope(2000.0));
    EXPECT_DOUBLE_EQ(0.0, table.getSegmentSlopeLogLinear(1000.0));
    EXPECT_DOUBLE_EQ(0.0, table.getSegmentSlopeLogLinear(2000.0));
    EXPECT_DOUBLE_EQ(0.0, table.getSegmentSlopeLogLinear(500.0));
    EXPECT_DOUBLE_EQ(0.0, table.getSegmentSlopeLogLinear(2600.0));
    double const d_phi = 1e-6;
    double const phi_out = phi_of(500.0);  // fully outside, FD stays outside
    EXPECT_DOUBLE_EQ(
        0.0, (effectiveAugmentationPrefactor(params, phi_out + d_phi) -
              effectiveAugmentationPrefactor(params, phi_out - d_phi)) /
                 (2.0 * d_phi));

    // Interior knot rho_d = 1500 (exact argument): LEFT-segment slope
    // (one-sided), consistent with getValue's lower_bound interval
    // selection. Derived in-file: (16-10)/500 = 0.012 (J/kg)/(kg/m^3).
    double const expected_left_slope = (16.0 - 10.0) / 500.0;
    EXPECT_NEAR(expected_left_slope, table.getSegmentSlope(1500.0),
                1e-12 * expected_left_slope);

    // Same interior knot, LOG-LINEAR pair — the methods the LIVE path
    // actually calls since 2026-08-26. Physics anchor (CLAUDE.md §3a):
    // analytical limit / derived identity; every expected number below is
    // computed in-file from the same structural knots constructed above.
    //
    // (1) Node preservation: getValueLogLinear returns the stored knot
    //     values. On THIS knot set that is the whole of what the three
    //     assertions below establish — the boundary knots 1000 and 2000
    //     come back through the endpoint-hold clamp, and at the interior
    //     knot 1500 the round-trip 10*exp(1*ln(16/10)) happens to land
    //     exactly on 16 (measured 2026-08-31 by evaluating the same
    //     expression standalone under both toolchains on this machine,
    //     Apple clang 21.0.0 — which builds this test — and Homebrew
    //     clang 22.1.8, arm64, -O0 and -O2 each), so they hold with OR
    //     without the t == 1 branch of
    //     logLinearValueOnSegment(). They therefore do not pin that
    //     branch; the knot pair whose round-trip does miss, and which
    //     does pin it, is asserted in the test
    //     InteriorKnotBitExactWhereRoundTripMisses below.
    double const K_interior_knot = 16.0;  // the interior knot value above
    EXPECT_EQ(K_interior_knot, table.getValueLogLinear(1500.0));
    EXPECT_EQ(10.0, table.getValueLogLinear(1000.0));
    EXPECT_EQ(30.0, table.getValueLogLinear(2000.0));

    // (2) One-sided LEFT-segment convention of the log-linear tangent
    //     (idx = lower_bound - 1), derived in-file from the same knots:
    //       left  [1000,1500]: dK/drho_d = K(1500)*ln(16/10)/(1500-1000)
    //       right [1500,2000]: dK/drho_d = K(1500)*ln(30/16)/(2000-1500)
    //     The convention must select the LEFT value. A lower_bound ->
    //     upper_bound regression would change the answer ONLY at an exact
    //     knot — the FD checks below (1250/1750) live strictly inside one
    //     segment and are blind to it — so the right-segment candidate is
    //     computed as well and asserted to lie far outside the tolerance
    //     used, which is what makes the assertion discriminating.
    //     Tolerance: the same 1e-12 relative used by the linear
    //     interior-knot assertion above (the table is called directly here,
    //     so there is no phi round-trip; 1e-12 is conservative).
    double const expected_left_loglinear_slope =
        K_interior_knot * std::log(16.0 / 10.0) / (1500.0 - 1000.0);
    double const right_loglinear_slope =
        K_interior_knot * std::log(30.0 / 16.0) / (2000.0 - 1500.0);
    double const loglinear_slope_tol = 1e-12 * expected_left_loglinear_slope;
    EXPECT_NEAR(expected_left_loglinear_slope,
                table.getSegmentSlopeLogLinear(1500.0), loglinear_slope_tol);
    EXPECT_GT(std::abs(right_loglinear_slope - expected_left_loglinear_slope),
              loglinear_slope_tol);

    // Interior of each segment: FD-vs-analytic under the LOG-LINEAR live
    // value. The value is exponential (no longer linear) in phi inside a
    // segment, so the central FD is no longer exact; its relative
    // truncation error is ~(ln(K_r/K_l)/(x_r-x_l) * rho_SR*d_phi)^2/6
    // ~ 1e-12 here (derived in-file, the same estimate as in
    // AnalyticPhiTangentMatchesFDLogLinear), far inside the unchanged
    // 1e-9 relative tolerance below.
    for (double const rho_d : {1250.0, 1750.0})
    {
        double const phi_c = phi_of(rho_d);
        double const analytic =
            effectiveAugmentationPrefactorPhiDerivative(params, phi_c);
        double const fd =
            (effectiveAugmentationPrefactor(params, phi_c + d_phi) -
             effectiveAugmentationPrefactor(params, phi_c - d_phi)) /
            (2.0 * d_phi);
        EXPECT_NEAR(fd, analytic, 1e-9 * std::abs(analytic));
    }
}

// ── Interior-knot bit-exactness, on a knot pair whose log/exp round-trip
// actually misses ─────────────────────────────────────────────────────────
// Physics anchor (CLAUDE.md §3a): analytical limit / interpolation identity
// — node preservation, K(x_i) = K_i exactly at every knot of the table,
// which is what the t == 1 (and t == 0) branch of
// AugmentationPrefactorTable::logLinearValueOnSegment() exists to make true
// bit-for-bit rather than to ~1 ULP.
//
// The knots are STRUCTURAL: this pair is used because it exercises the
// log/exp round-trip at t == 1, NOT because anything about the physics of
// dd900 is being asserted. Nothing below depends on K(900) being a
// calibrated value — only on the pair being one whose round-trip misses.
// Sources of the three numbers (CLAUDE.md §1.1 item 3, prior commits
// traceable to user-approved work):
//   K(900)  = 4367.227700212952 J/kg — the SUPERSEDED dd900 calibration
//             record, kept on record in the provenance block of
//             Tests/Data/RichardsMechanics/ANCHORS_MS33_ModelI/
//             ms33_modelI_dd900.prj (lines 107/123/412) and in commit
//             642a8f867a. The value that superseded it is deliberately not
//             used here: with the current 900-knot the segment round-trips
//             exactly and the test would stop discriminating (see below).
//   K(1400) = 46000.0, K(1600) = 104689.9129 J/kg — the shipped live table
//             <prefactors> at 7ec39ecf4c (e.g.
//             ms33_modelIII_gapswitch.prj:192-193).
//
// Why this knot set and not the 10/16/30 one used by
// AnalyticPhiTangentClampedEdgesAndKnots above: there
// 10*exp(1*ln(16/10)) == 16 and 16*exp(1*ln(30/16)) == 30 EXACTLY, and all
// four knots of the table shipped at 7ec39ecf4c round-trip exactly too
// (both measured 2026-08-31 by evaluating the same expressions standalone,
// Homebrew clang 22.1.8 arm64, -O0 and -O2). Bit-exactness assertions on
// those knots pass with or without the t == 1 branch, so they cannot pin
// it — the defect the branch removes is latent on the shipped table. On
// the 900 -> 1400 pair below the round-trip K_l*exp(1*ln(K_r/K_l)) returns
// 45999.999999999993 instead of 46000.0 (rel -1.58e-16; the same miss the
// landing audit recorded in ProcessLib/RichardsMechanics/DSM/AGENTS.md as
// "1400 low by 1 ULP, rel 1.6e-16"). Measured 2026-08-31, before this test
// was added to the suite, by running getValueLogLinear/locateSegment/
// logLinearValueOnSegment extracted verbatim into a standalone program and
// compiled twice, with and without the t == 1 branch: both assertions below
// fail without it and pass with it, while the 10/16/30 assertions and the
// shipped table's four knots pass either way. That standalone verdict was
// re-confirmed 2026-08-31 under BOTH toolchains present on this machine, at
// -O0 and -O2 each: Apple clang 21.0.0 (arm64-apple-darwin25.6.0,
// /usr/bin/c++ — the compiler that actually builds this test) and Homebrew
// clang 22.1.8. MEASURED through the built test binary 2026-08-31 18:06,
// Release, Apple clang 21.0.0: this test PASSES in testrunner
// (--gtest_filter='RichardsMechanics*', 46 tests, 44 passed, 2 pre-existing
// skips, 0 failures). The FAILING half — the pre-fix accessor going red —
// is measurable only in the standalone harness, since reproducing it
// in-binary would mean removing the branch from the shipped header; it is
// harness-measured, not a testrunner observation.
//
// The second assertion pins the BRANCH rather than a number: the same knot
// value reached through the endpoint-hold clamp of a two-knot table — a
// path that returns the stored value verbatim, with no logarithm in it —
// must be bit-identical to what the interior-knot path returns, and its
// expected side is produced by the table itself, so no expected literal is
// written down at all. The other candidate form — computing
// K_l*exp(1*ln(K_r/K_l)) inline and asserting that it DIFFERS from the
// accessor — discriminates equally well on this toolchain (in the same
// measurement the inline round-trip reproduced the pre-fix accessor value
// bit for bit, at -O0 and -O2), but is not used: once the branch is in
// place that comparison no longer probes the accessor at all, only whether
// the platform's libm round-trips exactly, so a differently-rounding libm
// would turn it red on correct code. The two assertions below instead
// degrade to tautologies on such a platform rather than failing.
//
// The table is called directly rather than through
// effectiveAugmentationPrefactor() because the phi round-trip
// rho_SR*(1 - phi) is not ulp-exact (same reason as in the AT-knot
// assertions of AnalyticPhiTangentClampedEdgesAndKnots above), which would
// defeat a bit-exact comparison.
TEST(RichardsMechanicsLiveKOfRhoD, InteriorKnotBitExactWhereRoundTripMisses)
{
    double const K_900 = 4367.227700212952;  // J/kg, superseded dd900 knot
    double const K_1400 = 46000.0;           // J/kg
    double const K_1600 = 104689.9129;       // J/kg
    // The 1600 knot is present only so that 1400 is an INTERIOR knot: in a
    // two-knot table 1400 is the top support point and is reached through
    // the endpoint-hold clamp, which never enters the log/exp path.
    AugmentationPrefactorTable const table(
        std::vector<double>{900.0, 1400.0, 1600.0},
        std::vector<double>{K_900, K_1400, K_1600});

    // Interior knot: lower_bound places rho_d = 1400 in the LEFT segment
    // [900, 1400] at t == 1, where this pair's round-trip lands 1 ULP below
    // K_r. Node preservation is bit-exact only if the stored knot value is
    // returned.
    EXPECT_EQ(K_1400, table.getValueLogLinear(1400.0));

    // Same knot, reached instead through the endpoint-hold clamp of a
    // two-knot table ending at it: that branch returns the stored value
    // verbatim. The interior-knot path must agree with it bit for bit.
    AugmentationPrefactorTable const table_clamped(
        std::vector<double>{900.0, 1400.0}, std::vector<double>{K_900, K_1400});
    EXPECT_EQ(table_clamped.getValueLogLinear(1400.0),
              table.getValueLogLinear(1400.0));

    // Boundary knots of the three-knot table, for completeness: unchanged
    // clamp convention, stored values held exactly.
    EXPECT_EQ(K_900, table.getValueLogLinear(900.0));
    EXPECT_EQ(K_1600, table.getValueLogLinear(1600.0));
}

// ── §8 NEW TEST (review 2026-06-14): assembled displacement-channel
// Jacobian consistency for exact + kinematic + live-K at a finite-eps_v
// compliant state ──────────────────────────────────────────────────────────
// Physics anchor (CLAUDE.md §3d): symmetry / derived identity — the analytic
// displacement (eps_v) tangents the assembly inserts into the Jacobian must
// equal the central finite difference of the residual quantities they
// linearize. NO Vinay expected value; tolerances derive from the FD step and
// the live-K table's piecewise-linear (C0) structure.
//
// Covers the displacement-channel tangent FORMULAE the fixed assembly inserts:
//   (A) H2/M1 — the exact-route mu_lR eps_v tangent
//       g_cut * pair.dmu_mech_deps_v (the dispatch M1 wires into K[p,u]). The
//       analytic uses the SAME effective K the bare mu_lR is built with (H2):
//       FD of the assembled exact mu_lR(eps_v) must match it. With the pre-H2
//       scalar-K passed into the pair, g_cut = bare_live / pair_bare_scalar
//       diverges from 1 under live K, breaking this identity.
//   (B) M2 — the live-K swelling-eigenstress eps_v tangent
//       d(delta_sigma_sw)/dK * dK/dphi * dphi/deps_v. This is the term M2 adds
//       to K[u,u]; pre-M2 the analytic side of this identity (the chain) was
//       absent from the Jacobian entirely. The leg FDs the residual increment
//       (computeSwellingStressIncrement) through the live-K phi channel and
//       checks the analytic chain reproduces it.
//
// SCOPE: this is a HELPER-LEVEL FD-vs-analytic identity on the tangent
// formulae (the assembly reconstructs these same expressions); it does NOT
// drive the global assembleWithJacobian (no run-level FEM harness exists in
// this unit-test directory). The analytic legs are reconstructed from the SAME
// public helpers the FEM assembly calls (effectiveAugmentationPrefactor
// [PhiDerivative], computeStrainedFilmEnergyPair, computeSwellingStress-
// Increment), so it is a genuine FD-vs-analytic check, not a self-comparison.
// (A) FAILS on pre-H2 code in any state where the macro-floor cutoff is active
// under live K; (B)'s analytic chain did not exist pre-M2.
TEST(RichardsMechanicsLiveKOfRhoD, AssembledDisplacementTangentExactKinematicLiveK)
{
    using KV = MathLib::KelvinVector::KelvinVectorType<2>;
    auto const& I2 =
        MathLib::KelvinVector::Invariants<MathLib::KelvinVector::
            kelvin_vector_dimensions(2)>::identity2;

    // Sample state: exact + kinematic + live K, finite (compressive) eps_v with
    // the dry density rho_d = rho_SR*(1-phi) in the table interior so dK/dphi is
    // a genuine (nonzero) segment slope. Structural constants (CLAUDE.md §1.2);
    // material values mirror the prior approved tests in this file.
    PotentialExchangeParameters params;
    params.enabled = true;
    params.film_pressure_coupling = true;
    params.film_strain_coupling = FilmStrainCouplingMode::Kinematic;
    params.film_energy_route = FilmEnergyRoute::Exact;
    params.film_strain_kappa = FilmStrainKappaMode::Aggregate;
    params.micro_potential_convention =
        MicroPotentialConvention::NegativeAttractive;
    params.hamaker_constant = 6.0e-20;
    params.specific_surface = 1000.0;
    params.micro_solid_density_reference = 2650.0;  // rho_SR
    // lambda = characteristic film thickness [m]; structural h-scale
    // h0 = n_l/(nS*rho_SR*Sa) ~ 1.9e-7 m (xi0 ~ 1), as in this file's other
    // augmentation tests. Must be > 0 when K > 0 (law guard).
    params.potential_augmentation_exponent =
        0.30 / (0.70 * 2650.0 * 1000.0);  // m
    params.potential_augmentation_prefactor = 20.0;  // structural scalar fallback
    // Live K table: K(rho_d) over a span around the sample rho_d, finite slope.
    params.potential_augmentation_prefactor_live_dry_density = true;
    params.potential_augmentation_prefactor_vs_dry_density =
        std::make_shared<AugmentationPrefactorTable const>(
            std::vector<double>{1000.0, 2000.0},
            std::vector<double>{10.0, 50.0});  // J/kg vs kg/m^3

    double const sign = microPotentialSignFactorFromParameters(params);
    double const rho_lR = 1100.0;     // micro liquid density scale
    double const rho_LR = 1000.0;     // bulk
    double const n_l = 0.30;
    double const n_l_prev = 0.27;
    double const biot = 1.0;
    double const K_drained = 1.5e8;   // Pa, prior approved test value
    double const eps_v = -0.02;       // compression
    double const eps_v_prev = 0.0;
    double const phi0 = 0.40;         // rho_d = 2650*0.6 = 1590 kg/m^3 (interior)
    double const active_nS = 1.0 - n_l;
    double const kappa = active_nS;   // Aggregate

    // ── (A) H2/M1 — exact-route mu_lR eps_v tangent ───────────────────────────
    // Assembled exact mu_lR(eps_v) = bare(n_l; K) + g_cut * pair.mu_mech, with
    // g_cut = bare/pair.mu_bare_pre (== 1 here, cutoff inactive). Phi held fixed
    // to isolate the explicit eps_v channel (the live-K phi channel is part (B)).
    {
        double const K_aug = effectiveAugmentationPrefactor(params, phi0);
        auto const bare_mu = computeVanDerWaalsMicroPotential(
            n_l, rho_lR, active_nS, params.micro_solid_density_reference,
            params.hamaker_constant, params.specific_surface, sign, K_aug,
            params.potential_augmentation_exponent, 0.0,
            params.micro_water_content_floor);
        auto const mu_exact = [&](double const e)
        {
            auto const pr = computeStrainedFilmEnergyPair(
                n_l, e, kappa, biot, K_drained, true, rho_lR, active_nS,
                params.micro_solid_density_reference, params.hamaker_constant,
                params.specific_surface, sign, K_aug,
                params.potential_augmentation_exponent, 0.0,
                params.micro_water_content_floor);
            double const g_cut = bare_mu.mu_lR / pr.mu_bare_pre;  // [-]
            return bare_mu.mu_lR + g_cut * pr.mu_mech;            // J/kg
        };
        auto const pr0 = computeStrainedFilmEnergyPair(
            n_l, eps_v, kappa, biot, K_drained, true, rho_lR, active_nS,
            params.micro_solid_density_reference, params.hamaker_constant,
            params.specific_surface, sign, K_aug,
            params.potential_augmentation_exponent, 0.0,
            params.micro_water_content_floor);
        double const g_cut0 = bare_mu.mu_lR / pr0.mu_bare_pre;  // [-]
        double const analytic = g_cut0 * pr0.dmu_mech_deps_v;  // J/kg per strain
        double const h = 1e-7;
        double const fd = (mu_exact(eps_v + h) - mu_exact(eps_v - h)) / (2 * h);
        EXPECT_NEAR(fd, analytic, 5e-5 * std::abs(analytic) + 1e-10)
            << "exact-route mu_lR eps_v tangent (H2/M1)";
    }

    // ── (B) M2 — live-K swelling-eigenstress eps_v tangent (through phi) ───────
    // The residual delta_sigma_sw uses K = effectiveAugmentationPrefactor(phi);
    // with the live-K table phi(eps_v) couples sigma_sw to displacement. FD the
    // residual increment w.r.t. eps_v THROUGH the live-K channel only (phi moved
    // by the PorosityFromMassBalance chain dphi/deps_v on the SAME eps_v step,
    // n_l and the explicit-eps_v drained line held), and compare to the M2
    // analytic chain d(delta_sigma_sw)/dK * dK/dphi * dphi/deps_v.
    {
        // PorosityFromMassBalance dphi/deps_v = (alpha-phi)/(1+w); here on a pure
        // strain step w = delta_eps_v, alpha = biot. Interior -> no clamp.
        double const alpha = biot;
        double const w_phi = eps_v - eps_v_prev;  // pure strain step
        double const dphi_deps_v = (alpha - phi0) / (1.0 + w_phi);  // [-]
        double const dK_dphi =
            effectiveAugmentationPrefactorPhiDerivative(params, phi0);
        ASSERT_NE(dK_dphi, 0.0);  // sample state must be in the table interior

        auto const sigma_inc = [&](double const K_aug) -> double
        {
            // Build a params copy with the live-K table OFF and the scalar set
            // to K_aug, so computeSwellingStressIncrement uses exactly K_aug
            // (isolating the K dependence; eps_v/n_l/p_conf held).
            auto p = params;
            p.potential_augmentation_prefactor_live_dry_density = false;
            p.potential_augmentation_prefactor = K_aug;
            // p_conf NaN -> drain dropped; the K dependence enters through Pi.
            KV const inc = computeSwellingStressIncrement<2>(
                n_l_prev, n_l, active_nS, rho_lR, rho_lR, rho_LR,
                MathLib::KelvinVector::KelvinMatrixType<2>::Identity() *
                    K_drained,
                p, biot, std::numeric_limits<double>::quiet_NaN(), eps_v,
                eps_v_prev, std::numeric_limits<double>::quiet_NaN());
            return inc.dot(I2) / I2.dot(I2);  // scalar on identity2 [Pa]
        };
        double const K0 = effectiveAugmentationPrefactor(params, phi0);
        double const dK = 1e-4 * std::max(1.0, std::abs(K0));
        double const d_sigma_dK =
            (sigma_inc(K0 + dK) - sigma_inc(K0 - dK)) / (2 * dK);  // Pa per J/kg
        double const analytic = d_sigma_dK * dK_dphi * dphi_deps_v;  // Pa/strain

        // FD of the residual increment w.r.t. eps_v through phi(eps_v) only:
        // move K by K(phi0 + dphi_deps_v*he) on an eps_v step he.
        double const he = 1e-6;
        double const K_plus = effectiveAugmentationPrefactor(
            params, phi0 + dphi_deps_v * he);
        double const K_minus = effectiveAugmentationPrefactor(
            params, phi0 - dphi_deps_v * he);
        double const fd = (sigma_inc(K_plus) - sigma_inc(K_minus)) / (2 * he);
        EXPECT_NEAR(fd, analytic,
                    5e-4 * std::abs(analytic) + 1e-3 * std::abs(d_sigma_dK))
            << "live-K swelling-eigenstress eps_v tangent (M2)";
    }
}
