// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <cmath>
#include <limits>

#include "ProcessLib/RichardsMechanics/PotentialExchangeParameters.h"
#include "ProcessLib/RichardsMechanics/ConstitutiveRelations/MicroWaterContent.h"
#include "ProcessLib/RichardsMechanics/ConstitutiveRelations/PotentialExchange.h"

using namespace ProcessLib::RichardsMechanics;

namespace
{
double centralDiff(auto&& f, double const x, double const h)
{
    return (f(x + h) - f(x - h)) / (2.0 * h);
}
}  // namespace

TEST(RichardsMechanics, PotentialExchangeYoungLaplaceMacroPotential)
{
    double const rho = 1000.0;
    double const ptol = 1.0;

    {
        auto const mu = computeYoungLaplaceMacroPotential(0.0, rho, ptol);
        EXPECT_TRUE(mu.saturated_branch);
        EXPECT_DOUBLE_EQ(mu.mu_LR, 0.0);
        EXPECT_DOUBLE_EQ(mu.dmu_LR_dpLR, 0.0);
    }

    {
        double const pLR = -ptol;
        auto const mu = computeYoungLaplaceMacroPotential(pLR, rho, ptol);
        EXPECT_FALSE(mu.saturated_branch);
        EXPECT_DOUBLE_EQ(mu.mu_LR, pLR / rho);
        EXPECT_DOUBLE_EQ(mu.dmu_LR_dpLR, 1.0 / rho);
    }

    {
        double const pLR = -ptol + 1e-6;
        auto const mu = computeYoungLaplaceMacroPotential(pLR, rho, ptol);
        EXPECT_TRUE(mu.saturated_branch);
        EXPECT_DOUBLE_EQ(mu.mu_LR, 0.0);
        EXPECT_DOUBLE_EQ(mu.dmu_LR_dpLR, 0.0);
    }

    {
        double const pLR = -1.0e4;
        auto const mu = computeYoungLaplaceMacroPotential(pLR, rho, ptol);
        EXPECT_FALSE(mu.saturated_branch);
        EXPECT_DOUBLE_EQ(mu.mu_LR, pLR / rho);
        EXPECT_DOUBLE_EQ(mu.dmu_LR_dpLR, 1.0 / rho);
        EXPECT_DOUBLE_EQ(mu.dmu_LR_drho_LR, -pLR / (rho * rho));
    }

    {
        // Derivative check away from the branch.
        double const pLR = -3.0e4;
        auto f = [rho, ptol](double const p)
        { return computeYoungLaplaceMacroPotential(p, rho, ptol).mu_LR; };
        auto const mu = computeYoungLaplaceMacroPotential(pLR, rho, ptol);
        EXPECT_NEAR(mu.dmu_LR_dpLR, centralDiff(f, pLR, 1e-3), 1e-11);
    }
}

TEST(RichardsMechanics, PotentialExchangeVanDerWaalsMicroPotential)
{
    double const n_l = 0.03;
    double const rho_lR = 1000.0;
    double const nS = 0.6;
    double const rho_SR = 2700.0;
    double const A = 6.0e-20;
    double const Sa = 1000.0;

    auto const mu =
        computeVanDerWaalsMicroPotential(n_l, rho_lR, nS, rho_SR, A, Sa);

    EXPECT_GT(mu.omega_l, 0.0);
    EXPECT_GT(mu.mu_lR, 0.0);
    EXPECT_LT(mu.dmu_lR_dnl, 0.0);
    EXPECT_GT(mu.dmu_lR_dnS, 0.0);
    EXPECT_GT(mu.dmu_lR_drho_SR, 0.0);
    EXPECT_LT(mu.dmu_lR_drho_lR, 0.0);

    auto f_nl = [=](double const x)
    {
        return computeVanDerWaalsMicroPotential(x, rho_lR, nS, rho_SR, A, Sa)
            .mu_lR;
    };
    auto f_nS = [=](double const x)
    {
        return computeVanDerWaalsMicroPotential(n_l, rho_lR, x, rho_SR, A, Sa)
            .mu_lR;
    };
    auto f_rho_lR = [=](double const x)
    {
        return computeVanDerWaalsMicroPotential(n_l, x, nS, rho_SR, A, Sa)
            .mu_lR;
    };
    auto f_rhoSR = [=](double const x)
    {
        return computeVanDerWaalsMicroPotential(n_l, rho_lR, nS, x, A, Sa)
            .mu_lR;
    };

    EXPECT_NEAR(mu.dmu_lR_dnl, centralDiff(f_nl, n_l, 1e-8), 1e-4);
    EXPECT_NEAR(mu.dmu_lR_drho_lR, centralDiff(f_rho_lR, rho_lR, 1e-4), 1e-8);
    EXPECT_NEAR(mu.dmu_lR_dnS, centralDiff(f_nS, nS, 1e-8), 1e-3);
    EXPECT_NEAR(mu.dmu_lR_drho_SR, centralDiff(f_rhoSR, rho_SR, 1e-4), 1e-6);
}

TEST(RichardsMechanics, PotentialExchangeVanDerWaalsMicroPotentialNegativeAttractiveConvention)
{
    double const n_l = 0.03;
    double const rho_lR = 1000.0;
    double const nS = 0.6;
    double const rho_SR = 2700.0;
    double const A = 6.0e-20;
    double const Sa = 1000.0;

    auto const mu = computeVanDerWaalsMicroPotential(
        n_l, rho_lR, nS, rho_SR, A, Sa, -1.0);

    EXPECT_GT(mu.omega_l, 0.0);
    EXPECT_LT(mu.mu_lR, 0.0);
    EXPECT_GT(mu.dmu_lR_dnl, 0.0);
    EXPECT_LT(mu.dmu_lR_dnS, 0.0);
    EXPECT_LT(mu.dmu_lR_drho_SR, 0.0);
    EXPECT_GT(mu.dmu_lR_drho_lR, 0.0);

    auto f_nl = [=](double const x)
    {
        return computeVanDerWaalsMicroPotential(
                   x, rho_lR, nS, rho_SR, A, Sa, -1.0)
            .mu_lR;
    };
    auto f_nS = [=](double const x)
    {
        return computeVanDerWaalsMicroPotential(
                   n_l, rho_lR, x, rho_SR, A, Sa, -1.0)
            .mu_lR;
    };
    auto f_rho_lR = [=](double const x)
    {
        return computeVanDerWaalsMicroPotential(
                   n_l, x, nS, rho_SR, A, Sa, -1.0)
            .mu_lR;
    };
    auto f_rhoSR = [=](double const x)
    {
        return computeVanDerWaalsMicroPotential(
                   n_l, rho_lR, nS, x, A, Sa, -1.0)
            .mu_lR;
    };

    EXPECT_NEAR(mu.dmu_lR_dnl, centralDiff(f_nl, n_l, 1e-8), 1e-4);
    EXPECT_NEAR(mu.dmu_lR_drho_lR, centralDiff(f_rho_lR, rho_lR, 1e-4), 1e-8);
    EXPECT_NEAR(mu.dmu_lR_dnS, centralDiff(f_nS, nS, 1e-8), 1e-3);
    EXPECT_NEAR(mu.dmu_lR_drho_SR, centralDiff(f_rhoSR, rho_SR, 1e-4), 1e-6);
}

TEST(RichardsMechanics, PotentialExchangeVdWAugmentationOnlyPositiveConvention)
{
    // h = n_l / (nS * rho_SR * Sa), xi = h / lambda
    // mu_lR_aug = K * exp(-xi)
    double const n_l    = 0.03;
    double const rho_lR = 1000.0;
    double const nS     = 0.6;
    double const rho_SR = 2700.0;
    double const Sa     = 1000.0;
    double const K      = 0.5;
    double const lambda = 1.0e-5;  // [m]
    double const A_tiny = std::numeric_limits<double>::min();

    auto const mu = computeVanDerWaalsMicroPotential(
        n_l, rho_lR, nS, rho_SR, A_tiny, Sa,
        /*potential_sign_factor=*/1.0, K, lambda);

    double const xi = n_l / (lambda * nS * rho_SR * Sa);
    double const mu_aug_expected = K * std::exp(-xi);

    EXPECT_NEAR(mu.mu_lR, mu_aug_expected,
                1e-12 * std::abs(mu_aug_expected));
    EXPECT_NEAR(mu.dmu_lR_dnl, -mu_aug_expected * xi / n_l,
                1e-10 * std::abs(mu.dmu_lR_dnl));
    EXPECT_NEAR(mu.dmu_lR_drho_lR, -mu.mu_lR / rho_lR,
                1e-12 * std::abs(mu.dmu_lR_drho_lR));
    EXPECT_NEAR(mu.dmu_lR_dnS, mu_aug_expected * xi / nS,
                1e-10 * std::abs(mu.dmu_lR_dnS));
    EXPECT_NEAR(mu.dmu_lR_drho_SR, mu_aug_expected * xi / rho_SR,
                1e-10 * std::abs(mu.dmu_lR_drho_SR));

    auto f_nl = [=](double const x)
    {
        return computeVanDerWaalsMicroPotential(
                   x, rho_lR, nS, rho_SR, A_tiny, Sa, 1.0, K, lambda)
            .mu_lR;
    };
    EXPECT_NEAR(mu.dmu_lR_dnl, centralDiff(f_nl, n_l, 1e-8), 1e-6);

    auto f_nS = [=](double const x)
    {
        return computeVanDerWaalsMicroPotential(
                   n_l, rho_lR, x, rho_SR, A_tiny, Sa, 1.0, K, lambda)
            .mu_lR;
    };
    EXPECT_NEAR(mu.dmu_lR_dnS, centralDiff(f_nS, nS, 1e-8), 1e-6);

    auto f_rhoSR = [=](double const x)
    {
        return computeVanDerWaalsMicroPotential(
                   n_l, rho_lR, nS, x, A_tiny, Sa, 1.0, K, lambda)
            .mu_lR;
    };
    EXPECT_NEAR(mu.dmu_lR_drho_SR, centralDiff(f_rhoSR, rho_SR, 1e-4), 1e-7);
}

TEST(RichardsMechanics, PotentialExchangeVdWAugmentationCombinedNegativeAttractive)
{
    double const n_l    = 0.03;
    double const rho_lR = 1000.0;
    double const nS     = 0.6;
    double const rho_SR = 2700.0;
    double const A      = 6.0e-20;
    double const Sa     = 1000.0;
    double const K      = 0.4;
    double const lambda = 1.0e-5;  // [m]

    constexpr double pi = 3.141592653589793238462643383279502884;
    double const prefactor_vdW = A * Sa * Sa * Sa / (6.0 * pi);
    double const mu_vdW_expected =
        -prefactor_vdW * (nS * nS * nS) * (rho_SR * rho_SR * rho_SR) /
        (n_l * n_l * n_l);
    double const xi = n_l / (lambda * nS * rho_SR * Sa);
    double const mu_aug_expected = -K * std::exp(-xi);

    auto const mu = computeVanDerWaalsMicroPotential(
        n_l, rho_lR, nS, rho_SR, A, Sa,
        /*potential_sign_factor=*/-1.0, K, lambda);

    EXPECT_NEAR(mu.mu_lR, mu_vdW_expected + mu_aug_expected,
                1e-10 * std::abs(mu.mu_lR));
    EXPECT_LT(mu.mu_lR, 0.0);
    EXPECT_GT(mu.dmu_lR_dnl, 0.0);
    EXPECT_NEAR(mu.dmu_lR_drho_lR, -mu.mu_lR / rho_lR,
                1e-12 * std::abs(mu.dmu_lR_drho_lR));

    auto f_nl = [=](double const x)
    {
        return computeVanDerWaalsMicroPotential(
                   x, rho_lR, nS, rho_SR, A, Sa, -1.0, K, lambda)
            .mu_lR;
    };
    auto f_nS = [=](double const x)
    {
        return computeVanDerWaalsMicroPotential(
                   n_l, rho_lR, x, rho_SR, A, Sa, -1.0, K, lambda)
            .mu_lR;
    };
    auto f_rhoSR = [=](double const x)
    {
        return computeVanDerWaalsMicroPotential(
                   n_l, rho_lR, nS, x, A, Sa, -1.0, K, lambda)
            .mu_lR;
    };

    // Use relative tolerances: the vdW term dominates and produces large
    // derivative values (~50000), so absolute CD truncation error exceeds 1e-6.
    EXPECT_NEAR(mu.dmu_lR_dnl,
                centralDiff(f_nl, n_l, 1e-8),
                1e-7 * std::abs(mu.dmu_lR_dnl));
    EXPECT_NEAR(mu.dmu_lR_dnS,
                centralDiff(f_nS, nS, 1e-8),
                1e-7 * std::abs(mu.dmu_lR_dnS));
    EXPECT_NEAR(mu.dmu_lR_drho_SR,
                centralDiff(f_rhoSR, rho_SR, 1e-4),
                1e-7 * std::abs(mu.dmu_lR_drho_SR));
}

TEST(RichardsMechanics, PotentialExchangeVdWAugmentationZeroCoefficientIsIdentity)
{
    double const n_l = 0.05;
    double const rho_lR = 998.0;
    double const nS = 0.55;
    double const rho_SR = 2650.0;
    double const A = 4.5e-20;
    double const Sa = 800.0;

    auto const mu_baseline =
        computeVanDerWaalsMicroPotential(n_l, rho_lR, nS, rho_SR, A, Sa);

    auto const mu_augmented = computeVanDerWaalsMicroPotential(
        n_l, rho_lR, nS, rho_SR, A, Sa,
        /*potential_sign_factor=*/1.0,
        /*vdw_augmentation_prefactor=*/0.0,
        /*vdw_augmentation_decay_length=*/1.0e-5);

    EXPECT_DOUBLE_EQ(mu_baseline.mu_lR, mu_augmented.mu_lR);
    EXPECT_DOUBLE_EQ(mu_baseline.dmu_lR_dnl, mu_augmented.dmu_lR_dnl);
    EXPECT_DOUBLE_EQ(mu_baseline.dmu_lR_drho_lR, mu_augmented.dmu_lR_drho_lR);
    EXPECT_DOUBLE_EQ(mu_baseline.dmu_lR_dnS, mu_augmented.dmu_lR_dnS);
    EXPECT_DOUBLE_EQ(mu_baseline.dmu_lR_drho_SR, mu_augmented.dmu_lR_drho_SR);
}

TEST(RichardsMechanics, PotentialExchangeSourceIoName)
{
    EXPECT_EQ(ioName(static_cast<MicroExchangeSourceTag*>(nullptr)),
              "micro_exchange_source");
}

TEST(RichardsMechanics, PotentialDrivenMassExchangeDataAndDerivatives)
{
    double const alpha_M = 2.5e-10;
    double const mu_LR = -8.0;
    double const mu_lR = 3.0;

    auto const ex = computePotentialDrivenMassExchange(alpha_M, mu_LR, mu_lR);
    EXPECT_DOUBLE_EQ(ex.rho_l_hat, alpha_M * (mu_LR - mu_lR));
    EXPECT_DOUBLE_EQ(ex.rho_L_hat, -ex.rho_l_hat);
    EXPECT_DOUBLE_EQ(ex.drho_l_hat_dmu_LR, alpha_M);
    EXPECT_DOUBLE_EQ(ex.drho_l_hat_dmu_lR, -alpha_M);
    EXPECT_DOUBLE_EQ(ex.drho_l_hat_dalpha_M, mu_LR - mu_lR);
}

TEST(RichardsMechanics, PotentialExchangeDirectExchangeAlgebra)
{
    double const alpha_bar = 2.5e-10;
    double const mu_LR = 0.0;
    double const mu_lR = -3.0;

    auto const ex = computePotentialDrivenMassExchange(alpha_bar, mu_LR,
                                                       mu_lR);
    EXPECT_DOUBLE_EQ(ex.rho_l_hat, alpha_bar * (mu_LR - mu_lR));
    EXPECT_DOUBLE_EQ(ex.rho_L_hat, -ex.rho_l_hat);
    EXPECT_GT(ex.rho_l_hat, 0.0);
}
