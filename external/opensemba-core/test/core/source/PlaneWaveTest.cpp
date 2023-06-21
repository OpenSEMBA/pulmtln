
#include <gtest/gtest.h>

#include "core/source/PlaneWave.h"
#include "core/math/function/Gaussian.h"

using namespace SEMBA;
using namespace Source;
using namespace Geometry;

class EMSourcePlaneWaveTest : public ::testing::Test {
    void SetUp() {

    }

    void TearDown() {

    }

protected:
    static const std::unique_ptr<Magnitude::Magnitude> buildMagnitude() {
        return std::make_unique<Magnitude::Magnitude>(Magnitude::Magnitude(
            new math::Function::Gaussian(
                math::Function::Gaussian::buildFromMaximumFrequency(
                    1e9,
                    1.0
                )
            )
        ));
    }
};

TEST_F(EMSourcePlaneWaveTest, PolarCoordinatesDirAndPolarization) {
    {
        math::CVecR3 dir(1.0, 0.0, 0.0);
        math::CVecR3 pol(0.0, 0.0, 1.0);
        PlaneWave pw(buildMagnitude(), {}, dir, pol);
        EXPECT_NEAR(math::Constants::pi_2, pw.getTheta(), 1e-3);
        EXPECT_NEAR(0.0,                   pw.getPhi(),   1e-3);
        EXPECT_NEAR(0.0,                   pw.getAlpha(), 1e-3);
        EXPECT_NEAR(0.0,                   pw.getBeta(),  1e-3);
    }

    {
        math::CVecR3 dir(-1.0, 0.0, 0.0);
        math::CVecR3 pol( 0.0, 1.0, 0.0);
        PlaneWave pw(buildMagnitude(), {}, dir, pol);
        EXPECT_NEAR(math::Constants::pi_2, pw.getTheta(), 1e-3);
        EXPECT_NEAR(math::Constants::pi, pw.getPhi(), 1e-3);
        EXPECT_NEAR(math::Constants::pi_2, pw.getAlpha(), 1e-3);
        EXPECT_NEAR(math::Constants::pi_2, pw.getBeta(), 1e-3);
    }

    {
        math::CVecR3 dir(1.0, 0.0, 0.0);
        math::CVecR3 pol(0.0, 1.0, 0.0);
        PlaneWave pw(buildMagnitude(), {}, dir, pol);
        EXPECT_NEAR(math::Constants::pi_2, pw.getTheta(), 1e-3);
        EXPECT_NEAR(0.0, pw.getPhi(), 1e-3);
        EXPECT_NEAR(math::Constants::pi_2, pw.getAlpha(), 1e-3);
        EXPECT_NEAR(math::Constants::pi_2, pw.getBeta(), 1e-3);
    }

    {
        math::CVecR3 dir(0.0, 0.0, -1.0);
        math::CVecR3 pol(1.0, 1.0, 0.0);
        PlaneWave pw(buildMagnitude(), {}, dir, pol);
        EXPECT_NEAR(math::Constants::pi, pw.getTheta(), 1e-3);
        EXPECT_NEAR(0.0, pw.getPhi(), 1e-3);
        EXPECT_NEAR(math::Constants::pi_2, pw.getAlpha(), 1e-3);
        EXPECT_NEAR(math::Constants::pi / 4.0, pw.getBeta(), 1e-3);
    }

    {
        math::CVecR3 dir(0.0, 0.0, -1.0);
        math::CVecR3 pol(1.0, 1.0, 0.0);
        PlaneWave pw(buildMagnitude(), {}, dir, pol);
        EXPECT_NEAR(math::Constants::pi, pw.getTheta(), 1e-3);
        EXPECT_NEAR(0.0, pw.getPhi(), 1e-3);
        EXPECT_NEAR(math::Constants::pi_2, pw.getAlpha(), 1e-3);
        EXPECT_NEAR(math::Constants::pi / 4.0, pw.getBeta(), 1e-3);
    }

    {
        math::CVecR3 dir(0.0, 0.0, -1.0);
        math::CVecR3 pol(-1.0, 1.0, 0.0);
        PlaneWave pw(buildMagnitude(), {}, dir, pol);
        EXPECT_NEAR(math::Constants::pi,         pw.getTheta(), 1e-3);
        EXPECT_NEAR(0.0,                         pw.getPhi(),   1e-3);
        EXPECT_NEAR(math::Constants::pi_2,       pw.getAlpha(), 1e-3);
        EXPECT_NEAR(math::Constants::pi*3.0/4.0, pw.getBeta(),  1e-3);
    }
}

