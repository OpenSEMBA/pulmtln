
#include <gtest/gtest.h>

#include "core/math/function/Gaussian.h"
#include "core/math/SpaceGenerator.h"

using namespace std;

using namespace semba;
using namespace math;
using namespace function;

TEST(MathFunctionGaussianTest, fourierDC) {
    const Real spread = 0.5;
    const Real delay  = 5;
    const Real amplitude = 1.0;
    Gaussian gaussian(spread, delay, amplitude);
    EXPECT_EQ(0.0, gaussian.getFourier(0.0).imag());
    EXPECT_EQ(spread*sqrt(2.0 * Constants::pi),
            gaussian.getFourier(0.0).real());
}

TEST(MathFunctionGaussianTest, fourierSimmetric) {
    const Real spread = 0.5;
    const Real delay  = 0.0;
    const Real amplitude = 1.0;
    Gaussian gaussian(spread, delay, amplitude);
    vector<Real> fq = linspace(pair<Real,Real>(-10,10), 100);
    for (size_t i = 0; i < fq.size(); ++i) {
        EXPECT_EQ(0.0, gaussian.getFourier(fq[i]).imag());
    }
}
