
#include <gtest/gtest.h>

#include "core/math/CartesianVector.h"

using namespace semba;
using namespace math;

TEST(MathRealTest, equal) {
    std::vector<Real> num;
    num.push_back((Real) -2.69386e-10);
    num.push_back((Real) 1e-10);

    for (UInt i = 0; i < num.size(); i++) {
        EXPECT_TRUE(math::equal(num[i], num[i]));
        EXPECT_TRUE(math::equal(num[i], num[i], 0.0, 0.25));
    }

    EXPECT_FALSE(math::equal(0.0, 0.17, 0.5, 0.01));
}
