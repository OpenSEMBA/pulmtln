#include <gtest/gtest.h>
#include "core/math/simplex/Triangle.h"

using namespace semba;
using namespace math;

template <typename T>
class MathSimplexTriangleTest : public ::testing::Test {

};

TEST(MathSimplexTriangleTest, BasicOperations) {
    static constexpr std::size_t n = 3;
    simplex::Triangle<n> tri;

    Real sum = 0.0;
    std::vector<Real> weights = tri.getWeights();
    for (size_t i = 0; i < weights.size(); ++i) {
        sum += weights[i];
    }
    EXPECT_NEAR(1.0, sum, 1e-8);
}
