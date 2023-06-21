
#include <gtest/gtest.h>
#include "core/math/simplex/Line.h"

#include <type_traits>

using namespace semba;
using namespace math;

template <typename T>
class MathSimplexLineTest : public ::testing::Test {

};

TEST(MathSimplexLineTest, BasicOperations) {
    static constexpr std::size_t n = 3;
    simplex::Line<n> lin;

    Real sum = 0.0;
    std::vector<Real> weights = lin.getWeights();
    for (size_t i = 0; i < weights.size(); ++i) {
        sum += weights[i];
    }
    EXPECT_NEAR(1.0, sum, 1e-8);
}
