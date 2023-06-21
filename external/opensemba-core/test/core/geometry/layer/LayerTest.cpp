#include <gtest/gtest.h>

#include "core/geometry/Layer.h"

using namespace semba;
using namespace geometry;

class GeometryLayersLayerTest : public ::testing::Test {
};

TEST_F(GeometryLayersLayerTest, comparison) {
    Layer lay1P(LayerId(1), "Patata");
    Layer lay1T(LayerId(1), "Tomate");
    Layer lay2T(LayerId(2), "Tomate");

    EXPECT_EQ(lay1P, lay1P);
    EXPECT_NE(lay1P, lay1T);
    EXPECT_NE(lay1P, lay2T);
}




