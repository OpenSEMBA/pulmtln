#include <gtest/gtest.h>

#include "core/geometry/element/Line2.h"

using namespace semba;

TEST(Line2Test, CanCreate) {
	geometry::CoordR3 vertexLeft{ geometry::CoordId(), math::CVecR3({0.0, 0.0, 0.0}) };
	geometry::CoordR3 vertexRight{ geometry::CoordId(), math::CVecR3({1.0, 0.0, 0.0}) };

	geometry::Layer lay{geometry::LayerId(), "My layer"};
	geometry::element::Model model{MatId()};

	geometry::LinR2 lin{ geometry::ElemId(), {&vertexLeft, &vertexRight}, &lay, &model };

	auto newLines = lin.splitByMiddle();

	EXPECT_EQ(
		*(newLines[0]->getVertices()[0]),
		*lin.getVertices()[0]
	);
	EXPECT_EQ(
		*(newLines[0]->getVertices()[1]),
		geometry::CoordR3(geometry::CoordId(), math::CVecR3(0.5, 0.0, 0.0))
	);
	EXPECT_EQ(&lay, newLines[0]->getLayer());
	EXPECT_EQ(&model, newLines[0]->getModel());

	EXPECT_EQ(
		*(newLines[1]->getVertices()[0]),
		geometry::CoordR3(geometry::CoordId(), math::CVecR3(0.5, 0.0, 0.0))
	);
	EXPECT_EQ(
		*(newLines[1]->getVertices()[1]),
		*lin.getVertices()[1]
	);
	EXPECT_EQ(&lay, newLines[1]->getLayer());
	EXPECT_EQ(&model, newLines[1]->getModel());

	EXPECT_EQ(newLines[0]->getLayer()->getName(), lay.getName());
}
