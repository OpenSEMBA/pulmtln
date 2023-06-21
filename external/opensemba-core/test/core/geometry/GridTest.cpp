#include <gtest/gtest.h>
#include "core/geometry/Grid.h"

using namespace semba;
using namespace geometry;
using namespace math;
using namespace Constants;

class GeometryGridTest : public ::testing::Test {

protected:
	virtual void SetUp() {
		min_ = semba::math::CVecR3(0.0, 0.0, 0.0);
		max_ = semba::math::CVecR3(1.0, 1.0, 1.0);
		step_ = semba::math::CVecR3(0.05, 0.05, 0.05);
		grid_ = semba::geometry::Grid3(semba::geometry::BoxR3(min_, max_), step_);

		std::vector<double> offset_;
		std::vector<double> off_;
		off_.resize(3);
		semba::math::CVecR3 offsetIni;

		offsetIni = grid_.getOrigin();

	}

	semba::geometry::Grid3 grid_;
	semba::math::CVecR3 min_, max_, step_;
};


TEST_F(GeometryGridTest, Initialization) 
{
	Grid3 grid{
		BoxR3(
			math::CVecR3(0.0, 0.0, 0.0),
			math::CVecR3(1, 1, 1)
		),
		math::CVecR3({0.5, 0.5, 0.5})
	};

	Grid3 init = grid;

	EXPECT_EQ(init, grid);
};


TEST_F(GeometryGridTest, BasicOperations) {
    EXPECT_EQ(grid_.getNumCells(), CVecI3(20, 20, 20));
    EXPECT_EQ(Grid3(BoxR3(min_,min_), CVecR3()).getNumCells(), CVecI3(1,1,1));
    EXPECT_EQ(Grid3(BoxR3(min_,max_), CVecR3()).getNumCells(), CVecI3(1,1,1));
}

TEST_F(GeometryGridTest, NaturalCells) {
    const double tol = step_(CartesianAxis::x) / 1.0e3;
    EXPECT_EQ(grid_.getCell(min_, false, 0.0), CVecI3(0,0,0));
    EXPECT_EQ(grid_.getCell(max_, false, tol), grid_.getNumCells());
}

TEST_F(GeometryGridTest, NumberOfCells) {
    Grid3 grid(BoxR3(min_, max_), CVecR3(0.75, 0.75, 0.75));
    EXPECT_EQ(grid.getFullDomainBoundingBox(), BoxR3(min_,CVecR3(1.5,1.5,1.5)));
    EXPECT_EQ(grid.getNumCells(), CVecI3(2,2,2));
}

TEST_F(GeometryGridTest, PosInRange) {
    std::vector<Real> posInRange = grid_.getPosInRange(x,0.17,0.27);
    EXPECT_EQ(2, posInRange.size());
    EXPECT_EQ(0.2, posInRange[0]);
    EXPECT_EQ(0.25, posInRange[1]);
}

TEST_F(GeometryGridTest, OutOfRange) {
    std::vector<Real> posInRange = grid_.getPosInRange(x,1.025,1.5);
    EXPECT_EQ(0, posInRange.size());
}

TEST_F(GeometryGridTest, EnlargeBound) {
        CVecR3 padding(0.15);
        CVecR3 cellSizes(0.05);
        std::pair<CVecR3,CVecR3> pad(padding, padding);
        std::pair<CVecR3,CVecR3> siz(cellSizes, cellSizes);
        CVecI3 numCellsBefore = grid_.getNumCells();
        grid_.enlarge(pad, siz);
        EXPECT_EQ(numCellsBefore + 6, grid_.getNumCells());
}

TEST_F(GeometryGridTest, EnlargeBound2) {
        CVecR3 padding(50.0);
        CVecR3 cellSizes(5);
        std::pair<CVecR3,CVecR3> pad(padding, padding);
        std::pair<CVecR3,CVecR3> siz(cellSizes, cellSizes);
        CVecI3 numCellsBefore = grid_.getNumCells();
        grid_.enlarge(pad, siz);
        EXPECT_EQ(numCellsBefore + 94, grid_.getNumCells());
}

TEST_F(GeometryGridTest, GetSteps) {
    EXPECT_NEAR(0.05, grid_.getStep(0,2), math::tolerance);
    EXPECT_NEAR(0.05, grid_.getStep(1,5), math::tolerance);
    EXPECT_NEAR(0.05, grid_.getStep(2,5), math::tolerance);
}

TEST_F(GeometryGridTest, Equality) {
    EXPECT_EQ(
        grid_, 
        semba::geometry::Grid3(semba::geometry::BoxR3(min_, max_), step_)
    );
}

TEST_F(GeometryGridTest, BoundCoordinates) {
    auto box = grid_.getFullDomainBoundingBox();

	using BoundPairType = std::pair<CartesianAxis, std::vector<CVecR3>>;
    using BoundArrayPairType = std::array<BoundPairType, 3>;

	std::array<std::pair<CartesianBound, BoundArrayPairType>, 2> expectedValues({
		std::pair<CartesianBound, BoundArrayPairType>(
			CartesianBound::L,
			{
				BoundPairType(
					CartesianAxis::x,
					std::vector<CVecR3>({
						CVecR3(0.0, 0.0, 0.0),
						CVecR3(0.0, 0.0, 1.0),
						CVecR3(0.0, 1.0, 1.0),
						CVecR3(0.0, 1.0, 0.0)
					})
				),
				BoundPairType(
					CartesianAxis::y,
					std::vector<CVecR3>({
						CVecR3(0.0, 0.0, 0.0),
						CVecR3(1.0, 0.0, 0.0),
						CVecR3(1.0, 0.0, 1.0),
						CVecR3(0.0, 0.0, 1.0)
					})
				),
				BoundPairType(
					CartesianAxis::z,
					std::vector<CVecR3>({
						CVecR3(0.0, 0.0, 0.0),
						CVecR3(0.0, 1.0, 0.0),
						CVecR3(1.0, 1.0, 0.0),
						CVecR3(1.0, 0.0, 0.0)
					})
				),
			}
		),
		std::pair<CartesianBound, BoundArrayPairType>(
			CartesianBound::U,
			{
				BoundPairType(
					CartesianAxis::x,
					std::vector<CVecR3>({
						CVecR3(1.0, 0.0, 0.0),
						CVecR3(1.0, 0.0, 1.0),
						CVecR3(1.0, 1.0, 1.0),
						CVecR3(1.0, 1.0, 0.0)
					})
				),
				BoundPairType(
					CartesianAxis::y,
					std::vector<CVecR3>({
						CVecR3(0.0, 1.0, 0.0),
						CVecR3(1.0, 1.0, 0.0),
						CVecR3(1.0, 1.0, 1.0),
						CVecR3(0.0, 1.0, 1.0)
					})
				),
				BoundPairType(
					CartesianAxis::z,
					std::vector<CVecR3>({
						CVecR3(0.0, 0.0, 1.0),
						CVecR3(0.0, 1.0, 1.0),
						CVecR3(1.0, 1.0, 1.0),
						CVecR3(1.0, 0.0, 1.0)
					})
				),
			}
		)
	});

    for (const auto& bound : { CartesianBound::L, CartesianBound::U }) {
        for (const auto& axis : { CartesianAxis::x, CartesianAxis::y, CartesianAxis::z }) {
            EXPECT_EQ(
				expectedValues[bound].second[axis].second,
                box.getPosOfBound(axis, bound)
            );
        }
    }
}
