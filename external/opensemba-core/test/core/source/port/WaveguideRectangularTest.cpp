
#include <gtest/gtest.h>

#include "core/source/port/WaveguideRectangular.h"
#include "core/geometry/element/Quadrilateral4.h"
#include "core/physicalModel/Group.h"
#include "core/physicalModel/Bound.h"
#include "core/math/function/Gaussian.h"

using namespace semba;
using namespace source;
using namespace physicalModel;
using namespace std;

class SourcePortWaveguideRectangularTest : public ::testing::Test {
    void SetUp() {
        geometry::BoxI3 plane(math::CVecI3(0,0,0), math::CVecI3(30,10,0));
        vector<geometry::BoxI3> quadBoxes = plane.chop();
        geometry::ElemId id(0);
        for (size_t i = 0; i < quadBoxes.size(); i++) {
            std::vector<const geometry::CoordI3*> coords;
            for (size_t j = 0; j < geometry::Qua4::sizeOfCoordinates; j++) {
                coords.push_back(
                    cG_.addAndAssignId(
                        std::make_unique<geometry::CoordI3>(
                            geometry::CoordId(),
                            quadBoxes[i].getPos()[j]
                            )
                    )->get()
                );
            }

            surfs.add(
                std::make_unique<geometry::QuaI4>(++id, coords.data())
            );
        }

        excMode = port::Waveguide::ExcitationMode::TE;
        mode = pair<size_t,size_t>(1,0);

        boundType.add( std::make_unique <Bound>(MatId(1), Bound::Type::pml));
        boundType.add( std::make_unique <Bound>(MatId(2), Bound::Type::pmc));
        boundType.add( std::make_unique <Bound>(MatId(3), Bound::Type::pec));

        bounds[0][0] = boundType.getId(MatId(1))->castTo<Bound>();
        bounds[1][1] = boundType.getId(MatId(1))->castTo<Bound>();
        bounds[2][0] = boundType.getId(MatId(1))->castTo<Bound>();
        bounds[0][1] = boundType.getId(MatId(1))->castTo<Bound>();
        bounds[1][0] = boundType.getId(MatId(1))->castTo<Bound>();
        bounds[2][1] = boundType.getId(MatId(1))->castTo<Bound>();
    }

protected:
    geometry::CoordI3Group cG_;
    geometry::element::Group<geometry::Surf> surfs;
    port::Waveguide::ExcitationMode excMode;
    pair<size_t,size_t> mode;
    port::Bound3 bounds;
    PMGroup boundType;

    static const std::unique_ptr<Magnitude::Magnitude> buildMagnitude() {
        return std::make_unique<Magnitude::Magnitude>(Magnitude::Magnitude(
            new math::function::Gaussian(
                math::function::Gaussian::buildFromMaximumFrequency(
                    1e9,
                    1.0
                )
            )
        ));
    }
};

//TEST_F(SourcePortWaveguideRectangularTest, withoutSymmetries) 
//{
//    Port::WaveguideRectangular wp{
//        buildMagnitude(), 
//        {surfs.getOf<Geometry::Elem>().front()->getId()},
//        excMode, 
//        mode
//    };
//
//    EXPECT_EQ(math::CVecR3(0.0), wp.getOrigin());
//    EXPECT_EQ(30.0, wp.getWidth());
//    EXPECT_EQ(10.0, wp.getHeight());
//
//    math::CVecR3 midPoint = wp.getOrigin();
//    midPoint(math::Constants::x) = wp.getWidth() / 2.0;
//
//    EXPECT_EQ(math::CVecR3(0.0,1.0,0.0), wp.getWeight(midPoint));
//
//    EXPECT_EQ(math::CVecR3(0.0), wp.getWeight(wp.getOrigin()));
//}