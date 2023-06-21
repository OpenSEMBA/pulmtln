#include <gtest/gtest.h>

#include "core/source/port/TEMCoaxial.h"
#include "core/geometry/element/Quadrilateral4.h"
#include "core/math/function/Gaussian.h"

using namespace semba;
using namespace source;
using namespace std;

class SourcePortTEMCoaxialTest : public ::testing::Test {
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

        excMode = port::TEM::voltage;
        innerRadius_ = 1.0;
        outerRadius_ = 3.0;
    }

protected:
    geometry::CoordI3Group cG_;
    geometry::element::Group<geometry::Surf> surfs;
    port::TEM::ExcitationMode excMode;
    math::Real innerRadius_, outerRadius_;

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

//TEST_F(SourcePortTEMCoaxialTest, basic) {
//    Port::TEMCoaxial port(
//        buildMagnitude(),
//        {surfs.getOf<Geometry::Elem>().front()->getId()},
//        excMode,
//        math::CVecR3(0.0), 
//        innerRadius_, 
//        outerRadius_
//    );
//
//    EXPECT_EQ(math::CVecR3(0.0,1.0,0.0),
//            port.getWeight(math::CVecR3(0.0,innerRadius_,0.0)).normalize());
//    EXPECT_EQ(math::CVecR3(1.0,0.0,0.0),
//            port.getWeight(math::CVecR3(innerRadius_,0.0,0.0)).normalize());
//}
