#include <gtest/gtest.h>

#include "core/physicalModel/wire/Extremes.h"
#include "core/physicalModel/multiport/Predefined.h"

namespace semba {
namespace physicalModel {

class ExtremesTest : public ::testing::Test {
public:

protected:

};

TEST_F(ExtremesTest, build_from_wire_and_multiport)
{
    wire::Wire wire(Id(2), "Cable", 1.0, 1.0, 1.0);
    multiport::Predefined mp1(Id(3), "short", multiport::Multiport::Type::shortCircuit);
    multiport::Predefined mp2(Id(4), "open", multiport::Multiport::Type::openCircuit);
    
    wire::Extremes extremes("extremes", wire, &mp1, &mp2);

    EXPECT_EQ(Id(2), extremes.getId());
    EXPECT_EQ("extremes", extremes.getName()); 

    EXPECT_EQ(mp1.getType(), extremes.getExtreme(0)->getType());
    EXPECT_EQ(mp2.getType(), extremes.getExtreme(1)->getType());

    wire::Extremes copied(extremes);
    EXPECT_EQ(Id(2),         copied.getId());
    EXPECT_EQ("extremes",    copied.getName());
    EXPECT_EQ(mp1.getType(), copied.getExtreme(0)->getType());
    EXPECT_EQ(mp2.getType(), copied.getExtreme(1)->getType());

}


}
}