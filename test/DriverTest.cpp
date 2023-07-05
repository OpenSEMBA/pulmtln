#include <gtest/gtest.h>

#include "TestUtils.h"

#include "constants.h"
#include "Driver.h"

using namespace pulmtln;

using json = nlohmann::json;

class DriverTest : public ::testing::Test {};

TEST_F(DriverTest, empty_coax)
{
	// Empty Coaxial case.
	const std::string CASE{ "empty_coax" };
	auto fn{ casesFolder() + CASE + "/" + CASE + ".pulmtln.in.json" };
	
	auto out{ Driver::loadFromFile(fn).getMTLPUL() };

	auto CExpected{ EPSILON0_SI * 2 * M_PI / log(0.05 / 0.025) };
	ASSERT_EQ(1, out.C.NumCols() * out.C.NumRows());
	EXPECT_EQ(CExpected, out.C(0, 0));
	
	auto LExpected{ EPSILON0_SI * MU0_SI / CExpected };
	ASSERT_EQ(1, out.L.NumCols() * out.L.NumRows());
	EXPECT_EQ(LExpected, out.L(0, 0));
}

