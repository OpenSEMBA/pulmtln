#include <gtest/gtest.h>

#include "TestUtils.h"

#include "Parameters.h"

using namespace pulmtln;

using json = nlohmann::json;

class ParametersTest : public ::testing::Test {};

TEST_F(ParametersTest, toJSON)
{
	Parameters p;

	p.C = mfem::DenseMatrix(2, 2);
	p.C(0, 0) = 1.0; p.C(0, 1) = 2.0;
	p.C(1, 0) = 3.0; p.C(1, 1) = 4.0;
	
	p.L = mfem::DenseMatrix(2, 2);
	p.L(0, 0) = 1.0; p.L(0, 1) = 2.0;
	p.L(1, 0) = 3.0; p.L(1, 1) = 4.0;

	auto j{ p.toJSON() };

	ASSERT_TRUE(j.contains("C"));
	ASSERT_TRUE(j.contains("L"));
	// TODO

	EXPECT_TRUE(false);
}

TEST_F(ParametersTest, saveToJSONFile)
{
	EXPECT_TRUE(false);
}