#include <gtest/gtest.h>

#include "TestUtils.h"

#include "Parameters.h"

using namespace pulmtln;

using json = nlohmann::json;

class ParametersTest : public ::testing::Test {
public:
	void SetUp() {
		p.C = mfem::DenseMatrix(2, 2);
		p.C(0, 0) = 1.0; p.C(0, 1) = 2.0;
		p.C(1, 0) = 3.0; p.C(1, 1) = 4.0;

		p.L = mfem::DenseMatrix(2, 2);
		p.L(0, 0) = 1.0; p.L(0, 1) = 2.0;
		p.L(1, 0) = 3.0; p.L(1, 1) = 4.0;
	}

protected:
	PULParameters p;
};

TEST_F(ParametersTest, toJSON)
{
	auto j{ p.toJSON() };
	
	auto r{ j };

	EXPECT_EQ(p, r);
}

TEST_F(ParametersTest, saveToJSONFile)
{
	std::string fn = "to_from_JSONFile_test.pulmtln.out.json";

	p.saveToJSONFile(fn);

	json j;
	std::ifstream ifs{ fn };
	j << ifs;

	PULParameters r{ j };

	EXPECT_EQ(p, r);

	
}