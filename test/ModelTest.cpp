#include <gtest/gtest.h>

#include "Model.h"
#include "Parser.h"
#include "TestUtils.h"

using namespace pulmtln;

class ModelTest : public ::testing::Test {};

std::string inputCaseFile(const std::string& name) 
{
	return casesFolder() + name + "/" + name + ".pulmtln.in.json";
}

TEST_F(ModelTest, empty_coax_is_not_open)
{
	Model m{ Parser{ inputCaseFile("empty_coax")}.readModel() };
	EXPECT_FALSE(m.isFullyOpen());
}

TEST_F(ModelTest, agrawal1981_is_not_open)
{
	Model m{ Parser{ inputCaseFile("agrawal1981")}.readModel() };
	EXPECT_FALSE(m.isFullyOpen());
}

TEST_F(ModelTest, three_wires_ribbon_is_open)
{
	Model m{ Parser{ inputCaseFile("three_wires_ribbon")}.readModel() };
	EXPECT_TRUE(m.isFullyOpen());
}


