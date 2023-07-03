#include <gtest/gtest.h>

#include "Parser.h"
#include "TestUtils.h"


using namespace pulmtln;


class ParserTest : public ::testing::Test {};

TEST_F(ParserTest, empty_coax)
{
	const std::string CASE{ "empty_coax" };
	Parser parser{ 
		readJSON(casesFolder() + CASE + "/" + CASE + ".pulmtln.in.json") 
	};

	auto opts{ parser.readSolverOptions() };
	EXPECT_EQ(3, opts.order);
	EXPECT_EQ(true, opts.exportParaViewSolution);

	auto model{ parser.readModel() };
	auto pecs{ model.getMaterialsOfType(MaterialType::PEC) };
	ASSERT_EQ(2, pecs.size());
	EXPECT_EQ(1, pecs.at("Conductor_0"));
	EXPECT_EQ(2, pecs.at("Conductor_1"));

	EXPECT_NE(0, model.getMesh()->GetNE());
	EXPECT_EQ(2, model.getMesh()->bdr_attributes.Size());
}

