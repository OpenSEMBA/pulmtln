#include <gtest/gtest.h>

#include "Parser.h"
#include "TestUtils.h"

using namespace pulmtln;

class ParserTest : public ::testing::Test {};

TEST_F(ParserTest, empty_coax)
{
	const std::string CASE{ "empty_coax" };
	Parser parser{ 
		casesFolder() + CASE + "/" + CASE + ".pulmtln.in.json" 
	};

	auto opts{ parser.readDriverOptions() };
	EXPECT_EQ(3, opts.solverOptions.order);
	EXPECT_EQ(true, opts.exportParaViewSolution);

	auto model{ parser.readModel() };
	const auto& pecs{ model.getMaterials().pecs };
	ASSERT_EQ(2, pecs.size());
	EXPECT_EQ(pecs[0].name, "Conductor_0");
	EXPECT_EQ(pecs[0].tag, 1);
	EXPECT_EQ(pecs[1].name, "Conductor_1");
	EXPECT_EQ(pecs[1].tag, 2);

	const auto& diel{ model.getMaterials().dielectrics };
	ASSERT_EQ(0, diel.size());

	EXPECT_NE(0, model.getMesh()->GetNE());
	EXPECT_EQ(2, model.getMesh()->bdr_attributes.Size());
}

TEST_F(ParserTest, partially_filled_coax)
{
	const std::string CASE{ "partially_filled_coax" };
	Parser parser{
		casesFolder() + CASE + "/" + CASE + ".pulmtln.in.json"
	};

	auto opts{ parser.readDriverOptions() };
	EXPECT_EQ(3, opts.solverOptions.order);
	EXPECT_EQ(true, opts.exportParaViewSolution);

	auto model{ parser.readModel() };
	const auto& pecs{ model.getMaterials().pecs };
	ASSERT_EQ(2, pecs.size());
	EXPECT_EQ(pecs[0].name, "Conductor_0");
	EXPECT_EQ(pecs[0].tag, 1);
	EXPECT_EQ(pecs[1].name, "Conductor_1");
	EXPECT_EQ(pecs[1].tag, 2);

	const auto& diels{ model.getMaterials().dielectrics };
	ASSERT_EQ(1, diels.size());
	EXPECT_EQ(diels[0].name, "Dielectric_1");
	EXPECT_EQ(diels[0].tag, 4);
	EXPECT_EQ(diels[0].relativePermittivity, 4.0);

	EXPECT_NE(0, model.getMesh()->GetNE());
	EXPECT_EQ(2, model.getMesh()->bdr_attributes.Size());
}

