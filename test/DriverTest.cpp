#include <gtest/gtest.h>

#include "Driver.h"
#include "TestUtils.h"

#include <nlohmann/json.hpp>

using namespace pulmtln;

using json = nlohmann::json;

class DriverTest : public ::testing::Test {};

TEST_F(DriverTest, empty_coax)
{
	// Empty Coaxial case.
	const std::string CASE{ "empty_coax" };


	// PhysicalGroups
	const std::map<std::string, int> matToAtt{
		{ "Conductor_0", 1 }, // Outer boundary
		{ "Conductor_1", 2 }, // Inner boundary
		{ "Vacuum", 3 } // Domain
	};

	auto fn{ casesFolder() + CASE + "/" + CASE + ".pulmtln.in.json" };
		
	std::map<int, double> domainToEpsr{};
	
	SolverOptions opts;
	opts.order = 3;

	ElectrostaticSolver s(mesh, bcs, domainToEpsr, opts);
	s.Solve();

	exportSolution(s, CASE);

	// Expected capacitance C = eps0 * 2 * pi / log(ro/ri)
	// Expected charge      QExpected = V0 * C 
	double QExpected{ V0 * epsilon0_ * 2 * M_PI / log(0.05 / 0.025) };

	const double rTol{ 5e-3 }; // 0.5% error.

	mfem::Array<int> bdrAttr(1);

	bdrAttr[0] = matToAtt.at("Conductor_1");
	EXPECT_LE(relError(QExpected, s.chargeInBoundary(bdrAttr)), rTol);

	bdrAttr[0] = matToAtt.at("Conductor_0");
	EXPECT_LE(relError(-QExpected, s.chargeInBoundary(bdrAttr)), rTol);
}

