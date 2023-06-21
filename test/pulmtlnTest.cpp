#include <gtest/gtest.h>

#include "solver/Solver.h"


using namespace pulmtln;
using namespace mfem;
using namespace fixtures::sources;

class pulmtlnTest : public ::testing::Test {
};

TEST_F(pulmtlnTest, parallel_plates)
{
	pulmtln::Solver solver{
	buildModel(14,1,Element::Type::TRIANGLE, 1.0, 1.0, BdrCond::PMC,BdrCond::PEC,BdrCond::PMC,BdrCond::PEC),
	probes,
	buildGaussianInitialField(E, 0.1, mfem::Vector({0.5,0.5})),
	SolverOptions{}
		.setTimeStep(1e-2)
		.setCentered()
		.setFinalTime(2.0)
		.setOrder(3)
	}; 

	auto normOld{ solver.getFields().getNorml2() };
	solver.run();

	double tolerance{ 1e-2 };
	EXPECT_NEAR(normOld, solver.getFields().getNorml2(), tolerance);

	//At the left boundary the electric field should be closed to zero and
	//the magnetic field reaches a maximum close to 1.0 or -1.0
	//(the wave splits in two and doubles at the boundary).
	EXPECT_NEAR(0.0, abs(solver.getPointProbe(0).findFrameWithMax().second), tolerance);
	EXPECT_NEAR(0.0, abs(solver.getPointProbe(1).findFrameWithMax().second), tolerance);
	EXPECT_NEAR(1.0, abs(solver.getPointProbe(2).findFrameWithMax().second), tolerance);
	EXPECT_NEAR(1.0, abs(solver.getPointProbe(3).findFrameWithMax().second), tolerance);
}
