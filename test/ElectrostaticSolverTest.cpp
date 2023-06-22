#include <gtest/gtest.h>

#include "solver/ElectrostaticSolver.h"

using namespace mfem;
using namespace pulmtln;

class ElectrostaticSolverTest : public ::testing::Test {
};

TEST_F(ElectrostaticSolverTest, parallel_plates)
{
	//    0 V
	//  +-----+
	//  |  3  |
	//  |4   2|
	//  |  1  |
	//  +-----+
	//    1 V
	auto mesh{ 
		Mesh::MakeCartesian2D(11, 11, Element::QUADRILATERAL, 1.0, 1.0)
	};
	BoundaryConditions bcs{{
		{1,    1.0}, // bottom boundary.
		{3,    0.0}, // top boundary.
	}};
	
	SolverOptions opts;
	opts.order = 3;

	ElectrostaticSolver s(mesh, bcs, opts);
	s.Assemble();
	s.Solve();

	//ParaViewDataCollection paraview_dc{ "PULMTLN", &model_.mesh };
	//electrostaticSolver_->writeParaViewFields(paraview_dc);

	const double tol{ 1e-5 };
	EXPECT_NEAR(0.0, s.computeTotalChargeFromRho(), tol);
	EXPECT_NEAR(0.0, s.computeTotalCharge(), tol);

	mfem::Array<int> bdrAttr(1);
	bdrAttr[0] = 1;
	EXPECT_NEAR( 1.0, s.computeChargeInBoundary(bdrAttr), tol);
	
	bdrAttr[0] = 3;
	EXPECT_NEAR(-1.0, s.computeChargeInBoundary(bdrAttr), tol);
	
	bdrAttr[0] = 2;
	EXPECT_NEAR( 0.0, s.computeChargeInBoundary(bdrAttr), tol);
}