#include <gtest/gtest.h>

#include "solver/ElectrostaticSolver.h"

using namespace mfem;
using namespace pulmtln;

mfem::Vector getBaricenterOfElement(Mesh& mesh, int i)
{
	Element* elem{ mesh.GetElement(i) };
	mfem::Vector center({ 0.0, 0.0 });
	mfem::Array<int> vertices(elem->GetNVertices());
	elem->GetVertices(vertices);
	for (auto v{ 0 }; v < vertices.Size(); v++) {
		mfem::Vector vertexPos(mesh.GetVertex(i), 2);
		center += vertexPos;
	}
	center /= vertices.Size();
	return center;
}


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
		Mesh::MakeCartesian2D(5, 5, Element::QUADRILATERAL, 1.0, 1.0)
	};
	BoundaryConditions bcs{{
		{1,    1.0}, // bottom boundary.
		{3,    0.0}, // top boundary.
	}};
	
	SolverOptions opts;
	opts.order = 3;

	const double tol{ 1e-5 };
	{
		ElectrostaticSolver s(mesh, bcs, {}, opts);
		s.Solve();

		//ParaViewDataCollection paraview_dc{ "Parallel_plates", &mesh };
		//s.writeParaViewFields(paraview_dc);

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
	{
		std::map<int, double> d2eps{
			{1, 2.0}
		};
		ElectrostaticSolver s2(mesh, bcs, d2eps, opts);
		s2.Solve();

		mfem::Array<int> bdrAttr(1);
		bdrAttr[0] = 1;
		EXPECT_NEAR(2.0, s2.computeChargeInBoundary(bdrAttr), tol);
	}
}

TEST_F(ElectrostaticSolverTest, two_materials)
{
	//         0 V,   Bdr 3
	//       +-------------+
	//       | eps_r=1, D1 |
	// Bdr 4 |-------------| Bdr 2
	//       | eps_r=4, D2 |
	//       +-------------+
	//         1 V,   Bdr 1
	auto mesh{
		Mesh::MakeCartesian2D(11, 10, Element::QUADRILATERAL, 1.0, 1.0)
	};
	for (auto i{ 0 }; i < mesh.GetNE(); ++i) {
		auto center{ getBaricenterOfElement(mesh, i) };
		if (center[1] < 0.5) {
			mesh.SetAttribute(i, 2);
		}
	}
	mesh.Finalize();
	mesh.SetAttributes();

	BoundaryConditions bcs{ {
		{1, 1.0}, // bottom boundary.
		{3, 0.0}, // top boundary.
	} };

	std::map<int, double> domainToEpsr{
		{2, 100.0}
	};

	SolverOptions opts;
	opts.order = 3;

	ElectrostaticSolver s(mesh, bcs, domainToEpsr, opts);
	s.Solve();

	ParaViewDataCollection paraview_dc{ "two_materials", &mesh };
	s.writeParaViewFields(paraview_dc);

	//const double tol{ 1e-5 };
	//EXPECT_NEAR(0.0, s.computeTotalChargeFromRho(), tol);
	//EXPECT_NEAR(0.0, s.computeTotalCharge(), tol);

	//mfem::Array<int> bdrAttr(1);
	//bdrAttr[0] = 1;
	//EXPECT_NEAR(1.0, s.computeChargeInBoundary(bdrAttr), tol);

	//bdrAttr[0] = 3;
	//EXPECT_NEAR(-1.0, s.computeChargeInBoundary(bdrAttr), tol);

	//bdrAttr[0] = 2;
	//EXPECT_NEAR(0.0, s.computeChargeInBoundary(bdrAttr), tol);
}