#include <gtest/gtest.h>

#include "ElectrostaticSolver.h"
#include "TestUtils.h"

using namespace mfem;
using namespace pulmtln;

mfem::Vector getBaricenterOfElement(Mesh& mesh, int e)
{
	Element* elem{ mesh.GetElement(e) };
	mfem::Vector center({ 0.0, 0.0 });
	mfem::Array<int> vs(elem->GetNVertices());
	elem->GetVertices(vs);
	for (auto i{ 0 }; i < vs.Size(); i++) {
		mfem::Vector vertexPos(mesh.GetVertex(vs[i]), 2);
		center += vertexPos;
	}
	center /= vs.Size();
	return center;
}


void exportSolution(ElectrostaticSolver& s, const std::string& caseName)
{
	ParaViewDataCollection paraview_dc{ outFolder() + caseName, s.getMesh() };
	s.writeParaViewFields(paraview_dc);
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
		Mesh::MakeCartesian2D(1, 5, Element::QUADRILATERAL, 1.0, 1.0)
	};
	AttrToValueMap bcs{ {
		{1,    1.0}, // bottom boundary.
		{3,    0.0}, // top boundary.
	} };

	SolverOptions opts;
	opts.order = 3;


	ElectrostaticSolver s(mesh, bcs, {}, opts);
	s.Solve();

	exportSolution(s, "parallel_plates");

	const double aTol{ 1e-5 };
	const double rTol{ 1e-5 }; // Error percentage of 0.001%

	EXPECT_NEAR(0.0, s.totalChargeFromRho(), aTol);
	EXPECT_NEAR(0.0, s.totalCharge(), aTol);

	EXPECT_LE(relError(1.0, s.chargeInBoundary(1)), rTol);

	EXPECT_LE(relError(-1.0, s.chargeInBoundary(3)), rTol);

	EXPECT_NEAR(0.0, s.chargeInBoundary(2), aTol);

}

TEST_F(ElectrostaticSolverTest, parallel_plates_epsr2)
{
	//    0 V
	//  +-----+
	//  |  3  |
	//  |4   2|
	//  |  1  |
	//  +-----+
	//    1 V
	auto mesh{
		Mesh::MakeCartesian2D(1, 5, Element::QUADRILATERAL, 1.0, 1.0)
	};
	AttrToValueMap bcs{ {
		{1,    1.0}, // bottom boundary.
		{3,    0.0}, // top boundary.
	} };

	SolverOptions opts;
	opts.order = 3;

	std::map<int, double> domain2eps{
		{1, 2.0}
	};
	ElectrostaticSolver s(mesh, bcs, domain2eps, opts);
	s.Solve();

	exportSolution(s, "Parallel_plates_epsr2");

	const double tol{ 1e-6 };
	EXPECT_LE(relError(2.0, s.chargeInBoundary(1)), tol);
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
		Mesh::MakeCartesian2D(1, 10, Element::QUADRILATERAL, 1.0, 1.0)
	};
	for (auto i{ 0 }; i < mesh.GetNE(); ++i) {
		auto center{ getBaricenterOfElement(mesh, i) };
		if (center[1] < 0.5) {
			mesh.SetAttribute(i, 2);
		}
	}
	mesh.Finalize();
	mesh.SetAttributes();

	AttrToValueMap bcs{ {
		{1, 1.0}, // bottom boundary.
		{3, 0.0}, // top boundary.
	} };

	std::map<int, double> domainAttributeToEpsr{
		{2, 4.0}
	};

	SolverOptions opts;
	opts.order = 3;

	ElectrostaticSolver s(mesh, bcs, domainAttributeToEpsr, opts);
	s.Solve();

	exportSolution(s, "two_materials");

	const double aTol{ 1e-4 };
	const double rTol{ 1e-5 };  // 0.001%

	EXPECT_NEAR(0.0, s.totalChargeFromRho(), aTol);
	EXPECT_NEAR(0.0, s.totalCharge(), aTol);

	EXPECT_LE(relError(1.6, s.chargeInBoundary(1)), rTol);
	EXPECT_LE(relError(- 1.6, s.chargeInBoundary(3)), rTol);

}

TEST_F(ElectrostaticSolverTest, empty_coax)
{
	// Coaxial case.
	const std::string CASE{ "empty_coax" };
	// PhysicalGroups
	const std::map<std::string, int> matToAtt{
		{ "Conductor_0", 1 }, // Outer boundary
		{ "Conductor_1", 2 }, // Inner boundary
		{ "Vacuum", 3 } // Domain
	};

	auto fn{ casesFolder() + CASE + "/" + CASE + ".msh" };
	auto mesh{ Mesh::LoadFromFile(fn.c_str()) };

	const double V{ 1.0 };
	AttrToValueMap bcs{{
		{1, 0.0}, // outer boundary
		{2, V},   // inner boundary
	}};
	
	std::map<int, double> domainToEpsr{};
	
	SolverOptions opts;
	opts.order = 3;

	ElectrostaticSolver s(mesh, bcs, domainToEpsr, opts);
	s.Solve();

	exportSolution(s, CASE);

	// Expected capacitance C = eps0 * 2 * pi / log(ro/ri)
	// Expected charge      QExpected = V0 * C 
	double QExpected{ V * EPSILON0_NATURAL * 2 * M_PI / log(0.05 / 0.025) };

	const double rTol{ 5e-3 }; // 0.5% error.

	EXPECT_LE(relError(
		QExpected, 
		s.chargeInBoundary(matToAtt.at("Conductor_1"))
	), rTol);

	EXPECT_LE(relError(
		-QExpected, 
		s.chargeInBoundary(matToAtt.at("Conductor_0"))
	), rTol);
}


TEST_F(ElectrostaticSolverTest, two_wires_coax)
{
	const std::string CASE{ "two_wires_coax" };

	// PhysicalGroups
	const std::map<std::string, int> matToAtt{
		{ "Conductor_0", 1 }, // Outer boundary
		{ "Conductor_1", 2 }, // Inner boundary
		{ "Conductor_2", 3 }, // Inner boundary
		{ "Vacuum",      4 } // Domain
	};

	auto fn{ casesFolder() + CASE + "/" + CASE + ".msh" };
	auto mesh{ Mesh::LoadFromFile(fn.c_str()) };

	const double V{ 1.0 }; // Voltage
	AttrToValueMap bcs{ {
		{1, 0.0}, // Conductor 0 bdr (GND).
		{2, V},   // Conductor 1 bdr.
		{3, 0.0}, // Conductor 2 bdr.
	} };

	std::map<int, double> domainToEpsr{};

	SolverOptions opts;
	opts.order = 3;

	ElectrostaticSolver s(mesh, bcs, domainToEpsr, opts);
	s.Solve();

	exportSolution(s, CASE);

	double CMat11Expected{  2.15605359 };
	double CMat12Expected{ -0.16413431 };

	const double rTol{ 2.5e-2 };

	EXPECT_LE(
		relError(
			CMat11Expected, 
			s.chargeInBoundary(matToAtt.at("Conductor_1")) / V
		), rTol);
	
	EXPECT_LE(
		relError(
			CMat12Expected, 
			s.chargeInBoundary(matToAtt.at("Conductor_2")) / V
		), rTol);
}