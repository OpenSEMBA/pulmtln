#include <gtest/gtest.h>

#include "ElectrostaticSolver.h"
#include "TestUtils.h"

#include <functional>

using namespace mfem;
using namespace pulmtln;

Vector getBaricenterOfElement(Mesh& mesh, int e)
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

void wireField(const Vector& pos, Vector& res)
{
	// Assumes a charge of 1 unit in the wire.
	double norm = pos.Norml2();
	res = pos;
	res /= (norm * 2.0 * M_PI);
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
	auto m{ Mesh::MakeCartesian2D(1, 5, Element::QUADRILATERAL, 1.0, 1.0) };

	SolverParameters params;
	params.dirichletBoundaries = {
		{
			{1,    1.0}, // bottom boundary.
			{3,    0.0}, // top boundary.
		} 
	};

	ElectrostaticSolver s{m, params};
	s.Solve();

	exportSolution(s, "parallel_plates");

	const double aTol{ 1e-5 };
	const double rTol{ 1e-5 }; // Error percentage of 0.001%

	EXPECT_NEAR(0.0, s.totalChargeFromRho(), aTol);
	EXPECT_NEAR(0.0, s.totalCharge(), aTol);

	EXPECT_LE(relError(1.0, s.chargeInBoundary(1)), rTol);

	EXPECT_LE(relError(-1.0, s.chargeInBoundary(3)), rTol);

	EXPECT_NEAR(0.0, s.chargeInBoundary(2), aTol);

	// Expected energy formula (1/2) * epsilon_0 * E^2 * A.
	// Area A = 1.0 * 1.0
	// Electric field is constant = 1.0 V/m.
	double expectedEnergy{ 0.5 * EPSILON0_NATURAL };
	EXPECT_LE(relError(expectedEnergy, s.totalEnergy()), rTol);

}

TEST_F(ElectrostaticSolverTest, parallel_plates_energy)
{
	//    0 V
	//  +-----+
	//  |  3  |
	//  |4   2|
	//  |  1  |
	//  +-----+
	//    1 V
	auto m{ Mesh::MakeCartesian2D(1, 5, Element::QUADRILATERAL, 1.0, 1.0) };

	SolverParameters params;
	params.dirichletBoundaries = {
		{
			{1,    1.0}, // bottom boundary.
			{3,    0.0}, // top boundary.
		}
	};

	ElectrostaticSolver s{ m, params };
	s.Solve();

	// Expected energy formula (1/2) * epsilon_0 * E^2 * A.
	// Area A = 1.0 * 1.0
	// Electric field is constant = 1.0 V/m.
	const double rTol{ 1e-5 }; // Error percentage of 0.001%
	double expectedEnergy{ 0.5 * EPSILON0_NATURAL };
	EXPECT_LE(relError(expectedEnergy, s.totalEnergy()), rTol);

}

TEST_F(ElectrostaticSolverTest, parallel_plates_neumann)
{
	//    Q = -1 
	//  +-----+
	//  |  3  |
	//  |4   2|
	//  |  1  |
	//  +-----+
	//    Q = 1 
	auto m{ Mesh::MakeCartesian2D(5, 5, Element::QUADRILATERAL, 1.0, 1.0) };

	SolverParameters params;
	params.neumannBoundaries = {
		{
			{1,  1.0}, 
			{3, -1.0}, 
		}
	};

	ElectrostaticSolver s{ m, params };
	s.Solve();

	exportSolution(s, "parallel_plates_neumann");


	const double rTol{ 5e-2 }; 
	EXPECT_LE(relError( 1.0, s.chargeInBoundary(1)), rTol);
    EXPECT_LE(relError(-1.0, s.chargeInBoundary(3)), rTol);

	const double aTol{ 1e-5 };
	EXPECT_NEAR(0.0, s.totalCharge(), aTol);
	EXPECT_NEAR(0.0, s.chargeInBoundary(2), aTol);
	EXPECT_NEAR(0.0, s.chargeInBoundary(4), aTol);

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
	auto m{ Mesh::MakeCartesian2D(1, 5, Element::QUADRILATERAL, 1.0, 1.0) };

	SolverParameters p;
	p.dirichletBoundaries = { {
		{1,    1.0}, // bottom boundary.
		{3,    0.0}, // top boundary.
	} };
	p.domainPermittivities = {{ {1, 2.0} }};

	ElectrostaticSolver s{m, p, SolverOptions{}};
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
	auto m{ Mesh::MakeCartesian2D(1, 10, Element::QUADRILATERAL, 1.0, 1.0) };
	for (auto i{ 0 }; i < m.GetNE(); ++i) {
		auto center{ getBaricenterOfElement(m, i) };
		if (center[1] < 0.5) {
			m.SetAttribute(i, 2);
		}
	}
	m.Finalize();
	m.SetAttributes();

	SolverParameters p;
	p.dirichletBoundaries = { {
		{1, 1.0}, // bottom boundary.
		{3, 0.0}, // top boundary.
	} };
	p.domainPermittivities = { {{2, 4.0}} };

	ElectrostaticSolver s{m, p};
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

	auto fn{ casesFolder() + CASE + "/" + CASE + ".msh" };
	auto m{ Mesh::LoadFromFile(fn.c_str()) };

	const double V{ 1.0 };
	SolverParameters p;
	p.dirichletBoundaries = {{
		{1, 0.0}, // outer boundary
		{2, V},   // inner boundary
	}};
	
	ElectrostaticSolver s{m, p};
	s.Solve();

	exportSolution(s, CASE);

	// Expected capacitance C = eps0 * 2 * pi / log(ro/ri)
	// Expected charge      QExpected = V0 * C 
	double QExpected{ V * EPSILON0_NATURAL * 2 * M_PI / log(0.05 / 0.025) };

	const double rTol{ 5e-3 }; // 0.5% error.

	EXPECT_LE(relError( QExpected, s.chargeInBoundary(2)), rTol); // Boundary 2 is the internal.
	EXPECT_LE(relError(-QExpected, s.chargeInBoundary(1)), rTol); // Boundary 1 is the external.
}

TEST_F(ElectrostaticSolverTest, empty_coax_neumann)
{
	// Coaxial case.
	const std::string CASE{ "empty_coax" };

	auto fn{ casesFolder() + CASE + "/" + CASE + ".msh" };
	auto m{ Mesh::LoadFromFile(fn.c_str()) };

	const double V{ 1.0 };
	SolverParameters p;
	p.neumannBoundaries = { {
		{2, 1.0 / (2 * M_PI * 25e-3)},   // Inner boundary, 1 unit of charge in total.
	} };
	p.dirichletBoundaries = { {
		{1, 0.0},   // Outer boundary. 
	} };

	ElectrostaticSolver s{ m, p };
	s.Solve();

	exportSolution(s, getCaseName());

	const double rTol{ 5e-3 }; // 0.5% error.
	EXPECT_LE(relError( 1.0, s.chargeInBoundary(2)), rTol); // Boundary 2 is the internal.
	EXPECT_LE(relError(-1.0, s.chargeInBoundary(1)), rTol); // Boundary 1 is the external.
}

TEST_F(ElectrostaticSolverTest, empty_coax_neumann_quadrupole)
{
	// Coaxial case.
	const std::string CASE{ "empty_coax" };

	auto fn{ casesFolder() + CASE + "/" + CASE + ".msh" };
	auto m{ Mesh::LoadFromFile(fn.c_str()) };

	SolverParameters p;
	
	std::vector<multipolarCoefficient> ab = {
		{0.0, 0.0}, // a0
		{0.0, 0.0}, // a1, b1
		{1.0, 0.0}, // a2, b2
	};
	
	ElectrostaticSolver s{ m, p };

	// Sets multipolar expansion over internal boundary.
	{
		std::function<double(const Vector&)> f =
			std::bind(&multipolarExpansion, std::placeholders::_1, ab);
		FunctionCoefficient fc(f);
		s.setNeumannCondition(2, fc);
	}

	// Sets opposite multipolar expansion over external boundary.
	{
		ab[2] = { -1.0, 0.0 };
		std::function<double(const Vector&)> f =
			std::bind(&multipolarExpansion, std::placeholders::_1, ab);
		FunctionCoefficient fc(f);
		s.setNeumannCondition(1, fc);
	}

			
	s.Solve();

	exportSolution(s, getTestCaseName());

	const double aTol{ 5e-4 }; 
	EXPECT_NEAR(0.0, s.chargeInBoundary(2), aTol); // Boundary 2 is the internal.
	EXPECT_NEAR(0.0, s.chargeInBoundary(1), aTol); // Boundary 1 is the external.
}

TEST_F(ElectrostaticSolverTest, wire_in_open_region)
{
	// Coaxial case.
	const std::string CASE{ "empty_coax" };

	auto fn{ casesFolder() + CASE + "/" + CASE + ".msh" };
	auto m{ Mesh::LoadFromFile(fn.c_str()) };
	
	SolverParameters p;
	p.dirichletBoundaries = { {
		{2, 1.0}
	} };
	p.openBoundaries = { 1 }; // Outer boundary.

	ElectrostaticSolver s{ m, p };
	s.Solve();

	exportSolution(s, getCaseName());

	double Q = s.chargeInBoundary(2);
	ConstantCoefficient chargeCoeff{ Q };
	VectorFunctionCoefficient exactField{2, wireField, &chargeCoeff};
	auto error = s.GetElectricField().ComputeL2Error(exactField);
	EXPECT_LE(error, 1.2);

	double CExpected = EPSILON0_NATURAL * 2 * M_PI / log(0.05 / 0.025);
	const double rTol{ 1e-2 }; // 1% error.

	auto U{ s.totalEnergy()};
	double CFromEnergy = 0.5 * std::pow(Q,2) / U;
	EXPECT_LE(relError(CExpected, CFromEnergy), rTol);

	double Qb = s.chargeInBoundary(1);
	EXPECT_LE(relError(Q, -Qb), rTol);

	double Vb = s.averagePotentialInBoundary(1);
	auto CFromVb = EPSILON0_NATURAL * Q / (1.0 - Vb);
	EXPECT_LE(relError(CExpected, CFromVb), rTol);

	exportSolution(s, getCaseName());

}

TEST_F(ElectrostaticSolverTest, two_wires_coax)
{
	const std::string CASE{ "two_wires_coax" };

	auto fn{ casesFolder() + CASE + "/" + CASE + ".msh" };
	auto m{ Mesh::LoadFromFile(fn.c_str()) };

	const double V{ 1.0 }; // Voltage
	SolverParameters p;
	p.dirichletBoundaries = { {
		{1, 0.0}, // Conductor 0 bdr (GND).
		{2, V},   // Conductor 1 bdr.
		{3, 0.0}, // Conductor 2 bdr.
	} };

	ElectrostaticSolver s{m, p};
	s.Solve();

	exportSolution(s, getCaseName());

	double C11Expected{  2.15605359 };
	double C12Expected{ -0.16413431 };

	const double rTol{ 2.5e-2 };

	EXPECT_LE(relError(C11Expected, s.chargeInBoundary(2) / V), rTol);
	EXPECT_LE(relError(C12Expected, s.chargeInBoundary(3) / V), rTol);
}

TEST_F(ElectrostaticSolverTest, two_wires_open_capacitance)
{
	const std::string CASE{ "two_wires_open" };
	auto fn{ casesFolder() + CASE + "/" + CASE + ".msh" };
	auto m{ Mesh::LoadFromFile(fn.c_str()) };

	const double V{ 1.0 }; // Voltage
	SolverParameters p;
	p.dirichletBoundaries = { {
		{1,  V}, // Conductor 1 bdr.
		{2, -V}, // Conductor 2 bdr.
	} };
	p.openBoundaries = { 3 };

	ElectrostaticSolver s{ m, p };
	s.Solve();

	exportSolution(s, getCaseName());

	// Capacity between two straight wires.
	double d = 50;
	double rw1 = 2;
	double rw2 = 2;

	double CExpected{ 
		2*M_PI*EPSILON0_NATURAL /
		std::acosh( (d*d - rw1*rw1 - rw2*rw2) / (2*rw1*rw2) )
	};
	
	double chargeInOpenBoundary{ s.chargeInBoundary(3) };
	EXPECT_LE(1e-6, std::abs(chargeInOpenBoundary));

	const double rTol{ 5e-4 };
	double CComputed{ s.chargeInBoundary(1) / (2*V) };
	EXPECT_LE(relError(CExpected, CComputed), rTol);

	double CComputedEnergy{ 2.0 * s.totalEnergy() / (4.0*V*V) };
	EXPECT_LE(relError(CExpected, CComputedEnergy), rTol);
}

TEST_F(ElectrostaticSolverTest, two_wires_open_charges)
{
	const std::string CASE{ "two_wires_open" };
	auto fn{ casesFolder() + CASE + "/" + CASE + ".msh" };
	auto m{ Mesh::LoadFromFile(fn.c_str()) };

	const double V{ 1.0 }; // Voltage
	SolverParameters p;
	p.dirichletBoundaries = { {
		{1,  V}, // Conductor 1 bdr.
		{2,  V}, // Conductor 2 bdr.
	} };
	p.openBoundaries = { 3 };

	ElectrostaticSolver s{ m, p };
	s.Solve();

	exportSolution(s, getTestCaseName());


	auto Q1{ s.chargeInBoundary(1) };
	auto Q2{ s.chargeInBoundary(2) };
	auto Qb{ s.chargeInBoundary(3) };

	auto Vb{ s.averagePotentialInBoundary(3) };
	auto Vd = V - Vb;

	EXPECT_NEAR(0.0, Q1 + Q2 + Qb, 1e-3);
	
	auto CFromEnergy{ 2.0 * s.totalEnergy() / (Vd * Vd) };
	
	auto C1b = Q1 / Vd;
	auto C2b = Q2 / Vd;
	auto CFromCharge{ C1b + C2b };

	EXPECT_LE(relError(CFromEnergy, CFromCharge), 1e-3);
}

TEST_F(ElectrostaticSolverTest, three_wires_ribbon_zero_net_charge)
{
	// Three wires ribbon open problem. 
	// Comparison with Clayton Paul's book:  
	// Analysis of multiconductor transmision lines. 2007.
	// Sec. 5.2.3, p. 187.

	const std::string CASE{ "three_wires_ribbon" };
	auto fn{ casesFolder() + CASE + "/" + CASE + ".msh" };
	auto m{ Mesh::LoadFromFile(fn.c_str()) };

	SolverParameters p;
	p.dirichletBoundaries = {
		{
			{1, 1.0}, // Conductor 0
			{2, 0.0}, // Conductor 1
			{3, 1.0}, // Conductor 2
		}
	};

	SolverOptions solverOpts;
	ElectrostaticSolver s{ m, p, solverOpts };
	s.Solve();

	// For debugging.
	ParaViewDataCollection pd{ outFolder() + CASE + "_zero_net_charge", s.getMesh() };
	s.writeParaViewFields(pd);

	auto Q0 = s.chargeInBoundary(1);
	auto Q1 = s.chargeInBoundary(2);
	auto Q2 = s.chargeInBoundary(3);
	auto Qb = s.chargeInBoundary(4);
	EXPECT_NEAR(0.0, Q0+Q1+Q2+Qb, 1e-4);
}