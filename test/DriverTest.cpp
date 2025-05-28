#include "gtest/gtest.h"

#include "TestUtils.h"

#include "constants.h"
#include "Driver.h"
#include "ElectrostaticSolver.h"

using namespace pulmtln;

using json = nlohmann::json;

class DriverTest : public ::testing::Test {};

TEST_F(DriverTest, empty_coax)
{
	// Empty Coaxial case.
	auto out{ Driver::loadFromFile(inputCase("empty_coax")).getMTLPUL() };

	auto CExpected{ EPSILON0_SI * 2 * M_PI / log(0.05 / 0.025) };

	const double rTol{ 0.005 };
	ASSERT_EQ(1, out.C.NumCols() * out.C.NumRows());
	EXPECT_LE(relError(CExpected, out.C(0, 0)), rTol);

	auto LExpected{ EPSILON0_SI * MU0_SI / CExpected };
	ASSERT_EQ(1, out.L.NumCols() * out.L.NumRows());
	EXPECT_LE(relError(LExpected, out.L(0, 0)), rTol);
}

TEST_F(DriverTest, partially_filled_coax)
{
	// Partially filled coax.
	// External radius -> r0 = 50 mm
	// Internal radius -> r1 = 25 mm 
	// Dielectric internal radius -> rI_dielectric = 25 mm
	// Dielectric external radius -> rO_dielectric = 35 mm
	// Dielectric permittivity -> eps_r = 4.0
	auto out{ Driver::loadFromFile(inputCase("partially_filled_coax")).getMTLPUL() };

	// Equivalent capacity is the series of the inner and outer capacitors.
	auto COut{ EPSILON0_SI * 2 * M_PI / log(0.050 / 0.035) };
	auto CIn{ 4.0 * EPSILON0_SI * 2 * M_PI / log(0.035 / 0.025) };
	auto CExpected = COut * CIn / (COut + CIn);

	const double rTol{ 0.002 }; // 0.2% Error.
	ASSERT_EQ(1, out.C.NumCols() * out.C.NumRows());
	EXPECT_LE(relError(CExpected, out.C(0, 0)), rTol);

	// Inductance can be calculated from free-space capacity (C0).
	auto C0 = EPSILON0_SI * 2 * M_PI / log(0.05 / 0.025);
	auto LExpected{ EPSILON0_SI * MU0_SI / C0 };

	ASSERT_EQ(1, out.L.NumCols() * out.L.NumRows());
	EXPECT_LE(relError(LExpected, out.L(0, 0)), rTol);
}

TEST_F(DriverTest, partially_filled_coax_by_domains)
{
	auto dr{ Driver::loadFromFile(inputCase("partially_filled_coax")) };
	auto globalOut{ dr.getMTLPUL() };

	auto out{ dr.getMTLPULByDomains() };
	EXPECT_EQ(1, out.domainTree.verticesSize());
	ASSERT_EQ(1, out.domainToPUL.size());


	// There are some minor differences in output. 
	// I do not know why but my guess is that it is due to initial seeds in
	// the iterative solver.
	const double rTol{ 1e-6 };
	ASSERT_EQ(1, out.domainToPUL.count(0));
	ASSERT_EQ(globalOut.L.NumRows(), out.domainToPUL.at(0).L.NumRows());
	ASSERT_EQ(globalOut.C.NumRows(), out.domainToPUL.at(0).C.NumRows());
	EXPECT_LE(relError(globalOut.L(0, 0), out.domainToPUL.at(0).L(0, 0)), rTol);
	EXPECT_LE(relError(globalOut.C(0, 0), out.domainToPUL.at(0).C(0, 0)), rTol);
}

TEST_F(DriverTest, two_wires_coax)
{
	const std::string CASE{ "two_wires_coax" };
	auto fn{ casesFolder() + CASE + "/" + CASE + ".pulmtln.in.json" };

	mfem::DenseMatrix CMatExpected(2, 2);
	CMatExpected(0, 0) = 2.15605359;
	CMatExpected(0, 1) = -0.16413431;
	CMatExpected(1, 0) = CMatExpected(0, 1);
	CMatExpected(1, 1) = CMatExpected(0, 0);
	CMatExpected *= EPSILON0_SI;

	const double rTol{ 2.5e-2 };

	auto out{ Driver::loadFromFile(fn).getMTLPUL() };

	const int N{ 2 };
	ASSERT_EQ(N, out.C.NumCols());
	ASSERT_EQ(N, out.C.NumRows());

	// Compares with analytical solution.
	for (int i{ 0 }; i < N; i++) {
		for (int j{ 0 }; j < N; j++) {
			EXPECT_LE(relError(CMatExpected(i, j), out.C(i, j)), rTol);
		}
	}

	// Checks matrix are symmetric.
	for (int i{ 0 }; i < N; i++) {
		for (int j{ 0 }; j < N; j++) {
			EXPECT_EQ(out.C(i, j), out.C(j, i));
			EXPECT_EQ(out.L(i, j), out.L(j, i));
		}
	}
}

TEST_F(DriverTest, two_wires_shielded_floating_potentials)
{
	const std::string CASE{ "two_wires_shielded" };
	auto fn{ casesFolder() + CASE + "/" + CASE + ".pulmtln.in.json" };
		
	const double rTol{ 2.5e-2 };

	auto fp{ Driver::loadFromFile(fn).getFloatingPotentials().electric };

	const int N{ 2 };
	ASSERT_EQ(N, fp.NumCols());
	ASSERT_EQ(N, fp.NumRows());

	// Solves problem and checks that charge is zero in the floating conductor.
	{
		auto m{ Mesh::LoadFromFile(casesFolder() + CASE + "/" + CASE + ".msh") };

		SolverInputs p;
		p.dirichletBoundaries = {
			{
				{1, 0.0},     // Conductor 0 bdr (GND).
				{2, fp(0,0)}, // Conductor 1 prescribed potential.
				{3, fp(0,1)}, // Conductor 2 floating potential.
			}
		};

		SolverOptions solverOpts;
		ElectrostaticSolver s{ m, p, solverOpts };
		s.Solve();

		// For debugging.
		ParaViewDataCollection pd{ outFolder() + CASE + "_floating", s.getMesh()};
		s.writeParaViewFields(pd);

		auto Q0 = s.getChargeInBoundary(1);
		auto Q1 = s.getChargeInBoundary(2);
		auto Q2 = s.getChargeInBoundary(3);
		EXPECT_NEAR(0.0, Q2, 1e-4);

	}
}

TEST_F(DriverTest, two_wires_open_floating_potentials)
{
	const std::string CASE{ "two_wires_open" };
	auto fn{ casesFolder() + CASE + "/" + CASE + ".pulmtln.in.json" };

	auto fp{ Driver::loadFromFile(fn).getFloatingPotentials().electric };

	const int N{ 2 };
	ASSERT_EQ(N, fp.NumCols());
	ASSERT_EQ(N, fp.NumRows());

	// Solves problem and checks that charge is zero in the floating conductor.
	{
		auto fn{ casesFolder() + CASE + "/" + CASE + ".msh" };
		auto m{ Mesh::LoadFromFile(fn) };

		SolverInputs p;
		p.dirichletBoundaries = {
			{
				{1, fp(0,0)}, // Conductor 0, prescribed.
				{2, fp(0,1)}, // Conductor 1, floating.
			}
		};
		p.openBoundaries = { 3 };

		ElectrostaticSolver s{ m, p };
		s.Solve();

		// For debugging.
		ParaViewDataCollection pd{ outFolder() + CASE + "_floating", s.getMesh() };
		s.writeParaViewFields(pd);

		auto Q0 = s.getChargeInBoundary(1);
		auto Q1 = s.getChargeInBoundary(2);
		auto Qb = s.getChargeInBoundary(3);

		EXPECT_NEAR(0.0, Q1, 1e-4); // Floating conductor,should not have charge.

	}
}

TEST_F(DriverTest, five_wires)
{
	// Five wires in round shield. 
	// Comparison with SACAMOS data (No Laplace).

	const std::string CASE{ "five_wires" };
	auto fn{ casesFolder() + CASE + "/" + CASE + ".pulmtln.in.json" };

	mfem::DenseMatrix couplingExpected(5, 5);
	double couplingExpectedData[25] = {
		 1.0000, -0.2574    , -0.2574    , -0.2574    , -0.2574    ,
		-0.2574,  1.0000    , -0.029017  , -3.5871E-05, -0.029017  ,
		-0.2574, -0.029017  ,  1         , -0.029017  , -3.5871E-05,
		-0.2574, -3.5871E-05, -0.029017	 ,  1         , -0.029017  ,
		-0.2574, -0.029017  , -3.5871E-05, -0.029017  ,  1.0000
	};
	couplingExpected.UseExternalData(couplingExpectedData, 5, 5);

	auto out{
		Driver::loadFromFile(fn).getMTLPUL().getCapacitiveCouplingCoefficients()
	};

	ASSERT_EQ(couplingExpected.NumRows(), out.NumRows());
	ASSERT_EQ(couplingExpected.NumCols(), out.NumCols());

	double rTol{ 0.05 };
	for (int i{ 0 }; i < couplingExpected.NumRows(); i++) {
		for (int j{ 0 }; j < couplingExpected.NumCols(); j++) {
			EXPECT_LE(std::abs(couplingExpected(i, j) - out(i, j)), rTol)
				<< "In C(" << i << ", " << j << ")";
		}
	}
}

TEST_F(DriverTest, three_wires_ribbon)
{
	// Three wires ribbon open problem. 
	// Comparison with Clayton Paul's book:  
	// Analysis of multiconductor transmision lines. 2007.
	// Sec. 5.2.3, p. 187.

	const std::string CASE{ "three_wires_ribbon" };
	auto fn{ casesFolder() + CASE + "/" + CASE + ".pulmtln.in.json" };

	double CExpectedData[4] = {
		 37.432, -18.716,
		 -18.716, 24.982
	};
	mfem::DenseMatrix CExpected(2, 2);
	CExpected.UseExternalData(CExpectedData, 2, 2);
	CExpected *= 1e-12;

	double LExpectedData[4] = { 
		0.74850, 0.50770,
		0.50770, 1.0154
	};
	mfem::DenseMatrix LExpected(2, 2);
	LExpected.UseExternalData(LExpectedData, 2, 2);
	LExpected *= 1e-6;

	auto out{ Driver::loadFromFile(fn).getMTLPUL() };

	// Tolerance is quite for this test. 
	// I guess that Paul's method is not very exact for this case.
	const double rTol{ 0.22 };

	ASSERT_EQ(CExpected.NumRows(), out.C.NumRows());
	ASSERT_EQ(CExpected.NumCols(), out.C.NumCols());
	for (int i{ 0 }; i < CExpected.NumRows(); i++) {
		for (int j{ 0 }; j < CExpected.NumCols(); j++) {
			EXPECT_LE(relError(CExpected(i, j), out.C(i, j)), rTol) <<
				"In C(" << i << ", " << j << ")";
		}
	}

	ASSERT_EQ(LExpected.NumRows(), out.L.NumRows());
	ASSERT_EQ(LExpected.NumCols(), out.L.NumCols());
	for (int i{ 0 }; i < LExpected.NumRows(); i++) {
		for (int j{ 0 }; j < LExpected.NumCols(); j++) {
			EXPECT_LE(relError(LExpected(i, j), out.L(i, j)), rTol) <<
				"In L(" << i << ", " << j << ")";
		}
	}
}

TEST_F(DriverTest, three_wires_ribbon_floating_potentials)
{
	// Three wires ribbon open problem. 
	// Comparison with Clayton Paul's book:  
	// Analysis of multiconductor transmision lines. 2007.
	// Sec. 5.2.3, p. 187.

	const std::string CASE{ "three_wires_ribbon" };
	auto fn{ casesFolder() + CASE + "/" + CASE + ".pulmtln.in.json" };
	
	auto fp = Driver::loadFromFile(fn).getFloatingPotentials().electric;

	// Solves problem and checks that charge is zero in the floating conductor.
	{
		auto fn{ casesFolder() + CASE + "/" + CASE + ".msh" };
		auto m{ Mesh::LoadFromFile(fn) };

		SolverInputs p;
		p.dirichletBoundaries = {
			{
				{1, fp(1,0)}, // Conductor 0 floating potential.
				{2, fp(1,1)}, // Conductor 1 prescribed potential.
				{3, fp(1,2)}, // Conductor 2 floating potential.
			}
		};
		p.openBoundaries = { 4 };
		p.domainPermittivities = {
			{
				{6, 3.5},
				{7, 3.5},
				{8, 3.5}
			}
		};

		SolverOptions solverOpts;
		ElectrostaticSolver s{ m, p, solverOpts };
		s.Solve();

		// For debugging.
		ParaViewDataCollection pd{ outFolder() + CASE + "_floating", s.getMesh() };
		s.writeParaViewFields(pd);

		auto Q0 = s.getChargeInBoundary(1);
		auto Q1 = s.getChargeInBoundary(2);
		auto Q2 = s.getChargeInBoundary(3);
		auto Qb = s.getChargeInBoundary(4);

		const double aTol{ 1e-3 };
		EXPECT_NEAR(0.0, Q0, aTol);
		EXPECT_NEAR(0.0, Q2, aTol);
		EXPECT_NEAR(0.0, Q0 + Q1 + Q2 + Qb, aTol);
	}
}

TEST_F(DriverTest, nested_coax)
{
	auto out{ Driver::loadFromFile(inputCase("nested_coax")).getMTLPUL() };


	auto C01{ EPSILON0_SI * 2.0 * M_PI / log(8.0 / 5.6) };
	auto C12{ EPSILON0_SI * 2.0 * M_PI / log(4.8 / 2.0) };

	double CExpectedData[4] = {
		  C01 + C12, -C12,
		 -C12,      C12
	};
	mfem::DenseMatrix CExpected(2, 2);
	CExpected.UseExternalData(CExpectedData, 2, 2);

	const double rTol{ 0.10 };

	ASSERT_EQ(CExpected.NumRows(), out.C.NumRows());
	ASSERT_EQ(CExpected.NumCols(), out.C.NumCols());
	for (int i{ 0 }; i < CExpected.NumRows(); i++) {
		for (int j{ 0 }; j < CExpected.NumCols(); j++) {
			EXPECT_LE(relError(CExpected(i, j), out.C(i, j)), rTol) <<
				"In C(" << i << ", " << j << ")";
		}
	}
}

TEST_F(DriverTest, nested_coax_by_domains)
{
	auto out{ Driver::loadFromFile(inputCase("nested_coax")).getMTLPULByDomains() };

	auto C01{ EPSILON0_SI * 2.0 * M_PI / log(8.0 / 5.6) };
	auto C12{ EPSILON0_SI * 2.0 * M_PI / log(4.8 / 2.0) };

	const double rTol{ 0.10 };

	EXPECT_EQ(2, out.domainTree.verticesSize());
	std::vector<std::pair<int, int>> domainConnections{ std::make_pair(0,1) };
	EXPECT_EQ(domainConnections, out.domainTree.getEdgesAsPairs());

	ASSERT_EQ(2, out.domainToPUL.size());

	ASSERT_EQ(1, out.domainToPUL.at(0).C.NumRows());
	EXPECT_LE(relError(C01, out.domainToPUL.at(0).C(0, 0)), rTol);

	ASSERT_EQ(1, out.domainToPUL.at(1).C.NumRows());
	EXPECT_LE(relError(C12, out.domainToPUL.at(1).C(0, 0)), rTol);
}

TEST_F(DriverTest, agrawal1981)
{
	// Agrawal, Ashok K. and Price, Harold J.
	// Experimental Characterization of Partially Degenerate Three-Conductor
	// Transmission Lines in the Time Domain. IEEE-TEMC. 1981.

	const std::string CASE{ "agrawal1981" };
	auto fn{ casesFolder() + CASE + "/" + CASE + ".pulmtln.in.json" };

	// P.U.L. Capacitances obtained from eigenvectors.
	const int nConductors = 3;
	double CExpectedData[nConductors * nConductors] = {
		  74.54, -34.57, -34.2,
		 -34.63,  73.87, -33.96,
		 -34.29, -34.0,   73.41
	};
	mfem::DenseMatrix CExpected(nConductors, nConductors);
	CExpected.UseExternalData(CExpectedData, nConductors, nConductors);
	CExpected *= 1e-12;

	auto out{ Driver::loadFromFile(fn).getMTLPUL() };

	const double rTol{ 0.10 };

	ASSERT_EQ(CExpected.NumRows(), out.C.NumRows());
	ASSERT_EQ(CExpected.NumCols(), out.C.NumCols());
	for (int i{ 0 }; i < CExpected.NumRows(); i++) {
		for (int j{ 0 }; j < CExpected.NumCols(); j++) {
			EXPECT_LE(relError(CExpected(i, j), out.C(i, j)), rTol) <<
				"In C(" << i << ", " << j << ")";
		}
	}
}

TEST_F(DriverTest, DISABLED_lansink2024_inner_multipole_boundaries_o1) // Disabling test until deciding on expected results
{
	const std::string CASE{ "lansink2024_inner" };
	auto fn{ casesFolder() + CASE + "/" + CASE + ".pulmtln.in.json" };
	
	// To comput floating potentials we are using first order  ABC.
	auto fp{ Driver::loadFromFile(fn).getFloatingPotentials().electric };

	SolverInputs p;
	p.dirichletBoundaries = { {
		{1,  fp(0,0)}, // Conductor 0 bdr.
		{2,  fp(0,1)}, // Conductor 1 bdr.
	} };

	auto m{ Mesh::LoadFromFile(casesFolder() + CASE + "/" + CASE + ".msh") };

	const int numberOfPoles = 3;
	for (int n = 0; n < numberOfPoles; n++) {

		ElectrostaticSolver s{ m, p };
		
		// Sets multipolar expansion over internal boundary.
		mfem::Vector origin({ 0.0, 0.0 });
		{
			std::vector<multipolarCoefficient> ab(n + 1);
			std::fill(ab.begin(), ab.end(), std::make_pair(0.0,0.0));
			ab.back() = { -1.0, 0.0 };
			std::function<double(const Vector&)> f =
				std::bind(&multipolarExpansion, std::placeholders::_1, ab, origin);
			FunctionCoefficient fc(f);
			s.setNeumannCondition(3, fc);
			s.Solve();

			std::stringstream ss;
			ss << n;
			ParaViewDataCollection pd(outFolder() + CASE + "_a"+ss.str(), s.getMesh());
			s.writeParaViewFields(pd);

			auto Q1{ s.getChargeInBoundary(1) };
			auto Q2{ s.getChargeInBoundary(2) };
			auto Qb{ s.getChargeInBoundary(3) };

			EXPECT_NEAR(0.0, Q1 + Q2 + Qb, 1e-3);
		}

		if (n > 0) {
			std::vector<multipolarCoefficient> ab(n + 1);
			std::fill(ab.begin(), ab.end(), std::make_pair(0.0, 0.0));
			ab.back() = { 0.0, -1.0 };
			std::function<double(const Vector&)> f =
				std::bind(&multipolarExpansion, std::placeholders::_1, ab, origin);
			FunctionCoefficient fc(f);
			s.setNeumannCondition(3, fc);
			s.Solve();

			std::stringstream ss;
			ss << n;
			ParaViewDataCollection pd(outFolder() + CASE + "_b" + ss.str(), s.getMesh());
			s.writeParaViewFields(pd);

			auto Q1{ s.getChargeInBoundary(1) };
			auto Q2{ s.getChargeInBoundary(2) };
			auto Qb{ s.getChargeInBoundary(3) };

			EXPECT_NEAR(0.0, Q1 + Q2 + Qb, 1e-3);
		}

	}
}

TEST_F(DriverTest, lansink2024_floating_potentials)
{
	// From:
	// Rotgerink, J.L. et al. (2024, September).
	// Numerical Computation of In - cell Parameters for Multiwire Formalism in FDTD.
	// In 2024 International Symposium on Electromagnetic Compatibility
	// EMC Europe(pp. 334 - 339). IEEE.
	
	const std::string CASE{ "lansink2024" };
	
	auto dr{ Driver::loadFromFile(casesFolder() + CASE + "/" + CASE + ".pulmtln.in.json") };

	auto fp{ dr.getFloatingPotentials().electric }; 
	auto inCell{ dr.getInCellPotentials() };

	auto m{ Mesh::LoadFromFile(casesFolder() + CASE + "/" + CASE + ".msh") };

	SolverInputs p;
	p.dirichletBoundaries = {
		{
			{1, fp(0,0)}, // Conductor 0 floating potential.
			{2, fp(0,1)}, // Conductor 1 prescribed potential.
		}
	};
	p.openBoundaries = { 3 };

	ElectrostaticSolver s{ m, p };
	s.Solve();

	auto Q0 = s.getChargeInBoundary(1);
	auto Q1 = s.getChargeInBoundary(2);
	auto Qb = s.getChargeInBoundary(3);

	// For debugging.
	ParaViewDataCollection pd{ outFolder() + CASE + "_floating", s.getMesh() };
	s.writeParaViewFields(pd);

	// Expectations.
	const double aTol{ 1e-3 };
	EXPECT_NEAR(0.0, Q1, aTol);
	EXPECT_NEAR(0.0, Q0 + Q1 + Qb, aTol);

	const double a0 = inCell.electric.at(0).ab[0].first;
	EXPECT_NEAR(Q0, a0, 1e-4);
}

TEST_F(DriverTest, lansink2024_fdtd_in_cell_parameters_around_conductor_1)
{
	// From:
	// Rotgerink, J.L. et al. (2024, September).
	// Numerical Computation of In - cell Parameters for Multiwire Formalism in FDTD.
	// In 2024 International Symposium on Electromagnetic Compatibility
	// EMC Europe(pp. 334 - 339). IEEE.

	const std::string CASE{ "lansink2024_fdtd_cell" };

	auto inCell{
		Driver::loadFromFile(
			casesFolder() + CASE + "/" + CASE + ".pulmtln.in.json"
		).getInCellPotentials()
	};

	const double rTol = 0.03;

	// In this test case inner region coincides with fdtd-cell.
	// In-cell capacitances.
	{
		auto computedC00 = Driver::getInCellCapacitanceUsingInnerRegion(inCell, 0, 0);
		auto expectedC00 = 14.08e-12; // C11 for floating in paper. Table 1.
		EXPECT_NEAR(0.0, relError(expectedC00, computedC00), rTol);
	}

	{
		auto computedC01 = Driver::getInCellCapacitanceUsingInnerRegion(inCell, 0, 1);
		auto expectedC01 = 43.99e-12; // C12 for floating in paper. Table 1.
		EXPECT_NEAR(0.0, relError(expectedC01, computedC01), rTol);
	}

	// In-cell inductances
	{
		auto computedL00 = Driver::getInCellInductanceUsingInnerRegion(inCell, 0, 0);
		auto expectedL00 = 791e-9; // L11 for floating in paper. Table 1.
		EXPECT_NEAR(0.0, relError(expectedL00, computedL00), rTol);
	}

	{
		auto computedL01 = Driver::getInCellInductanceUsingInnerRegion(inCell, 0, 1);
		auto expectedL01 = 253e-9; // L12 for floating in paper. Table 1.
		EXPECT_NEAR(0.0, relError(expectedL01, computedL01), rTol);
	}
}

TEST_F(DriverTest, lansink2024_single_wire_in_cell_parameters)
{
	// From:
	// Rotgerink, J.L. et al. (2024, September).
	// Numerical Computation of In - cell Parameters for Multiwire Formalism in FDTD.
	// In 2024 International Symposium on Electromagnetic Compatibility
	// EMC Europe(pp. 334 - 339). IEEE.
	// VALUES IN TABLE 3 ARE WRONG AND HAVE BEEN CORRECTED.

	const std::string CASE{ "lansink2024_single_wire" };

	auto inCell{
		Driver::loadFromFile(
			casesFolder() + CASE + "/" + CASE + ".pulmtln.in.json"
		).getInCellPotentials()
	};

	const double rTol = 0.05;

	// In this test case inner region coincides with fdtd-cell.
	// In-cell capacitances.
	{
		auto computedC00 = Driver::getInCellCapacitanceUsingInnerRegion(inCell, 0, 0);
		auto expectedC00 = 49.11e-12; // C11 with insulation. Table 3. 
		                              // Paper has a mistake, this is the correct value.
		EXPECT_NEAR(0.0, relError(expectedC00, computedC00), rTol);
	}

	// In-cell inductances
	{
		auto computedL00 = Driver::getInCellInductanceUsingInnerRegion(inCell, 0, 0);
		auto expectedL00 = 320e-9; // L11 with insulation. Table 3. 
								   // Paper has a mistake, this is the correct value.
		EXPECT_NEAR(0.0, relError(expectedL00, computedL00), rTol);
	}

}