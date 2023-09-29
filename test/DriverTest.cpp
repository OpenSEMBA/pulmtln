#include "gtest/gtest.h"

#include "TestUtils.h"

#include "constants.h"
#include "Driver.h"

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
	auto COut{       EPSILON0_SI * 2 * M_PI / log(0.050 / 0.035) };
	auto CIn{  4.0 * EPSILON0_SI * 2 * M_PI / log(0.035 / 0.025) };
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
	ASSERT_EQ(2, out.C.NumCols());
	ASSERT_EQ(2, out.C.NumRows());
	for (int i{ 0 }; i < 2; i++) {
		for (int j{ 0 }; j < 2; j++) {
			EXPECT_LE(relError(CMatExpected(i, j), out.C(i, j)), rTol);
		}
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
	// Comparison with Paul's book: 
	// Analysis of multiconductor transmision lines. 2007.
	// Sec. 5.2.3, p. 187.

	const std::string CASE{ "three_wires_ribbon" };
	auto fn{ casesFolder() + CASE + "/" + CASE + ".pulmtln.in.json" };

	double CExpectedData[4] = {
		 37.432, -18.716,
		-18.716,  24.982
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
	
	// Tolerance is quite high probably because open region is not far enough.
	const double rTol{ 0.21 }; 
	
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

TEST_F(DriverTest, nested_coax)
{
	// Coaxial inside a coaxial

	const std::string CASE{ "nested_coax" };
	auto fn{ casesFolder() + CASE + "/" + CASE + ".pulmtln.in.json" };

	auto C01{ EPSILON0_SI * 2.0 * M_PI / log(8.0 / 5.6) };
	auto C12{ EPSILON0_SI * 2.0 * M_PI / log(4.8 / 2.0) };

	double CExpectedData[4] = {
		  C01+C12, -C12,
		 -C12,      C12
	};
	mfem::DenseMatrix CExpected(2, 2);
	CExpected.UseExternalData(CExpectedData, 2, 2);

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

TEST_F(DriverTest, agrawal1981)
{
	// Agrawal, Ashok K. and Price, Harold J.
	// Experimental Characterization of Partially Degenerate Three-Conductor
	// Transmission Lines in the Time Domain. IEEE-TEMC. 1981.

	const std::string CASE{ "agrawal1981" };
	auto fn{ casesFolder() + CASE + "/" + CASE + ".pulmtln.in.json" };

	// P.U.L. Capacitances obtained from eigenvectors.
	const int nConductors = 3;
	double CExpectedData[nConductors*nConductors] = {
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
