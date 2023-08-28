#include <gtest/gtest.h>

#include "TestUtils.h"

#include "constants.h"
#include "Driver.h"

using namespace pulmtln;

using json = nlohmann::json;

class DriverTest : public ::testing::Test {};

TEST_F(DriverTest, empty_coax)
{
	// Empty Coaxial case.
	const std::string CASE{ "empty_coax" };
	auto fn{ casesFolder() + CASE + "/" + CASE + ".pulmtln.in.json" };
	
	auto out{ Driver::loadFromFile(fn).getMTLPUL() };

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

	const std::string CASE{ "partially_filled_coax" };
	auto fn{ casesFolder() + CASE + "/" + CASE + ".pulmtln.in.json" };

	auto out{ Driver::loadFromFile(fn).getMTLPUL() };

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
	
	double rTol{ 0.15 };
	for (int i{ 0 }; i < 5; i++) {
		for (int j{ 0 }; j < 5; j++) {
			EXPECT_LE(relError(couplingExpected(i, j), out(i, j)), rTol) << 
				"In C(" << i << ", " << j << ")";
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
		 37.8189, -18.0249,
		-18.0249,  26.2148
	};
	mfem::DenseMatrix CExpected(2, 2);
	CExpected.UseExternalData(CExpectedData, 2, 2);
	CExpected *= 1e-12;

	auto out{ 
		Driver::loadFromFile(fn).getMTLPUL().C 
	};

	ASSERT_EQ(CExpected.NumRows(), out.NumRows());
	ASSERT_EQ(CExpected.NumCols(), out.NumCols());
	
	double rTol{ 0.05 };
	for (int i{ 0 }; i < CExpected.NumRows(); i++) {
		for (int j{ 0 }; j < CExpected.NumCols(); j++) {
			EXPECT_LE(relError(CExpected(i, j), out(i, j)), rTol) << 
				"In C(" << i << ", " << j << ")";
		}
	}
}
