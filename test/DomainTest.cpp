#include "gtest/gtest.h"

#include "TestUtils.h"

#include "Domain.h"
#include "Parser.h"

using namespace pulmtln;

class DomainTest : public ::testing::Test {};

TEST_F(DomainTest, build_domains_for_empty_coax)
{
	auto model{ Parser{ inputCase("empty_coax") }.readModel() };
	
	auto domains{ Domain::buildDomains(model) };

	ASSERT_EQ(1, domains.size());
	EXPECT_EQ(1, domains.count(0));

	EXPECT_EQ(0, domains.at[0].ground);

	ASSERT_EQ(1, domains.at[0].conductorIds.size());
	EXPECT_EQ(1, domains.at[0].conductorIds[0]);

	EXPECT_EQ(model.getMesh()->GetNE(), domains.at[0].elements.size());
}


//TEST_F(DriverTest, five_wires) // TODO
//{
//	const std::string CASE{ "five_wires" };
//	auto fn{ casesFolder() + CASE + "/" + CASE + ".pulmtln.in.json" };
//
//}

//TEST_F(DriverTest, three_wires_ribbon) // TODO
//{
//	// Three wires ribbon open problem. 
//	// Comparison with Paul's book: 
//	// Analysis of multiconductor transmision lines. 2007.
//	// Sec. 5.2.3, p. 187.
//
//	const std::string CASE{ "three_wires_ribbon" };
//	auto fn{ casesFolder() + CASE + "/" + CASE + ".pulmtln.in.json" };
//
//	double CExpectedData[4] = {
//		 37.432, -18.716,
//		-18.716,  24.982
//	};
//	mfem::DenseMatrix CExpected(2, 2);
//	CExpected.UseExternalData(CExpectedData, 2, 2);
//	CExpected *= 1e-12;
//
//	double LExpectedData[4] = {
//		0.74850, 0.50770,
//		0.50770, 1.0154
//	};
//	mfem::DenseMatrix LExpected(2, 2);
//	LExpected.UseExternalData(LExpectedData, 2, 2);
//	LExpected *= 1e-6;
//
//	auto out{ Driver::loadFromFile(fn).getMTLPUL() };
//	
//	// Tolerance is quite high probably because open region is not far enough.
//	const double rTol{ 0.21 }; 
//	
//	ASSERT_EQ(CExpected.NumRows(), out.C.NumRows());
//	ASSERT_EQ(CExpected.NumCols(), out.C.NumCols());
//	for (int i{ 0 }; i < CExpected.NumRows(); i++) {
//		for (int j{ 0 }; j < CExpected.NumCols(); j++) {
//			EXPECT_LE(relError(CExpected(i, j), out.C(i, j)), rTol) << 
//				"In C(" << i << ", " << j << ")";
//		}
//	}
//
//	ASSERT_EQ(LExpected.NumRows(), out.L.NumRows());
//	ASSERT_EQ(LExpected.NumCols(), out.L.NumCols());
//	for (int i{ 0 }; i < LExpected.NumRows(); i++) {
//		for (int j{ 0 }; j < LExpected.NumCols(); j++) {
//			EXPECT_LE(relError(LExpected(i, j), out.L(i, j)), rTol) <<
//				"In L(" << i << ", " << j << ")";
//		}
//	}
//}

//TEST_F(DriverTest, nested_coax) // TODO
//{
//	// Coaxial inside a coaxial
//
//	const std::string CASE{ "nested_coax" };
//	auto fn{ casesFolder() + CASE + "/" + CASE + ".pulmtln.in.json" };
//
//	auto C01{ EPSILON0_SI * 2.0 * M_PI / log(8.0 / 5.6) };
//	auto C12{ EPSILON0_SI * 2.0 * M_PI / log(4.8 / 2.0) };
//
//	double CExpectedData[4] = {
//		  C01+C12, -C12,
//		 -C12,      C12
//	};
//	mfem::DenseMatrix CExpected(2, 2);
//	CExpected.UseExternalData(CExpectedData, 2, 2);
//
//	auto out{ Driver::loadFromFile(fn).getMTLPUL() };
//
//	const double rTol{ 0.10 };
//
//	ASSERT_EQ(CExpected.NumRows(), out.C.NumRows());
//	ASSERT_EQ(CExpected.NumCols(), out.C.NumCols());
//	for (int i{ 0 }; i < CExpected.NumRows(); i++) {
//		for (int j{ 0 }; j < CExpected.NumCols(); j++) {
//			EXPECT_LE(relError(CExpected(i, j), out.C(i, j)), rTol) <<
//				"In C(" << i << ", " << j << ")";
//		}
//	}
//}

//TEST_F(DriverTest, agrawal1981) // TODO
//{
//	// Agrawal, Ashok K. and Price, Harold J.
//	// Experimental Characterization of Partially Degenerate Three-Conductor
//	// Transmission Lines in the Time Domain. IEEE-TEMC. 1981.
//
//	const std::string CASE{ "agrawal1981" };
//	auto fn{ casesFolder() + CASE + "/" + CASE + ".pulmtln.in.json" };
//
//	// P.U.L. Capacitances obtained from eigenvectors.
//	const int nConductors = 3;
//	double CExpectedData[nConductors*nConductors] = {
//		  74.54, -34.57, -34.2,
//		 -34.63,  73.87, -33.96,
//		 -34.29, -34.0,   73.41
//	};
//	mfem::DenseMatrix CExpected(nConductors, nConductors);
//	CExpected.UseExternalData(CExpectedData, nConductors, nConductors);
//	CExpected *= 1e-12;
//
//	auto out{ Driver::loadFromFile(fn).getMTLPUL() };
//
//	const double rTol{ 0.10 };
//
//	ASSERT_EQ(CExpected.NumRows(), out.C.NumRows());
//	ASSERT_EQ(CExpected.NumCols(), out.C.NumCols());
//	for (int i{ 0 }; i < CExpected.NumRows(); i++) {
//		for (int j{ 0 }; j < CExpected.NumCols(); j++) {
//			EXPECT_LE(relError(CExpected(i, j), out.C(i, j)), rTol) <<
//				"In C(" << i << ", " << j << ")";
//		}
//	}
//}
