#include <gtest/gtest.h>

#include "solver/Driver.h"


using namespace mfem;
using namespace pulmtln;

class pulmtlnTest : public ::testing::Test {
};

TEST_F(pulmtlnTest, parallel_plates)
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
	std::map<int, double> dirichletBCs{
		{1,    1.0}, // bottom boundary.
		{3,    0.0}, // top boundary.
	};
	Model model{ mesh, dirichletBCs };

	SolverOptions opts;
	opts.order = 3;

	pulmtln::Driver driver{ model, opts	};
	
	const double tol{ 1e-5 };

	auto totalChargeRho{ 
		driver.getElectrostaticSolver().computeTotalChargeFromRho() 
	};
	EXPECT_NEAR(0.0, totalChargeRho, tol);
	auto totalCharge{
		driver.getElectrostaticSolver().computeTotalCharge()
	};
	EXPECT_NEAR(0.0, totalCharge, tol);

	mfem::Array<int> attr(1);

	attr[0] = 1;
	auto chargeAtBottom{ 
		driver.getElectrostaticSolver().computeChargeInBoundary(attr)
	};
	EXPECT_NEAR( 1.0, chargeAtBottom, tol);
	
	attr[0] = 3;
	auto chargeAtTop{
		driver.getElectrostaticSolver().computeChargeInBoundary(attr)
	};
	EXPECT_NEAR(-1.0, chargeAtTop, tol);
	
	attr[0] = 2;
	auto chargeAtRight{
		driver.getElectrostaticSolver().computeChargeInBoundary(attr)
	};
	EXPECT_NEAR( 0.0, chargeAtRight, tol);
}

TEST_F(pulmtlnTest, surface_integral)
{
	auto mesh{
		Mesh::MakeCartesian2D(10, 10, Element::QUADRILATERAL, 1.0, 1.0)
	};
	RT_FESpace hdiv(&mesh, 4, mesh.Dimension());
	
	GridFunction d(&hdiv);
	VectorConstantCoefficient xU(Vector({1.0, 0.0}));
	d.ProjectCoefficient(xU);
	
	{
		LinearForm surf_int(&hdiv);
		surf_int.AddBoundaryIntegrator(new VectorFEBoundaryFluxLFIntegrator);
		surf_int.Assemble();
		
		auto res{ surf_int(d) };

		EXPECT_NEAR(0.0, res, 1e-12);
	}

	{
		LinearForm left_bdr_int(&hdiv);
		
		Array<int> attr(1);
		attr[0] = 4;

		Array<Coefficient*> coefs(1);
		ConstantCoefficient one{ 1.0 };
		coefs[0] = &one;

		PWCoefficient pwcoeff{ attr, coefs };

		left_bdr_int.AddBoundaryIntegrator(
			new VectorFEBoundaryFluxLFIntegrator(pwcoeff));
		left_bdr_int.Assemble();
		
		auto res{ left_bdr_int(d) };

		EXPECT_NEAR(-1.0, res, 1e-12);
	}

	{
		LinearForm right_bdr_int(&hdiv);

		Array<int> attr(1);
		attr[0] = 2;

		Array<Coefficient*> coefs(1);
		ConstantCoefficient one{ 1.0 };
		coefs[0] = &one;

		PWCoefficient pwcoeff{ attr, coefs };

		right_bdr_int.AddBoundaryIntegrator(
			new VectorFEBoundaryFluxLFIntegrator(pwcoeff));
		right_bdr_int.Assemble();

		auto res{ right_bdr_int(d) };

		EXPECT_NEAR(1.0, res, 1e-12);
	}
}