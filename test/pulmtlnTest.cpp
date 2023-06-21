#include <gtest/gtest.h>

#include "solver/Driver.h"


using namespace mfem;
using namespace pulmtln;

class pulmtlnTest : public ::testing::Test {
};

TEST_F(pulmtlnTest, parallel_plates)
{
	auto mesh{ 
		Mesh::MakeCartesian2D(10, 10, Element::QUADRILATERAL, 1.0, 1.0)
	};
	std::map<int, double> dirichletBCs{
		{1, 1.0}, // bottom boundary.
		{3, 0.0}, // top boundary.
	};
	std::map<int, double> neumannBCs{
		{2, 1.0}, // right boundary.
		{4, 0.0}, // left boundary.
	};
	Model model{ mesh, dirichletBCs, neumannBCs };

	SolverOptions opts;
	opts.order = 2;

	pulmtln::Driver driver{ model, opts	};

}
