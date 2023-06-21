#include <gtest/gtest.h>

#include "solver/Solver.h"


using namespace mfem;
using namespace pulmtln;

class pulmtlnTest : public ::testing::Test {
};

TEST_F(pulmtlnTest, parallel_plates)
{
	auto m{
		Mesh::MakeCartesian2D
	};
	SolverOptions opts;
	opts.order = 2;

	pulmtln::Solver solver {
		m,
		opts,
	};

}
