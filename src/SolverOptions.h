#pragma once

namespace pulmtln {

struct SolverOptions {
	int order{3};     // Basis function order
	bool exportParaViewSolution{ true };
};

}