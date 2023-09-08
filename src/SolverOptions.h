#pragma once

#include <string>

namespace pulmtln {

struct SolverOptions {
	int order{3};     // Basis function order
	
	bool printIterations{ false };

	bool exportParaViewSolution{ true };
	bool exportVisItSolution{ false };
	bool exportMatrices{ true };
	std::string exportFolder{ "" };
};

}