#pragma once

namespace pulmtln {

struct SolverOptions {
	int order = 1;     // Basis function order
    int maxit = 100;
    int serial_ref_levels = 0;
    int parallel_ref_levels = 0;
    bool visualization = true;
};

}