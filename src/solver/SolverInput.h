#pragma once

#include "components/Problem.h"
#include "SolverOptions.h"

namespace pulmtln {

struct SolverInput {
    Problem problem;
    SolverOptions options;
};

}
