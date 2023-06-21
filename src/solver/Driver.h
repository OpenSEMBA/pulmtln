#pragma once

#include "SolverOptions.h"
#include "Model.h"

namespace pulmtln {

class Driver {
public:
    Driver(const Model&, const SolverOptions& opts);
private:
    Model model_;
    SolverOptions opts_;
};

}