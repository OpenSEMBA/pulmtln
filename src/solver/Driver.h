#pragma once

#include "SolverOptions.h"
#include "Model.h"
#include "ElectrostaticSolver.h"

namespace pulmtln {

class Driver {
public:
    Driver(const Model&, const SolverOptions& opts);

    const ElectrostaticSolver& getElectrostaticSolver() const;
private:
    Model model_;
    SolverOptions opts_;
    std::unique_ptr<ElectrostaticSolver> electrostaticSolver_;
};

}