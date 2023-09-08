#pragma once

#include "SolverOptions.h"
#include "Model.h"
#include "Parameters.h"

namespace pulmtln {

class Driver {
public:
    Driver(const Model& model, const SolverOptions& opts);
    
    Parameters getMTLPUL() const;

    static Driver loadFromFile(const std::string& filename);

private:
    Model model_;
    SolverOptions opts_;
};

}