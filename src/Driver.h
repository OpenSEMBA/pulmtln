#pragma once

#include "SolverOptions.h"
#include "Model.h"

#include <nlohmann/json.hpp>

namespace pulmtln {

struct MTLPULParameters {
    mfem::DenseMatrix L, C; // Stored in SI units.
};

class Driver {
public:
    Driver(const Model& model, const SolverOptions& opts);
    
    MTLPULParameters getMTLPUL() const;

    static Driver loadFromFile(const std::string& filename);

private:
    Model model_;
    SolverOptions opts_;
};

}