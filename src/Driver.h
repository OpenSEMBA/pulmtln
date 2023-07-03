#pragma once

#include "SolverOptions.h"
#include "Model.h"

#include <nlohmann/json.hpp>

namespace pulmtln {

struct MTLPULParameters {
    mfem::DenseMatrix L, C; // Stored in SI.
};

class Driver {
public:
    Driver(const nlohmann::json& input);

    MTLPULParameters getMTLPUL() const;
private:
    Model model_;
    SolverOptions opts_;
};

}