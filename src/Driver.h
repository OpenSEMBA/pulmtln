#pragma once

#include "SolverOptions.h"
#include "Model.h"

#include <nlohmann/json.hpp>

namespace pulmtln {

class Driver {
public:
    Driver(const nlohmann::json& input);

    nlohmann::json getPULMTLN() const;
private:
    Model model_;
    SolverOptions opts_;
};

}