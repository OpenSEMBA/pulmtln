#pragma once

#include "SolverOptions.h"
#include "Model.h"
#include "Parameters.h"

namespace pulmtln {

class Driver {
public:
    Driver(Model&& model, const DriverOptions& opts);
    
    PULParameters getMTLPUL() const;
    PULParametersByDomain getMTLPULByDomains() const;

    static Driver loadFromFile(const std::string& filename);

private:
    Model model_;
    DriverOptions opts_;
};

}