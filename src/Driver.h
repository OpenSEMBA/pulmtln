#pragma once

#include "SolverOptions.h"
#include "Model.h"
#include "Parameters.h"

namespace pulmtln {

struct Domain {
    using Id = int;
    using ConductorId = int;

    ConductorId ground;
    std::vector<ConductorId> conductorIds;
    std::set<int> elements;
};

class Driver {
public:


    Driver(const Model& model, const DriverOptions& opts);
    
    PULParameters getMTLPUL() const;
    PULParametersByDomain getMTLPULByDomain() const;

    static Driver loadFromFile(const std::string& filename);

private:
    Model model_;
    DriverOptions opts_;

    std::map<Domain::Id, Domain> domains_;

};

}