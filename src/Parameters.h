#pragma once

#include <mfem.hpp>
#include <nlohmann/json.hpp>

#include "Domain.h"

namespace pulmtln {

struct PULParameters {
    PULParameters() = default;
    PULParameters(const nlohmann::json&);

    bool operator==(const PULParameters&) const;

    mfem::DenseMatrix getCapacitiveCouplingCoefficients() const;

    nlohmann::json toJSON() const;

    void saveToJSONFile(const std::string& filename) const;
    
    mfem::DenseMatrix L, C; // Stored in SI units.
};

struct PULParametersByDomain {
    DomainTree domainTree;
    std::map<Domain::Id, PULParameters> domainToPUL;
};

}