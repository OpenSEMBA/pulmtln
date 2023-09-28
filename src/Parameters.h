#pragma once

#include <mfem.hpp>
#include <nlohmann/json.hpp>

#include "DirectedGraph.h"

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
    using DomainId = VertexId;

    DirectedGraph domainGraph;
    std::map<DomainId, PULParameters> domainToPUL;
};

}