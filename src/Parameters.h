#pragma once

#include <nlohmann/json.hpp>

#include "FES.h"
#include "Domain.h"
#include "multipolarExpansion.h"

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

struct FloatingPotentials {
    mfem::DenseMatrix electric, magnetic;
};

struct InCellParameters {
    struct FieldParameters {
        double innerRegionAveragePotential;
        mfem::Vector expansionCenter;
        multipolarCoefficients ab;
    };
    double innerRegionRadius;

    // Electric and magnetic potentials multipolar expansions for each conductor.
    std::map<std::string, FieldParameters> electric, magnetic;
};

}