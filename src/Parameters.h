#pragma once

#include <mfem.hpp>
#include <nlohmann/json.hpp>

namespace pulmtln {

struct Parameters {
    Parameters() = default;
    Parameters(const nlohmann::json&);

    bool operator==(const Parameters&) const;

    mfem::DenseMatrix getCapacitiveCouplingCoefficients() const;

    nlohmann::json toJSON() const;

    void saveToJSONFile(const std::string& filename) const;
    
    mfem::DenseMatrix L, C; // Stored in SI units.
};

}