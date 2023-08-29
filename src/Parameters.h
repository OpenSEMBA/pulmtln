#pragma once

#include <mfem.hpp>
#include <nlohmann/json.hpp>

namespace pulmtln {

struct Parameters {
    mfem::DenseMatrix L, C; // Stored in SI units.

    mfem::DenseMatrix getCapacitiveCouplingCoefficients() const;

    nlohmann::json toJSON() const;
    void saveToJSONFile(const std::string& filename) const;
};

}