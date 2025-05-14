#pragma once

#include "FES.h"

namespace pulmtln {

using multipolarCoefficient = std::pair<double, double>;
using multipolarCoefficients = std::vector<multipolarCoefficient>;

static double multipolarExpansion(
    const mfem::Vector& position,
    const multipolarCoefficients& ab,
    const mfem::Vector expansionCenter)
{
    // 2D multipolar expansion from:
    // TSOGTGEREL GANTUMUR, MULTIPOLE EXPANSIONS IN THE PLANE. 
    // 2016-04-16 lecture notes. 

    mfem::Vector rVec{ position };
    rVec -= expansionCenter;

    double r{ rVec.Norml2() };
    double phi{ std::atan(rVec(1) / rVec(0)) };

    double res{ 0.0 };
    for (int n{ 0 }; n < ab.size(); ++n) {
        const auto& an = ab[n].first;
        const auto& bn = ab[n].second;
        if (n == 0) {
            res -= an * std::log(r);
            assert(bn == 0.0); // b0 should always be zero.
        }
        else {
            res += (an * std::cos(n * phi) + bn * std::sin(n * phi)) / std::pow(r, n);
        }
    }

    return res;
}

}