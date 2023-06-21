#pragma once

#include <stdexcept>
#include <vector>

#include "core/math/CartesianVector.h"
#include "core/math/function/Polynomial.h"

namespace semba::math::simplex {

class Simplex {
public:
    Simplex();
    virtual ~Simplex();
    virtual const function::Polynomial<Real>& getLagr(
        const std::size_t i) const = 0;
    virtual const function::Polynomial<Real>& getDLagr(
        const std::size_t i,
        const std::size_t f) const = 0;
    virtual std::vector<Real> getWeights() const = 0;

protected:

    virtual size_t nodeIndex(const std::size_t node,
                             const std::size_t coordinate) const = 0;

    function::Polynomial<Real> silvesterPol(const std::size_t m,
                                            const std::size_t n) const;
    void lagrangePolynomials(function::Polynomial<Real>* lagr,
                             const std::size_t n,
                             const std::size_t np,
                             const std::size_t nsc) const;

    Real integrate(const function::Polynomial<Real> pol,
                   const std::size_t dimension,
                   const Real sizeFactor) const;
    static std::size_t factorial(std::size_t n);
};

}
