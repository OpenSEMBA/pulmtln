#pragma once

#include "core/math/matrix/Static.h"

#include "Simplex.h"

#include <algorithm>

namespace semba::math::simplex {

template <size_t N>
class Line : public Simplex {
public:
    static const std::size_t faces = 2;
    static const std::size_t dimension = 1;
    static const std::size_t nsc = 2;
    static const std::size_t nfp = 1;
    static constexpr std::size_t np = (N + 1);
    static constexpr Real sizeFactor = 1.0;

    typedef CartesianVector<size_t,nsc> Index;

    Line();
    std::size_t vertex(const std::size_t) const;
    std::size_t sideVertex(const std::size_t f, const std::size_t i) const;
    std::size_t sideNode( const std::size_t f, const std::size_t i) const;

    std::size_t nodeIndex(const std::size_t i, const std::size_t j) const;

    const function::Polynomial<Real>& getLagr(
            const std::size_t node) const;
    const function::Polynomial<Real>& getDLagr(
            const std::size_t node, const std::size_t simplex) const;

    std::vector<Real> getWeights() const;

private:
    std::array<Index,np> indices;
    matrix::Static<Int,faces,nfp> sideNodes;

    function::Polynomial<Real> lagr[np];
    function::Polynomial<Real> dLagr[np][faces];

    std::array<CartesianVector<Real,nsc>, np> nodePositions;
    std::array<Real,np>                         weights;

    matrix::Static<Int, 1, (N+1)> RMatrix(const std::size_t s) const;
    matrix::Static<Int, (N+1), (N+1)> PMatrix(const std::size_t s) const;
};


template <size_t N>
Line<N>::Line() {

    for (std::size_t i = 0; i < indices.size(); i++) {
        indices[i](0) = N - i;
        indices[i](1) = i;
    }

    matrix::Static<Int, np, 1> nList;
    for (std::size_t i = 0; i < np; i++) {
        nList(i, 0) = i;
    }

    for (std::size_t f = 0; f < faces; f++) {
        matrix::Static<Int, nfp, 1> aux = RMatrix(f) * nList;
        for (std::size_t i = 0; i < nfp; i++) {
            sideNodes(f, i) = aux(i, 0);
        }
    }

    lagrangePolynomials(lagr, N, np, nsc);
    for (std::size_t i = 0; i < np; i++) {
        for (std::size_t s = 0; s < nsc; s++) {
            dLagr[i][s] = lagr[i];
            dLagr[i][s].derive(s);
        }
    }

    CartesianVector<Real, nsc> aux;
    for (std::size_t i = 0; i < np; i++) {
        aux = indices[i];
        nodePositions[i] = aux / (Real)N;
    }
    for (std::size_t i = 0; i < np; i++) {
        weights[i] = integrate(lagr[i], dimension, sizeFactor) / sizeFactor;
    }
};

template <size_t N>
inline std::size_t Line<N>::nodeIndex(const std::size_t i,
    const std::size_t j) const {
    return indices[i](j);
}

template <size_t N>
inline std::size_t Line<N>::vertex(const std::size_t vertexNum) const {
    return sideNode(vertexNum, 0);
}

template <size_t N>
inline std::size_t Line<N>::sideNode(const std::size_t face,
    const std::size_t num) const {
    return sideNodes(face, num);
}

template <size_t N>
inline const function::Polynomial<Real>& Line<N>::getLagr(
    const std::size_t i) const {
    return lagr[i];
}

template <size_t N>
inline const function::Polynomial<Real>& Line<N>::getDLagr(
    const std::size_t i,
    const std::size_t f) const {
    return dLagr[i][f];
}

template <size_t N>
inline std::vector<Real> Line<N>::getWeights() const {
    std::vector<Real> res(np);
    std::copy_n(weights.begin(), np, res.begin());
    return res;
}

template <std::size_t N>
matrix::Static<Int, (N + 1), (N + 1)> Line<N>::PMatrix(const std::size_t s) const {
    matrix::Static<Int, np, np> res;
    if (s == 0) {
        res.eye();
    }
    else {
        res.zeros();
        for (std::size_t i = 0; i < np; i++) {
            res(i, np - i - 1) = (Int)1;
        }
    }
    return res;
}

template <size_t N>
matrix::Static<Int, 1, (N + 1)> Line<N>::RMatrix(const std::size_t s) const {
    matrix::Static<Int, nfp, np> Raux;
    Raux.zeros();
    Raux(0, 0) = (Int)1;
    matrix::Static<Int, np, np> P;
    P = PMatrix(s);
    matrix::Static<Int, nfp, np> res = Raux * P;
    return res;
}

} 


