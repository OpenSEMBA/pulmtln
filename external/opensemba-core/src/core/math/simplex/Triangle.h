#pragma once

#include <algorithm>

#include "core/math/CartesianVector.h"
#include "core/math/matrix/Static.h"
#include "core/math/function/Polynomial.h"

#include "Line.h"


namespace semba::math::simplex {

template <size_t N>
class Triangle : public Simplex {
public:
    static const std::size_t faces = 3;
    static const std::size_t dimension = 2;
    static const std::size_t vertices = 3;
    static const std::size_t np  = ((N + 1) * (N + 2) / 2);
    static const std::size_t nfp = Line<N>::np;
    static constexpr Real sizeFactor{ 1.0 / 2.0 };

    using MatIntNpNp = matrix::Static<Int, np, np>;
    using MatIntNfpNp = matrix::Static<Int, nfp, np>;

    static const std::size_t nsc = 3;
    using Index = CartesianVector<size_t, nsc>;

    Triangle();

    std::size_t vertex(const std::size_t i) const;
    std::size_t sideVertex(const std::size_t f, const std::size_t i) const;
    std::size_t sideNode  (const std::size_t f, const std::size_t i) const;

    std::size_t nodeIndex(const std::size_t i, const std::size_t j) const;

    math::CVecR3 coordinate(const std::size_t i) const;

    const function::Polynomial<Real>& getLagr(const std::size_t i) const;
    const function::Polynomial<Real>& getDLagr(const std::size_t i,
                                               const std::size_t f) const;
    std::vector<Real> getWeights() const;

    static matrix::Dynamic<Int> PMatrix(const std::size_t n,
                                        const std::size_t s);

private:
    Index indices[np];
    matrix::Static<Int,faces,nfp> sideNodes;

    function::Polynomial<Real> lagr[np];
    function::Polynomial<Real> dLagr[np][faces];

    CartesianVector<Real,nsc> nodePositions[np];
    std::array<Real,np>         weights;

    MatIntNfpNp RMatrix(const std::size_t s) const;

    static size_t numberOfNodes(size_t order);
};


template <size_t N>
Triangle<N>::Triangle() {
    matrix::Dynamic<Int> ini(np, nsc);
    for (std::size_t i = 0; i <= N; i++) {
        for (std::size_t j = numberOfNodes(std::size_t(i - 1)); j < np; j++) {
            ini(j, 0) = N - i;
        }
    }

    matrix::Dynamic<Int> ord(np, nsc);
    for (std::size_t i = 0; i < nsc; i++) {
        ord = PMatrix(N, i) * ini;
        for (std::size_t j = 0; j < np; j++) {
            indices[j](i) = ord(j, 0);
        }
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

    for (std::size_t i = 0; i < np; i++) {
        CartesianVector<Real, nsc> aux = indices[i];
        nodePositions[i] = aux / (Real)N;
    }
    for (std::size_t i = 0; i < np; i++) {
        weights[i] = integrate(lagr[i], dimension, sizeFactor) / sizeFactor;
    }
}

template <size_t N>
inline std::vector<Real> Triangle<N>::getWeights() const {
    std::vector<Real> res(np);
    std::copy_n(weights.begin(), np, res.begin());
    return res;
}

template <size_t N>
inline std::size_t Triangle<N>::nodeIndex(const std::size_t i,
    const std::size_t j) const {
    return indices[i](j);
}

template <size_t N>
const function::Polynomial<Real>& Triangle<N>::getLagr(
    const std::size_t i) const {
    return lagr[i];
}

template <size_t N>
const function::Polynomial<Real>& Triangle<N>::getDLagr(
    const std::size_t i,
    const std::size_t f) const {
    return dLagr[i][f];
}

template<size_t N>
inline std::size_t Triangle<N>::vertex(std::size_t i) const {
    return sideNode(i, 0);
}

template<size_t N>
inline std::size_t Triangle<N>::sideVertex(const std::size_t f,
    const std::size_t i) const {
    return sideNode(f, i);
}

template <size_t N>
inline std::size_t Triangle<N>::sideNode(
    const std::size_t face, const std::size_t num) const {
    return sideNodes(face, num);
}

template <size_t N>
CartesianVector<Real, 3> Triangle<N>::coordinate(
    const std::size_t i) const {
    CartesianVector<Real, 3> res;
    res = indices[i];
    res /= (Real)N;
    return res;
}

template <size_t N>
std::size_t Triangle<N>::numberOfNodes(const std::size_t order) {
    size_t res = 1;
    for (std::size_t i = 1; i < nsc; i++) {
        res *= (order + i);
    }
    res /= factorial(nsc - 1);
    return res;
}

template <size_t N>
matrix::Dynamic<Int> Triangle<N>::PMatrix(std::size_t order,
    std::size_t s) {
    std::size_t np = numberOfNodes(order);
    std::size_t nfp = order + 1;
    matrix::Dynamic<Int> res(np, np);
    matrix::Dynamic<Int> original(nfp, nfp), rotated(nfp, nfp);
    std::size_t originalNum = 1;
    for (std::size_t i = 0; i < nfp; i++) {
        for (std::size_t j = 0; j <= i; j++) {
            original(i, j) = originalNum;
            originalNum++;
        }
    }

    std::size_t rotatedNum = 1;
    if (s == 0) {
        res.eye();
    }
    else if (s == 1 || s == 2) {
        for (std::size_t i = 0; i <= nfp; i++) {
            Int j = Int(nfp - 1);
            while (j >= Int(nfp - i)) {
                rotated(j, j - nfp + i) = rotatedNum;
                rotatedNum = rotatedNum + 1;
                j--;
            }
        }
        for (std::size_t i = 0; i < nfp; i++) {
            for (std::size_t j = 0; j <= i; j++) {
                res(rotated(i, j) - 1, original(i, j) - 1) = 1;
            }
        }
        if (s == 1) {
            res = res * res;
        }
    }
    return res;
}

template <size_t N>
typename Triangle<N>::MatIntNfpNp Triangle<N>::RMatrix(
    const std::size_t s) const {

    matrix::Static<Int, nfp, 1> nodeVec;
    std::size_t last = 0;
    for (std::size_t i = 0; i < N + 1; i++) {
        last += i;
        nodeVec(i, 0) = last;
    }

    matrix::Static<Int, nfp, np> Raux;
    for (std::size_t i = 0; i < nfp; i++) {
        Raux(i, nodeVec(i, 0)) = 1;
    }

    matrix::Dynamic<Int> P = PMatrix(N, s);
    matrix::Static<Int, np, np> Ps;
    Ps = P;
    matrix::Static<Int, nfp, np> res = Raux * Ps * Ps;
    return res;
}

}

