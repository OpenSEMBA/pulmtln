#pragma once

#include <iostream>
#include <complex>
#include <stdexcept>
#include <cassert>

#include "core/math/Types.h"
#include "core/math/Constants.h"
#include "Real.h"

namespace semba {
namespace math {

template <class T, std::size_t D>
class CartesianVector {
public:
    // TODO: Remove plain array
    T val[D];
    CartesianVector();
    CartesianVector<T,D>(const T val_);
    CartesianVector<T,D>(T val_[D]);
    CartesianVector<T,D>(const T val_[D]);
    CartesianVector<T,D>(const T, const T, const T);
    CartesianVector<T,D>(const CartesianVector<T,D>&,
                   const CartesianVector<T,D>&);
    template<class U>
    CartesianVector<T,D>(const CartesianVector<U,D>&);
    virtual ~CartesianVector();

    CartesianVector<T,D>& operator= (const T);

    template<class U>
    CartesianVector<T,D>& operator= (const CartesianVector<U,D>&);

    CartesianVector<T,D>& operator+=(const T param);
    CartesianVector<T,D>& operator+=(const CartesianVector<T,D>&);
    CartesianVector<T,D>& operator-=(const T param);
    CartesianVector<T,D>& operator-=(const CartesianVector<T,D>&);
    CartesianVector<T,D>& operator*=(const T param);
    CartesianVector<T,D>& operator*=(const CartesianVector<T,D>&);
    CartesianVector<T,D>& operator/=(const T param);

    CartesianVector<T,D>  operator+(const T param) const;
    CartesianVector<T,D>  operator+(const CartesianVector<T,D>& param) const;
    CartesianVector<T,D>& operator-();
    CartesianVector<T,D>  operator-(const T param) const;
    CartesianVector<T,D>  operator-(const CartesianVector<T,D>& param) const;
    CartesianVector<T,D>  operator*(const T param) const;
    CartesianVector<T,D>  operator*(const CartesianVector<T,D>& param) const;
    CartesianVector<T,D>  operator/(const T param) const;
    CartesianVector<T,D>  operator^(const CartesianVector<T,D>& param) const;

    T dot(const CartesianVector<T,D>& param) const;
    T getMax() const;

    bool operator==(const CartesianVector<T,D>& param) const;
    bool operator!=(const CartesianVector<T,D>& param) const;
    bool isContainedInPlane() const;
    bool isContainedInPlane(
            const Constants::CartesianPlane plane) const;
    bool isCoplanar(const CartesianVector<T,D>& param) const;

    bool                     isInCartesianAxis() const;
    Constants::CartesianAxis getCartesianAxis() const;

    T& operator() (std::size_t pos);
    T  operator() (std::size_t pos) const;

    T& operator[] (std::size_t pos);
    T  operator[] (std::size_t pos) const;

    CartesianVector<T,D>& setAsBinary(std::size_t number);
    CartesianVector<T,D>& setWithMinimalComponents(const CartesianVector<T,D>& rhs);

    Real norm() const;

    CartesianVector<T,D>& abs();
    CartesianVector<T,D>& normalize();
    CartesianVector<T,D>& setPlusInfty();
    CartesianVector<T,D>& setMinusInfty();

    CartesianVector<T,D>& cyclicPermutation(const std::size_t n=1);

    std::string toStr() const;
};

template<std::size_t D>
CartesianVector<Real,D> operator+(const CartesianVector<Int ,D>& lhs,
                            const CartesianVector<Real,D>& rhs);
template<std::size_t D>
CartesianVector<Real,D> operator/(const CartesianVector<Int,D>& lhs,
                            const Real rhs);
template<class T, std::size_t D>
bool operator< (const CartesianVector<T,D>& lhs,
                const CartesianVector<T,D>& rhs);
template<class T, std::size_t D>
bool operator<=(const CartesianVector<T,D>& lhs,
                const CartesianVector<T,D>& rhs);
template<class T, std::size_t D>
bool operator> (const CartesianVector<T,D>& lhs,
                const CartesianVector<T,D>& rhs);
template<class T, std::size_t D>
bool operator>=(const CartesianVector<T,D>& lhs,
                const CartesianVector<T,D>& rhs);

template <class T, std::size_t D>
std::ostream& operator<<(std::ostream& os, const CartesianVector<T,D>& vec) {
    return os << vec.toStr();
}


template <class T, std::size_t D>
CartesianVector<T, D>::CartesianVector() {
    for (std::size_t i = 0; i < D; i++) {
        val[i] = T(0);
    }
}

template <class T, std::size_t D>
CartesianVector<T, D>::CartesianVector(const T val_) {
    for (std::size_t i = 0; i < D; i++) {
        val[i] = val_;
    }
}

template<class T, std::size_t D>
CartesianVector<T, D>::CartesianVector(T val_[D]) {
    for (std::size_t i = 0; i < D; i++) {
        val[i] = val_[i];
    }
}

template<class T, std::size_t D>
CartesianVector<T, D>::CartesianVector(const T val_[D]) {
    for (std::size_t i = 0; i < D; i++) {
        val[i] = val_[i];
    }
}


template <class T, std::size_t D>
CartesianVector<T, D>::CartesianVector(const T x, const T y, const T z) {
    assert(D == 3);
    val[0] = x;
    val[1] = y;
    val[2] = z;
}

template <class T, std::size_t D>
CartesianVector<T, D>::CartesianVector(const CartesianVector<T, D>& begin,
    const CartesianVector<T, D>& end) {
    for (std::size_t i = 0; i < D; i++) {
        val[i] = end.val[i] - begin.val[i];
    }
}

template <class T, std::size_t D> template<class U>
CartesianVector<T, D>::CartesianVector(const CartesianVector<U, D>& param) {
    for (std::size_t i = 0; i < D; i++) {
        val[i] = (T)param.val[i];
    }
}

template<class T, std::size_t D>
CartesianVector<T, D>::~CartesianVector() {

}

template <class T, std::size_t D>
CartesianVector<T, D>& CartesianVector<T, D>::operator=(const T param) {
    for (std::size_t i = 0; i < D; i++) {
        val[i] = param;
    }
    return *this;
}

template <class T, std::size_t D> template<class U>
CartesianVector<T, D>& CartesianVector<T, D>::operator=(const CartesianVector<U, D>& param) {
    for (std::size_t i = 0; i < D; i++) {
        val[i] = (T)param.val[i];
    }
    return *this;
}
//
//template <class T, std::size_t D>
//Cartesian<T,D>& Cartesian<T,D>::operator=(
//        const Cartesian<std::size_t,D>& param) {
//    for (std::size_t i = 0; i < D; i++) {
//        val[i] = (T) param.val[i];
//    }
//    return *this;
//}

template <class T, std::size_t D>
CartesianVector<T, D>& CartesianVector<T, D>::operator+=(const T param) {
    for (std::size_t i = 0; i < D; i++) {
        val[i] += param;
    }
    return *this;
}

template <class T, std::size_t D>
inline CartesianVector<T, D>& CartesianVector<T, D>::operator+=(
    const CartesianVector<T, D>& param) {
    for (std::size_t i = 0; i < D; i++) {
        val[i] += param.val[i];
    }
    return *this;
}

template <class T, std::size_t D>
CartesianVector<T, D>& CartesianVector<T, D>::operator-=(const T param) {
    for (std::size_t i = 0; i < D; i++) {
        val[i] -= param;
    }
    return *this;
}

template <class T, std::size_t D>
CartesianVector<T, D>& CartesianVector<T, D>::operator-=(const CartesianVector<T, D>& param) {
    for (std::size_t i = 0; i < D; i++) {
        val[i] -= param.val[i];
    }
    return *this;
}

template <class T, std::size_t D>
inline CartesianVector<T, D>& CartesianVector<T, D>::operator*=(const T param) {
    for (std::size_t i = 0; i < D; i++) {
        val[i] *= param;
    }
    return *this;
}

template <class T, std::size_t D>
inline CartesianVector<T, D>& CartesianVector<T, D>::operator*=(const CartesianVector<T, D>& param) {
    for (std::size_t i = 0; i < D; i++) {
        val[i] *= param.val[i];
    }
    return *this;
}

template <class T, std::size_t D>
inline CartesianVector<T, D>& CartesianVector<T, D>::operator/=(const T param) {
    for (std::size_t i = 0; i < D; i++) {
        val[i] /= param;
    }
    return *this;
}

template <class T, std::size_t D>
CartesianVector<T, D> CartesianVector<T, D>::operator+(const T param) const {
    CartesianVector<T, D> res;
    for (std::size_t i = 0; i < D; i++) {
        res.val[i] = val[i] + param;
    }
    return res;
}

template <class T, std::size_t D>
CartesianVector<T, D> CartesianVector<T, D>::operator+(const CartesianVector<T, D>& param) const {
    CartesianVector<T, D> res;
    for (std::size_t i = 0; i < D; i++) {
        res.val[i] = val[i] + param.val[i];
    }
    return res;
}

template <class T, std::size_t D>
CartesianVector<T, D>& CartesianVector<T, D>::operator-() {
    for (std::size_t i = 0; i < D; i++) {
        val[i] = -val[i];
    }
    return *this;
}

template <class T, std::size_t D>
CartesianVector<T, D> CartesianVector<T, D>::operator-(const CartesianVector<T, D>& param) const {
    CartesianVector<T, D> res;
    for (std::size_t i = 0; i < D; i++) {
        res.val[i] = val[i] - param.val[i];
    }
    return res;
}

template <class T, std::size_t D>
CartesianVector<T, D> CartesianVector<T, D>::operator-(const T param) const {
    CartesianVector<T, D> res;
    for (std::size_t i = 0; i < D; i++) {
        res.val[i] = val[i] - param;
    }
    return res;
}

template <class T, std::size_t D>
CartesianVector<T, D> operator-(const T& lhs, const CartesianVector<T, D>& rhs) {
    CartesianVector<Real, D> res;
    for (std::size_t i = 0; i < D; i++) {
        res.val[i] = lhs - rhs.val[i];
    }
    return res;
}

template <class T, std::size_t D>
CartesianVector<T, D> CartesianVector<T, D>::operator*(const T param) const {
    CartesianVector<T, D>  res;
    for (std::size_t i = 0; i < D; i++) {
        res.val[i] = val[i] * param;
    }
    return res;
}

template <class T, std::size_t D> inline
CartesianVector<T, D> CartesianVector<T, D>::operator*(const CartesianVector<T, D>& param) const {
    CartesianVector<T, D>  res;
    for (std::size_t i = 0; i < D; i++) {
        res.val[i] = val[i] * param.val[i];
    }
    return res;
}

template <class T, std::size_t D>
CartesianVector<T, D> CartesianVector<T, D>::operator/(const T param) const {
    CartesianVector<T, D>  res;
    for (std::size_t i = 0; i < D; i++) {
        res.val[i] = val[i] / param;
    }
    return res;
}

template <class T, std::size_t D>
CartesianVector<T, D> CartesianVector<T, D>::operator^(const CartesianVector<T, D>& param) const {
    // PURPOSE: Computes vectorial product.
    assert(D == 3);
    CartesianVector<T, D> res;
    res.val[0] = val[1] * param.val[2] - param.val[1] * val[2];
    res.val[1] = val[2] * param.val[0] - param.val[2] * val[0];
    res.val[2] = val[0] * param.val[1] - param.val[0] * val[1];
    return res;
}

template <class T, std::size_t D>
inline T CartesianVector<T, D>::dot(const CartesianVector<T, D>& param) const {
    T res = T();
    for (std::size_t i = 0; i < D; i++) {
        res += val[i] * param.val[i];
    }
    return res;
}

template <class T, std::size_t D>
inline T CartesianVector<T, D>::getMax() const {
    T res = val[0];
    for (std::size_t i = 1; i < D; i++) {
        if (val[i] > res) {
            res = val[i];
        }
    }
    return res;
}

template <class T, std::size_t D>
CartesianVector<T, D>& CartesianVector<T, D>::setAsBinary(std::size_t number) {
    assert(number < pow(2, D));
    for (std::size_t d = 0; d < D; d++) {
        val[D - d - 1] = number % 2;
        number /= 2;
    }
    return *this;
}

template <class T, std::size_t D>
CartesianVector<T, D>& CartesianVector<T, D>::setWithMinimalComponents(
    const CartesianVector<T, D>& rhs) {
    for (std::size_t d = 0; d < D; d++) {
        if (val[d] > rhs.val[d]) {
            val[d] = rhs.val[d];
        }
    }
    return *this;
}

template <class T, std::size_t D>
bool CartesianVector<T, D>::operator==(const CartesianVector<T, D>& param) const {
    return equal((*this - param).norm(), 0.0, (*this + param).norm());
}

template <class T, std::size_t D>
inline bool CartesianVector<T, D>::operator!=(const CartesianVector<T, D>& param) const {
    return !(*this == param);
}

template <class T, std::size_t D>
bool CartesianVector<T, D>::isContainedInPlane() const {
    return (this->isContainedInPlane(Constants::xy) ||
        this->isContainedInPlane(Constants::yz) ||
        this->isContainedInPlane(Constants::zx));
}

template <class T, std::size_t D>
bool CartesianVector<T, D>::isContainedInPlane(
    const Constants::CartesianPlane plane) const {
    assert(D == 3);
    if (std::is_same<T, std::complex<Real>>::value) {
        switch (plane) {
        case Constants::xy:
            if (equal(std::abs(val[2]), 0.0)) {
                return true;
            }
            break;
        case Constants::yz:
            if (equal(std::abs(val[0]), 0.0)) {
                return true;
            }
            break;
        case Constants::zx:
            if (equal(std::abs(val[1]), 0.0)) {
                return true;
            }
            break;
        }
        return false;
    }
    else {
        switch (plane) {
        case Constants::xy:
            if (equal(std::abs(val[2]), 0.0)) {
                return true;
            }
            break;
        case Constants::yz:
            if (equal(std::abs(val[0]), 0.0)) {
                return true;
            }
            break;
        case Constants::zx:
            if (equal(std::abs(val[1]), 0.0)) {
                return true;
            }
            break;
        }
        return false;
    }
}

template<class T, std::size_t D>
inline bool CartesianVector<T, D>::isCoplanar(const CartesianVector<T, D>& param) const {
    return (*this - param).isContainedInPlane();
}

template <class T, std::size_t D>
inline T& CartesianVector<T, D>::operator() (std::size_t pos) {
    assert(pos >= 0 && pos < D);
    return val[pos];
}

template <class T, std::size_t D>
inline T CartesianVector<T, D>::operator() (std::size_t pos) const {
    assert(pos >= 0 && pos < D);
    return val[pos];
}

template <class T, std::size_t D>
inline T& CartesianVector<T, D>::operator[] (std::size_t pos) {
    assert(pos >= 0 && pos < D);
    return val[pos];
}

template <class T, std::size_t D>
inline T CartesianVector<T, D>::operator[] (std::size_t pos) const {
    assert(pos >= 0 && pos < D);
    return val[pos];
}

template <class T, std::size_t D>
inline Real CartesianVector<T, D>::norm() const {
    Real sum = 0;
    for (std::size_t i = 0; i < D; i++) {
        sum += (Real)std::abs(val[i]) * std::abs(val[i]);
    }
    return sqrt(sum);
}

template <class T, std::size_t D>
inline bool CartesianVector<T, D>::isInCartesianAxis() const {
    try {
        getCartesianAxis();
    }
    catch (const std::logic_error&) {
        return false;
    }
    return true;
}

template <class T, std::size_t D>
inline Constants::CartesianAxis CartesianVector<T, D>::getCartesianAxis() const {
    assert(D <= 3);
    for (size_t dir = 0; dir < D; dir++) {
        CartesianVector<T, D> axis;
        axis(dir) = (T)1.0;
        Real dotProduct = std::abs(this->dot(axis));
        Real norm = this->norm();
        if (equal(dotProduct, norm, 1e-2)) {
            return Constants::CartesianAxis(dir);
        }
    }
    throw std::logic_error("Vector is not in Cartesian Axis");
}


template <class T, std::size_t D>
CartesianVector<T, D>& CartesianVector<T, D>::cyclicPermutation(const std::size_t n) {
    CartesianVector<T, D> valAux(0.0);
    for (std::size_t i = 0; i < D; i++) {
        valAux.val[(i + n) % D] = val[i];
    }
    *this = valAux;

    return *this;
}

template <class T, std::size_t D> inline
CartesianVector<T, D>& CartesianVector<T, D>::abs() {
    for (std::size_t i = 0; i < D; i++) {
        val[i] = (T)std::fabs(val[i]);
    }
    return *this;
}

template <class T, std::size_t D>
CartesianVector<T, D>& CartesianVector<T, D>::normalize() {
    Real nor = norm();
    for (std::size_t i = 0; i < D; i++) {
        val[i] /= (T)nor;
    }
    return *this;
}

template <class T, std::size_t D> inline
CartesianVector<T, D>& CartesianVector<T, D>::setPlusInfty() {
    for (std::size_t i = 0; i < D; i++) {
        if (std::numeric_limits<T>::has_infinity) {
            val[i] = std::numeric_limits<T>::infinity();
        }
        else {
            val[i] = std::numeric_limits<T>::max();
        }
    }
    return *this;
}

template <class T, std::size_t D> inline
CartesianVector<T, D>& CartesianVector<T, D>::setMinusInfty() {
    for (std::size_t i = 0; i < D; i++) {
        if (std::numeric_limits<T>::has_infinity) {
            val[i] = -std::numeric_limits<T>::infinity();
        }
        else {
            val[i] = std::numeric_limits<T>::min();
        }
    }
    return *this;
}

template <class T, std::size_t D>
std::string CartesianVector<T, D>::toStr() const {
    std::stringstream ss;
    ss << "(";
    for (std::size_t i = 0; i < D; i++) {
        ss << val[i];
        if (i < D - 1) {
            ss << " , ";
        }
    }
    ss << ")";
    return ss.str();
}

template<std::size_t D>
CartesianVector<Real, D> operator+(const CartesianVector<Int, D>& lhs,
    const CartesianVector<Real, D>& rhs) {
    CartesianVector<Real, D> res;
    for (std::size_t i = 0; i < D; i++) {
        res(i) = lhs(i) + rhs(i);
    }
    return res;
}

template<std::size_t D>
CartesianVector<Real, D> operator/(const CartesianVector<Int, D>& lhs,
    const Real rhs) {
    CartesianVector<Real, D> res;
    for (std::size_t i = 0; i < D; i++) {
        res(i) = (Real)lhs(i) / rhs;
    }
    return  res;
}

template<class T, std::size_t D>
bool operator< (const CartesianVector<T, D>& lhs,
    const CartesianVector<T, D>& rhs) {
    for (std::size_t i = 0; i < D; i++) {
        if (lower(lhs(i), rhs(i))) {
            return true;
        }
        if (greater(lhs(i), rhs(i))) {
            return false;
        }
    }
    return false;
}

template<class T, std::size_t D>
bool operator<=(const CartesianVector<T, D>& lhs,
    const CartesianVector<T, D>& rhs) {
    return !(rhs < lhs);
}

template<class T, std::size_t D>
bool operator> (const CartesianVector<T, D>& lhs,
    const CartesianVector<T, D>& rhs) {
    return rhs < lhs;
}

template<class T, std::size_t D>
bool operator>=(const CartesianVector<T, D>& lhs,
    const CartesianVector<T, D>& rhs) {
    return !(lhs < rhs);
}


template<std::size_t D>
CartesianVector<Real, D> round(const CartesianVector<Real, D>& vec) {
    CartesianVector<Real, D> res;
    for (std::size_t i = 0; i < D; i++) {
        res(i) = round(vec(i));
    }
    return  res;
}

template<std::size_t D>
CartesianVector<Real,D> round(const CartesianVector<Real,D>& vec);

typedef CartesianVector<Real,2> CVecR2;
typedef CartesianVector<Real,3> CVecR3;
typedef CartesianVector<Real,4> CVecR4;
typedef CartesianVector<Int ,2> CVecI2;
typedef CartesianVector<Int ,3> CVecI3;

typedef CartesianVector<std::complex<Real>,3> CVecC3;

} 
} 

