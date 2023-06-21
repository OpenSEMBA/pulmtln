#pragma once

#include <algorithm>
#include <exception>
#include <utility>

#include "core/math/CartesianVector.h"
#include "core/geometry/Grid.h"

namespace semba {
namespace geometry {

template<std::size_t D> class Grid;

template <class T, std::size_t D>
class Box {
    typedef math::CartesianVector<T,D> CVecTD;
public:
    Box();
    Box(const std::pair<CVecTD,CVecTD>& boundsMinMax);
    Box(const CVecTD& min, const CVecTD& max);

	template <typename T2>
	Box<T,D>(const Box<T2, D>& rhs) {
		for (size_t i = 0; i < D; ++i) {
			min_(i) = (T)rhs.getMin()(i);
			max_(i) = (T)rhs.getMax()(i);
		}
	}

    virtual ~Box();

    Box<T,D>& operator= (const Box<T,D>& lBoxMin);
	
	template <typename T2>
	Box<T,D>& operator= (const Box<T2, D>& rhs) {
		for (size_t i = 0; i < D; ++i) {
			min_(i) = (T) rhs.getMin()(i);
			max_(i) = (T) rhs.getMax()(i);
		}
		return *this;
	}

    bool operator>(const Box<T,D>& lBoxMin) const;
    bool operator<(const Box<T,D>& lBoxMax) const;
    bool operator>=(const Box<T,D>& lBoxMin) const;
    bool operator<=(const Box<T,D>& lBoxMax) const;
    bool operator==(const Box<T,D>& lBoxMin) const;

    void operator+=(const Box<T,D>& lBoxSource);
    void operator+(const Box<T,D>& lBoxSource);

    void operator<<(const CVecTD& p);
    void operator<<(const Box<T,D>& p);

    bool isPoint  () const { return (numberOfDifferentCoords() == 0); }
    bool isLine   () const { return (numberOfDifferentCoords() == 1); }
    bool isSurface() const { return (numberOfDifferentCoords() == 2); }
    bool isVolume () const { return (numberOfDifferentCoords() == 3); }

    math::Constants::CartesianAxis getDirection() const;
    math::Constants::CartesianAxis getNormal   () const;

    Box<T,D> intersect(const Box<T,D>& rhs) const;
    bool isIntersected (const Box<T,D>& lBox) const;
    bool isInnerPoint(const CVecTD& point) const;

    inline CVecTD getMin() const;
    inline CVecTD getMax() const;
    inline CVecTD getLength() const;

    std::vector<CVecTD> getPos() const;
    std::vector<CVecTD> getPosOfBound(math::Constants::CartesianAxis  d,
                                      math::Constants::CartesianBound p) const;
    Box<T,D> getBoundAsBox(math::Constants::CartesianAxis  d,
                           math::Constants::CartesianBound p) const;
    CVecTD getBound(math::Constants::CartesianBound p) const;
    std::vector<Box<T,D>> chop(const CVecTD step = CVecTD(1,1,1)) const;
    std::vector<Box<T,D>> chop(const Grid<D>& grid) const;

    void set(const std::pair<CVecTD,CVecTD>& boundsMinMax);
    Box<T,D>& setInfinity();
    void scale(const math::Real factor);

private:
    CVecTD min_, max_;

    void setDefaultValues() {
        min_.setPlusInfty();
        max_.setMinusInfty();
    };

    std::size_t numberOfDifferentCoords() const;
};

template <class T, std::size_t D>
std::ostream& operator<<(std::ostream& os, const Box<T,D>& rhs) {
   return os << "Min: " << rhs.getMin() << ", Max: " << rhs.getMax();
}


template<class T, std::size_t D>
Box<T, D>::Box() {
    setDefaultValues();
}

template<class T, std::size_t D>
Box<T, D>::Box(const std::pair<CVecTD, CVecTD>& bounds) {
    set(bounds);
}

template<class T, std::size_t D>
Box<T, D>::Box(const CVecTD& minB, const CVecTD& maxB) {
    set(std::make_pair(minB, maxB));
}

template<class T, std::size_t D>
Box<T, D>::~Box() {
}

template<class T, std::size_t D>
Box<T, D>& Box<T, D>::operator=(const Box<T, D>& rhs) {
    if (&rhs == this) {
        return *this;
    }
    min_ = rhs.min_;
    max_ = rhs.max_;
    return *this;
}

template<class T, std::size_t D>
bool Box<T, D>::operator>(const Box<T, D>& rhs) const {
    for (std::size_t i = 0; i < D; i++) {
        if ((max_(i) <= rhs.max_(i)) || (min_(i) >= rhs.min_(i))) {
            return false;
        }
    }
    return true;
}

template<class T, std::size_t D>
bool Box<T, D>::operator<(const Box<T, D>& lBoxMax) const {
    return lBoxMax > *this;
}

template<class T, std::size_t D>
bool
Box<T, D>::operator <=(const Box<T, D>& rhs) const {
    for (std::size_t i = 0; i < D; i++) {
        if ((max_(i) > rhs.max_(i)) || (min_(i) < rhs.min_(i))) {
            return false;
        }
    }
    return true;
}

template<class T, std::size_t D>
bool Box<T, D>::operator==(const Box<T, D>& rhs) const {
    if (max_ != rhs.max_) { return false; }
    if (min_ != rhs.min_) { return false; }
    return true;
}

template<class T, std::size_t D>
void Box<T, D>::operator+=(const Box<T, D>& rhs) {
    for (std::size_t i = 0; i < D; i++) {
        if (min_(i) > rhs.min_(i)) {
            min_(i) = rhs.min_(i);
        }
        if (max_(i) < rhs.max_(i)) {
            max_(i) = rhs.max_(i);
        }
    }
}

template<class T, std::size_t D>
void Box<T, D>::operator << (const CVecTD& p) {
    for (std::size_t i = 0; i < D; i++) {
        if (min_(i) > p(i)) {
            min_(i) = p(i);
        }
        if (max_(i) < p(i)) {
            max_(i) = p(i);
        }
    }
}

template<class T, std::size_t D>
void Box<T, D>::operator<<(const Box<T, D>& p) {
    *this << p.min_;
    *this << p.max_;
}

template<class T, std::size_t D>
math::Constants::CartesianAxis Box<T, D>::getDirection() const {
    if (!isLine()) {
        throw std::logic_error("Box object is not a line");
    }
    math::Constants::CartesianAxis res = math::Constants::x;
    for (std::size_t d = 0; d < D; d++) {
        if (math::notEqual(max_(d), min_(d))) {
            res = math::Constants::CartesianAxis(d);
            break;
        }
    }
    return res;
}

template<class T, std::size_t D>
math::Constants::CartesianAxis Box<T, D>::getNormal() const {
    if (!isSurface()) {
        throw std::logic_error("Box object is not a surface");
    }
    assert(D == 3);
    math::Constants::CartesianAxis res = math::Constants::x;
    for (std::size_t d = 0; d < D; d++) {
        if (math::equal(max_(d), min_(d))) {
            res = math::Constants::CartesianAxis(d);
            break;
        }
    }
    return res;
}

template<class T, std::size_t D>
inline Box<T, D> Box<T, D>::intersect(const Box<T, D>& rhs) const {
    Box<T, D> res;
    for (std::size_t d = 0; d < D; d++) {
        if (min_(d) > rhs.min_(d)) {
            res.min_(d) = min_(d);
        }
        else {
            res.min_(d) = rhs.min_(d);
        }
        if (max_(d) < rhs.max_(d)) {
            res.max_(d) = max_(d);
        }
        else {
            res.max_(d) = rhs.max_(d);
        }
    }
    return res;
}

template<class T, std::size_t D>
bool Box<T, D>::isIntersected(const Box<T, D>& rhs) const {
    for (std::size_t i = 0; i < D; i++) {
        if (max_(i) < rhs.min_(i)) {
            return false;
        }
        if (rhs.max_(i) < min_(i)) {
            return false;
        }
    }
    return true;
}

template<class T, std::size_t D>
bool Box<T, D>::isInnerPoint(const CVecTD& point) const {
    bool isInner = true;
    for (std::size_t dir = 0; dir < D; dir++) {
        isInner &= (point(dir) <= getMax()(dir));
        isInner &= (point(dir) >= getMin()(dir));
    }
    return isInner;
}

template<class T, std::size_t D>
inline math::CartesianVector<T, D> Box<T, D>::getMin() const {
    return min_;
}

template<class T, std::size_t D>
inline math::CartesianVector<T, D> Box<T, D>::getMax() const {
    return max_;
}

template<class T, std::size_t D>
inline math::CartesianVector<T, D> Box<T, D>::getLength() const {
    return (max_ - min_);
}

template<class T, std::size_t D>
inline std::vector<math::CartesianVector<T, D> > Box<T, D>::getPos() const {
    assert(D == 3);
    std::vector<CVecTD> res;
    // TODO: Generalize this...
    if (isPoint()) {
        res.push_back(min_);
        return res;
    }
    else if (isLine()) {
        res.resize(2);
        res[0] = min_;
        res[1] = min_;
        for (std::size_t d = 0; d < D; d++) {
            if (math::notEqual(max_(d), min_(d))) {
                res[0](d) = min_(d);
                res[1](d) = max_(d);
            }
        }
    }
    else if (isSurface()) {
        res.resize(4);
        switch (getNormal()) {
        case math::Constants::x:
            res[0] = CVecTD(min_(math::Constants::x),
                min_(math::Constants::y),
                min_(math::Constants::z));
            res[1] = CVecTD(min_(math::Constants::x),
                max_(math::Constants::y),
                min_(math::Constants::z));
            res[2] = CVecTD(min_(math::Constants::x),
                max_(math::Constants::y),
                max_(math::Constants::z));
            res[3] = CVecTD(min_(math::Constants::x),
                min_(math::Constants::y),
                max_(math::Constants::z));
            break;
        case math::Constants::y:
            res[0] = CVecTD(min_(math::Constants::x),
                min_(math::Constants::y),
                min_(math::Constants::z));
            res[1] = CVecTD(min_(math::Constants::x),
                min_(math::Constants::y),
                max_(math::Constants::z));
            res[2] = CVecTD(max_(math::Constants::x),
                min_(math::Constants::y),
                max_(math::Constants::z));
            res[3] = CVecTD(max_(math::Constants::x),
                min_(math::Constants::y),
                min_(math::Constants::z));
            break;
        case math::Constants::z:
            res[0] = CVecTD(min_(math::Constants::x),
                min_(math::Constants::y),
                min_(math::Constants::z));
            res[1] = CVecTD(max_(math::Constants::x),
                min_(math::Constants::y),
                min_(math::Constants::z));
            res[2] = CVecTD(max_(math::Constants::x),
                max_(math::Constants::y),
                min_(math::Constants::z));
            res[3] = CVecTD(min_(math::Constants::x),
                max_(math::Constants::y),
                min_(math::Constants::z));
            break;
        }
    }
    else if (isVolume()) {
        res.resize(8);
        res[0] = CVecTD(min_(math::Constants::x),
            min_(math::Constants::y),
            min_(math::Constants::z));
        res[1] = CVecTD(max_(math::Constants::x),
            min_(math::Constants::y),
            min_(math::Constants::z));
        res[2] = CVecTD(max_(math::Constants::x),
            max_(math::Constants::y),
            min_(math::Constants::z));
        res[3] = CVecTD(min_(math::Constants::x),
            max_(math::Constants::y),
            min_(math::Constants::z));
        res[4] = CVecTD(min_(math::Constants::x),
            min_(math::Constants::y),
            max_(math::Constants::z));
        res[5] = CVecTD(max_(math::Constants::x),
            min_(math::Constants::y),
            max_(math::Constants::z));
        res[6] = CVecTD(max_(math::Constants::x),
            max_(math::Constants::y),
            max_(math::Constants::z));
        res[7] = CVecTD(min_(math::Constants::x),
            max_(math::Constants::y),
            max_(math::Constants::z));
    }
    return res;
}

template<class T, std::size_t D>
inline std::vector<math::CartesianVector<T, D> > Box<T, D>::getPosOfBound(
    math::Constants::CartesianAxis d,
    math::Constants::CartesianBound p) const {

    assert(D == 3);
    std::size_t rX = d;
    std::size_t rY = (d + 1) % D;
    std::size_t rZ = (d + 2) % D;
    std::vector<CVecTD> res;
    // TODO: Generalize this...
    {
        CVecTD aux;
        aux(rX) = getBound(p)(rX);
        aux(rY) = getBound(math::Constants::L)(rY);
        aux(rZ) = getBound(math::Constants::L)(rZ);
        res.push_back(aux);
    }
    {
        CVecTD aux;
        aux(rX) = getBound(p)(rX);
        aux(rY) = getBound(math::Constants::L)(rY);
        aux(rZ) = getBound(math::Constants::U)(rZ);
        res.push_back(aux);
    }
    {
        CVecTD aux;
        aux(rX) = getBound(p)(rX);
        aux(rY) = getBound(math::Constants::U)(rY);
        aux(rZ) = getBound(math::Constants::U)(rZ);
        res.push_back(aux);
    }
    {
        CVecTD aux;
        aux(rX) = getBound(p)(rX);
        aux(rY) = getBound(math::Constants::U)(rY);
        aux(rZ) = getBound(math::Constants::L)(rZ);
        res.push_back(aux);
    }
    return res;
}

template<class T, std::size_t D>
Box<T, D> Box<T, D>::getBoundAsBox(math::Constants::CartesianAxis d,
    math::Constants::CartesianBound p) const {
    std::vector<math::CartesianVector<T, D>> pos = getPosOfBound(d, p);
    assert(pos.size() == 4);
    Box<T, D> res(pos[0], pos[1]);
    for (std::size_t i = 2; i < pos.size(); i++) {
        res << pos[i];
    }
    assert(res.isSurface());
    return res;
}

template<class T, std::size_t D>
inline math::CartesianVector<T, D> Box<T, D>::getBound(
    math::Constants::CartesianBound p) const {
    if (p == math::Constants::L) {
        return getMin();
    }
    else {
        return getMax();
    }
}

template<class T, std::size_t D>
std::vector<Box<T, D>> Box<T, D>::chop(const CVecTD origStep) const {
    static_assert(D == 3, "Chop can't be used for Boxes with D != 3");
    CVecTD length = getLength();
    CVecTD step = origStep;
    for (std::size_t d = 0; d < D; d++) {
        if (length(d) < origStep(d)) {
            step(d) = length(d);
        }
        if (origStep(d) <= (T)0) {
            step(d) = length(d);
        }
    }
    math::CartesianVector<math::Real, 3> minR, maxR, stepR;
    for (std::size_t d = 0; d < D; d++) {
        stepR(d) = step(d);
        if (stepR(d) == 0.0) {
            minR(d) = (T)min_(d);
            maxR(d) = (T)max_(d);
        }
        else {
            minR(d) = (T)min_(d);
            minR(d) = floor(minR(d) / stepR(d)) * stepR(d);
            maxR(d) = (T)max_(d);
            maxR(d) = floor(maxR(d) / stepR(d)) * stepR(d);
        }
    }
    Grid<D> grid(Box<math::Real, 3>(minR, maxR), stepR);
    return chop(grid);
}

template<class T, std::size_t D>
std::vector<Box<T, D>> Box<T, D>::chop(const Grid<D>& grid) const {
    static_assert(D == 3, "Chop can't be used for Boxes with D != 3");
    // Prepares subgrid with the size of the box preserving grid positions.
    std::vector<math::Real> pos[D];
    for (std::size_t d = 0; d < D; d++) {
        pos[d] = grid.getPosInRange(d, min_(d), max_(d));
        if (min_(d) != pos[d].front()) {
            std::vector<math::Real> aux(1, min_(d));
            aux.insert(aux.end(), pos[d].begin(), pos[d].end());
            pos[d] = aux;
        }
        if (max_(d) != pos[d].back()) {
            pos[d].push_back(max_(d));
        }
    }
    Grid<D> subGrid;
    subGrid.setPos(pos);
    //
    math::CartesianVector<math::Int, D> numBoxes = subGrid.getNumCells();
    std::vector<Box<T, D>> res;
    res.reserve(numBoxes(math::Constants::x) *
        numBoxes(math::Constants::y) *
        numBoxes(math::Constants::z));
    for (math::Int i = 0; i < numBoxes(math::Constants::x); i++) {
        for (math::Int j = 0; j < numBoxes(math::Constants::y); j++) {
            for (math::Int k = 0; k < numBoxes(math::Constants::z); k++) {
                CVecTD minP = subGrid.getPos(CVecTD(i, j, k));
                CVecTD maxP = subGrid.getPos(CVecTD(i + 1, j + 1, k + 1));
                CVecTD minT, maxT;
                for (std::size_t d = 0; d < D; d++) {
                    if (minP(d) >= min_(d)) {
                        minT(d) = (T)minP(d);
                    }
                    else {
                        minT(d) = min_(d);
                    }
                    if (maxP(d) <= max_(d)) {
                        maxT(d) = (T)maxP(d);
                    }
                    else {
                        maxT(d) = max_(d);
                    }
                }
                res.push_back(Box<T, D>(minT, maxT));
            }
        }
    }
    return res;
}

template<class T, std::size_t D>
void Box<T, D>::set(const std::pair<CVecTD, CVecTD>& minMax) {
    for (std::size_t d = 0; d < D; d++) {
        min_(d) = std::min(minMax.first(d), minMax.second(d));
        max_(d) = std::max(minMax.first(d), minMax.second(d));
    }
}

template<class T, std::size_t D>
inline Box<T, D>& Box<T, D>::setInfinity() {
    for (std::size_t j = 0; j < D; j++) {
        min_(j) = -std::numeric_limits<T>::infinity();
        max_(j) = std::numeric_limits<T>::infinity();
    }
    return *this;
}

template<class T, std::size_t D>
void Box<T, D>::scale(const math::Real factor) {
    if (std::is_floating_point<T>::value) {
        min_ *= factor;
        max_ *= factor;
    }
}

template<class T, std::size_t D>
std::size_t Box<T, D>::numberOfDifferentCoords() const {
    std::size_t res = 0;
    for (std::size_t d = 0; d < D; d++) {
        if (math::notEqual(max_(d), min_(d))) {
            res++;
        }
    }
    return res;
}

typedef Box<math::Real,3> BoxR3;
typedef Box<math::Int ,2> BoxI2;
typedef Box<math::Int, 3> BoxI3;

} 
} 


