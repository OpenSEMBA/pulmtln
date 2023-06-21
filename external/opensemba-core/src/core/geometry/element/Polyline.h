#pragma once

#include <vector>

#include "Line.h"

namespace semba {
namespace geometry {
namespace element {

class PolylineBase : public virtual LineBase {
public:
    PolylineBase() {};
    virtual ~PolylineBase() {};
};

template<class T>
class Polyline : public virtual Line<T>,
                 public virtual PolylineBase {
public:
    Polyline(const Id id,
             const std::vector<const coordinate::Coordinate<T,3>*>& v,
             const Layer* lay = nullptr,
             const Model* mat = nullptr);
    Polyline(const Polyline<T>& rhs);
    virtual ~Polyline() = default;
    
    virtual std::unique_ptr<Base> clone() const override {
        return std::make_unique<Polyline>(*this);
    }

    inline std::size_t numberOfCoordinates() const { return v_.size(); }

    const coordinate::Coordinate<T,3>* getV    (const std::size_t i) const;
    const coordinate::Coordinate<T,3>* getSideV(const std::size_t f,
                                                const std::size_t i) const;

    const coordinate::Coordinate<T,3>* getVertex    (
            const std::size_t i) const;
    const coordinate::Coordinate<T,3>* getSideVertex(
            const std::size_t f,
            const std::size_t i) const;

    void setV(const std::size_t i, const coordinate::Coordinate<T,3>* coord);

    std::unique_ptr<ElemI> toStructured(const CoordI3Group&,
        const Grid3&,
        const math::Real = Grid3::tolerance) const;
    std::unique_ptr<ElemR> toUnstructured(const CoordR3Group&,
        const Grid3&) const;

private:
    std::vector<const coordinate::Coordinate<T,3>*> v_;
};


template<class T>
Polyline<T>::Polyline(const Id id,
    const std::vector<const coordinate::Coordinate<T, 3>*>& v,
    const Layer* lay,
    const Model* mat)
    : Identifiable<Id>(id),
    Elem(lay, mat) {
    v_ = v;
}

template<class T>
Polyline<T>::Polyline(const Polyline<T>& rhs)
    : Identifiable<Id>(rhs),
    Elem(rhs) {
    v_ = rhs.v_;
}

template<class T>
const coordinate::Coordinate<T, 3>* Polyline<T>::getV(
    const std::size_t i) const {
    return v_[i];
}

template<class T>
const coordinate::Coordinate<T, 3>* Polyline<T>::getSideV(
    const std::size_t f,
    const std::size_t i) const {
    if (f == 0) {
        return v_.front();
    }
    return v_.back();
}

template<class T>
const coordinate::Coordinate<T, 3>* Polyline<T>::getVertex(
    const std::size_t i) const {
    if (i == 0) {
        return v_.front();
    }
    return v_.back();
}

template<class T>
const coordinate::Coordinate<T, 3>* Polyline<T>::getSideVertex(
    const std::size_t f,
    const std::size_t i) const {
    if (f == 0) {
        return v_.front();
    }
    return v_.back();
}

template<class T>
void Polyline<T>::setV(const std::size_t i,
    const coordinate::Coordinate<T, 3>* coord) {

    assert(i < numberOfCoordinates());
    v_[i] = coord;
}

template<class T>
std::unique_ptr<ElemI> Polyline<T>::toStructured(
    const CoordI3Group& cG,
    const Grid3& grid, const math::Real tol) const {
    throw std::logic_error("Polyline::toStructured operation not permitted");
}

template<class T>
std::unique_ptr<ElemR> Polyline<T>::toUnstructured(
    const CoordR3Group& cG,
    const Grid3& grid) const {
    throw std::logic_error("Polyline::toUnstructured operation not permitted");
}

} 

typedef element::PolylineBase         Polylin;
typedef element::Polyline<math::Real> PolylinR;
typedef element::Polyline<math::Int > PolylinI;

} 
} 