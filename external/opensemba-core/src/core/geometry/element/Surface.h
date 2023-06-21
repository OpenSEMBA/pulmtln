

#pragma once

#include "Element.h"

namespace semba {
namespace geometry {
namespace element {

class SurfaceBase : public virtual Base {
public:
    SurfaceBase() = default;
    virtual ~SurfaceBase() = default;

    virtual math::CVecR3 getNormal() const = 0;
};

template<class T>
class Surface : public virtual Element<T>,
                public virtual SurfaceBase {
public:
    Surface() = default;
    virtual ~Surface() = default;

    bool isRectangular() const;
    bool isContainedInPlane() const;
    bool isContainedInPlane(const math::Constants::CartesianPlane plane) const;

    virtual math::CVecR3 getNormal() const override; 
};

template<class T>
bool Surface<T>::isRectangular() const {
    if (this->numberOfCoordinates() != 4 || this->numberOfFaces() != 4) {
        return false;
    }
    for (std::size_t f = 0; f < 4; f++) {
        math::CartesianVector<T, 3> p0 = this->getSideVertex(f % 4, 0)->pos();
        math::CartesianVector<T, 3> p1 = this->getSideVertex(f % 4, 1)->pos();
        math::CartesianVector<T, 3> p2 = this->getSideVertex((f + 1) % 4, 1)->pos();
        math::Real sProd = (math::Real)(p2 - p1).dot(p1 - p0);
        if (math::greater(sProd, 0.0, 1.0)) {
            return false;
        }
    }
    return true;
}

template<class T>
bool Surface<T>::isContainedInPlane() const {
    return (isContainedInPlane(math::Constants::xy) ||
        isContainedInPlane(math::Constants::yz) ||
        isContainedInPlane(math::Constants::zx));
}

template<class T>
bool Surface<T>::isContainedInPlane(const math::Constants::CartesianPlane plane) const {
    for (std::size_t i = 1; i < this->numberOfCoordinates(); i++) {
        if (!(*this->getV(i) - *this->getV(0)).isContainedInPlane(plane)) {
            return false;
        }
    }
    return true;
}

template<class T>
math::CVecR3 Surface<T>::getNormal() const
{
    math::CVecR3 v0 = this->getVertex(1)->pos() - this->getVertex(0)->pos();
    math::CVecR3 v1 = this->getVertex(2)->pos() - this->getVertex(0)->pos();
    return (v0 ^ v1).normalize();
}

} 

typedef element::SurfaceBase         Surf;
typedef element::Surface<math::Real> SurfR;
typedef element::Surface<math::Int > SurfI;

} 
} 
