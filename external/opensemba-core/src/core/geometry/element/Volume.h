#pragma once

#include <exception>

#include "core/math/simplex/Simplex.h"
#include "Surface.h"
#include "Element.h"

namespace semba {
namespace geometry {
namespace element {

class VolumeBase : public virtual Base {
public:
    VolumeBase() = default;
    virtual ~VolumeBase() = default;
};

template<class T>
class Volume : public virtual Element<T>,
               public virtual VolumeBase {
public:
    Volume() = default;
    virtual ~Volume() = default;

    bool isLocalFace(const std::size_t f,
                     const Surface<T>& surf) const;
    virtual bool isCurvedFace(const std::size_t face) const = 0;
    virtual bool isFaceContainedInPlane(
            const std::size_t face,
            const math::Constants::CartesianPlane plane) const;

    math::CartesianVector<T,3> getSideNormal(const std::size_t f) const;

    virtual math::Real getAreaOfFace(const std::size_t face) const = 0;
    virtual math::Real getVolume() const = 0;
    Box<T,3> getBoundOfFace(const std::size_t face) const;
    std::size_t getFaceNumber(const Surface<T>*) const;
};

template<class T>
bool Volume<T>::isLocalFace(const std::size_t f,
    const Surface<T>& surf) const {
    return getSideNormal(f) == surf.getNormal();
}

template<class T>
bool Volume<T>::isFaceContainedInPlane(
    const std::size_t face,
    const math::Constants::CartesianPlane plane) const {
    Box<T, 3> box = getBoundOfFace(face);
    math::CartesianVector<T, 3> vec = box.getMax() - box.getMin();
    return vec.isContainedInPlane(plane);
}

template<class T>
math::CartesianVector<T, 3> Volume<T>::getSideNormal(
    const std::size_t f) const {
    math::CartesianVector<T, 3> vec1, vec2, res;
    vec1 = *this->getSideVertex(f, 1) - *this->getSideVertex(f, 0);
    vec2 = *this->getSideVertex(f, 2) - *this->getSideVertex(f, 0);
    res = (vec1 ^ vec2).normalize();
    return res;
}

template<class T>
Box<T, 3> Volume<T>::getBoundOfFace(const std::size_t face) const {
    Box<T, 3> res;
    for (std::size_t i = 0; i < this->numberOfSideCoordinates(); i++) {
        res << this->getSideV(face, i)->pos();
    }
    return res;
}

template<class T>
std::size_t Volume<T>::getFaceNumber(const Surface<T>* surf) const {
    // Checks each face. Order is not important.
    for (std::size_t f = 0; f < this->numberOfFaces(); f++) {
        std::size_t vPresent = 0;
        for (std::size_t i = 0; i < surf->numberOfVertices(); i++) {
            for (std::size_t j = 0; j < surf->numberOfVertices(); j++) {
                if (surf->getVertex(j) == this->getSideVertex(f, i)) {
                    vPresent++;
                }
            }
            if (vPresent == surf->numberOfVertices()) {
                return f;
            }
        }
    }
    throw std::logic_error("Surface not found");;
}

} 

typedef element::VolumeBase         Vol;
typedef element::Volume<math::Real> VolR;
typedef element::Volume<math::Int > VolI;

} 
} 


