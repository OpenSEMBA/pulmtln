#pragma once

#include <array>

#include "Volume.h"

namespace semba {
namespace geometry {
namespace element {

class Hexahedron8Base : public virtual VolumeBase {
public:
    static const std::size_t sizeOfCoordinates = 8;

    virtual ~Hexahedron8Base() = default;

    inline bool isQuadratic() const { return false; }

    inline std::size_t numberOfFaces      () const { return 6; }
    inline std::size_t numberOfVertices   () const { return 8; }
    inline std::size_t numberOfCoordinates() const { return sizeOfCoordinates; }

    inline std::size_t numberOfSideVertices   (const std::size_t f = 0) const {
        return 4;
    }
    inline std::size_t numberOfSideCoordinates(const std::size_t f = 0) const {
        return 4;
    }
};

template<class T>
class Hexahedron8 : public virtual Volume<T>,
              public virtual Hexahedron8Base {
public:
    Hexahedron8(const Id id,
                const std::array<const coordinate::Coordinate<T,3>*, 8> v,
				const Layer* lay = nullptr,
                const Model* mat = nullptr);
    Hexahedron8(const Id id,
        const coordinate::Coordinate<T, 3>* v[8],
        const Layer* lay = nullptr,
        const Model* mat = nullptr);
    Hexahedron8(coordinate::Group<coordinate::Coordinate<T,3> >&,
                const Id id,
                const Box<T,3>& box,
                const Layer* lay = nullptr,
                const Model* mat = nullptr);
    
    Hexahedron8(const Hexahedron8<T>&) = default;
    Hexahedron8(Hexahedron8<T>&&) = default;
    Hexahedron8& operator=(const Hexahedron8<T>&) = default;
    Hexahedron8& operator=(Hexahedron8<T>&&) = default;
    virtual ~Hexahedron8() = default;

    virtual std::unique_ptr<Base> clone() const override {
        return std::make_unique<Hexahedron8>(*this);
    }

    bool isStructured(const Grid3&, const math::Real = Grid3::tolerance) const override;

    bool isRegular() const;
    inline bool isCurvedFace(const std::size_t f) const override{ return false; }

    const coordinate::Coordinate<T,3>* getV(const std::size_t i) const override;
    const coordinate::Coordinate<T,3>* getSideV(
        const std::size_t f, const std::size_t i) const override;

    const coordinate::Coordinate<T,3>* getVertex(const std::size_t i) const override;
    const coordinate::Coordinate<T,3>* getSideVertex(
        const std::size_t f, const std::size_t i) const override;

    std::vector<const coordinate::Coordinate<T,3>*> getVertices() const;
    std::vector<const coordinate::Coordinate<T,3>*> getSideVertices(
            const std::size_t face) const;

    math::Real getAreaOfFace(const std::size_t face) const override;
    math::Real getVolume() const override;

    void setV(const std::size_t i, const coordinate::Coordinate<T,3>*) override;

    std::unique_ptr<ElemI> toStructured(
        const coordinate::Group<CoordI3>&, const Grid3&, const math::Real = Grid3::tolerance) const override;
    std::unique_ptr<ElemR> toUnstructured(
        const coordinate::Group<CoordR3>&, const Grid3&) const override;

private:
    std::array<const coordinate::Coordinate<T,3>*, 8> v_;

    const static math::Real tolerance;
};

template<class T>
const math::Real Hexahedron8<T>::tolerance = 1e-15;

template<class T>
Hexahedron8<T>::Hexahedron8(const Id id,
    const std::array<const coordinate::Coordinate<T, 3>*,8> v,
    const Layer* lay, const Model* mat): Identifiable<Id>(id),
    Elem(lay, mat),
    v_{v}
{
}


template<class T>
Hexahedron8<T>::Hexahedron8(
    const Id id,
    const coordinate::Coordinate<T, 3>* v[8],
    const Layer* lay,
    const Model* mat)
{
    std::array<const coordinate::Coordinate<T, 3>*, 8> vArr;
    std::copy(v, v + 8, vArr.begin());
    *this = Hexahedron8<T>(id, vArr, lay, mat);
}

template<class T>
Hexahedron8<T>::Hexahedron8(
    coordinate::Group<coordinate::Coordinate<T, 3> >& cG,
    const Id id,
    const Box<T, 3>& box,
    const Layer* lay,
    const Model* mat)
    : Identifiable<Id>(id),
    Elem(lay, mat) {

    if (!box.isVolume()) {
        throw geometry::Error::Box::NotVolume();
    }
    std::vector<math::CartesianVector<T, 3> > pos = box.getPos();
    for (std::size_t i = 0; i < numberOfCoordinates(); i++) {
        v_[i] = cG.addPos(pos[i])->get();
    }
}

template<class T>
bool Hexahedron8<T>::isStructured(const Grid3& grid,
    const math::Real tol) const {
    if (!this->vertexInCell(grid, tol)) {
        return false;
    }
    if (!this->getBound().isVolume()) {
        return false;
    }
    if (!this->vertexInBound()) {
        return false;
    }
    return true;
}

template<class T>
bool Hexahedron8<T>::isRegular() const {
    // Checks that all edges are aligned with one of the axis.
    static const math::CartesianVector<T, 3> xAxe(1.0, 0.0, 0.0);
    static const math::CartesianVector<T, 3> yAxe(0.0, 1.0, 0.0);
    static const math::CartesianVector<T, 3> zAxe(1.0, 0.0, 1.0);
    for (std::size_t f = 0; f < numberOfFaces(); f++) {
        math::CartesianVector<T, 3> first, second;
        math::CartesianVector<T, 3> inc;
        for (std::size_t i = 0; i < numberOfSideVertices(); i++) {
            first = *getSideV(f, i);
            second = *getSideV(f, (i + 1) % numberOfSideVertices());
            inc = (second - first).normalize();
            if (!((std::abs(inc.dot(xAxe)) - 1.0) <= tolerance ||
                (std::abs(inc.dot(yAxe)) - 1.0) <= tolerance ||
                (std::abs(inc.dot(zAxe)) - 1.0) <= tolerance)) {

                return false;
            }
        }
    }
    return true;
}

template<class T>
const coordinate::Coordinate<T, 3>* Hexahedron8<T>::getV(
    const std::size_t i) const {
    return v_[i];
}

template<class T>
const coordinate::Coordinate<T, 3>* Hexahedron8<T>::getSideV(
    const std::size_t f, const std::size_t i) const {
    assert(f < numberOfFaces());
    assert(i < numberOfSideCoordinates());
    static const std::array<
        std::array<std::size_t, 4>, 6>
        index = { {{0,1,2,3},
                  {0,3,7,4},
                  {0,4,5,1},
                  {1,5,6,2},
                  {2,6,7,3},
                  {4,7,6,5}} };
    return v_[index[f][i]];
}

template<class T>
const coordinate::Coordinate<T, 3>* Hexahedron8<T>::getVertex(
    const std::size_t i) const {
    return v_[i];
}

template<class T>
const coordinate::Coordinate<T, 3>* Hexahedron8<T>::getSideVertex(
    const std::size_t f,
    const std::size_t i) const {
    return getSideV(f, i);
}

template<class T>
math::Real Hexahedron8<T>::getAreaOfFace(const std::size_t f) const {
    math::CartesianVector<T, 3> v1, v2;
    v1 = getSideV(f, 1)->pos() - getSideV(f, 0)->pos();
    v2 = getSideV(f, 2)->pos() - getSideV(f, 0)->pos();
    return ((math::Real)0.5 * (v1 ^ v2).norm());
}

template<class T>
math::Real Hexahedron8<T>::getVolume() const {
    throw std::logic_error("Hexahedron8::getVolume() not implemented");
}

template<class T>
void Hexahedron8<T>::setV(const std::size_t i,
    const coordinate::Coordinate<T, 3>* coord) {
    assert(i < numberOfCoordinates());
    v_[i] = coord;
}

template<class T>
std::unique_ptr<ElemI> Hexahedron8<T>::toStructured(const CoordI3Group& cG,
    const Grid3& grid,
    const math::Real tol) const {
    return std::make_unique<Hexahedron8<math::Int>>(this->getId(),
        this->vertexToStructured(cG, grid, tol).data(),
        this->getLayer(),
        this->getModel());
}

template<class T>
std::unique_ptr<ElemR> Hexahedron8<T>::toUnstructured(const CoordR3Group& cG,
    const Grid3& grid) const {
    return std::make_unique<Hexahedron8<math::Real>>(this->getId(),
        this->vertexToUnstructured(cG, grid).data(),
        this->getLayer(),
        this->getModel());
}


} 

typedef element::Hexahedron8Base         Hex8;
typedef element::Hexahedron8<math::Real> HexR8;
typedef element::Hexahedron8<math::Int > HexI8;

} 
} 

