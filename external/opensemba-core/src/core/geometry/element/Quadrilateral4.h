#pragma once

#include "Quadrilateral.h"

namespace semba {
namespace geometry {
namespace element {

class Quadrilateral4Base : public virtual SurfaceBase {
public:
    static const std::size_t sizeOfCoordinates = 4;

    Quadrilateral4Base() {}
    virtual ~Quadrilateral4Base() {}

    std::size_t numberOfCoordinates() const { return sizeOfCoordinates; }

    std::size_t numberOfSideCoordinates(const std::size_t f = 0) const { 
        return 2; 
    }
};

template<class T>
class Quadrilateral4: public virtual Quadrilateral<T>,
                      public virtual Quadrilateral4Base {
public:
    Quadrilateral4(
        const Id id,
        const coordinate::Coordinate<T,3>* coords[4],
        const Layer* lay = nullptr,
        const Model* mat = nullptr);
	Quadrilateral4(const Id id,
		std::array<const coordinate::Coordinate<T, 3>*, 4> v,
		const Layer* lay = nullptr,
		const Model* mat = nullptr);
    Quadrilateral4(const Quadrilateral4<T>& rhs) = default;
    virtual ~Quadrilateral4() = default;

    virtual std::unique_ptr<Base> clone() const override {
        return std::make_unique<Quadrilateral4<T>>(*this);
    }

    bool isStructured(const Grid3&, const math::Real = Grid3::tolerance) const override; 

    const coordinate::Coordinate<T,3>* getV(const std::size_t i) const override;
    const coordinate::Coordinate<T,3>* getSideV(
        const std::size_t f, const std::size_t i) const override;

    const coordinate::Coordinate<T,3>* getVertex(const std::size_t i) const override;
    const coordinate::Coordinate<T,3>* getSideVertex(
        const std::size_t f, const std::size_t i) const override;

    void setV(const std::size_t i, const coordinate::Coordinate<T,3>*) override;

    std::unique_ptr<ElemI> toStructured(
        const CoordI3Group&, const Grid3&, const math::Real = Grid3::tolerance) const override;
    std::unique_ptr<ElemR> toUnstructured(
        const CoordR3Group&, const Grid3&) const override;

private:
    std::array<const coordinate::Coordinate<T,3>*, 4> v_;
};

template<class T>
Quadrilateral4<T>::Quadrilateral4(const Id id,
    const coordinate::Coordinate<T, 3>* coords[4],
    const Layer* lay,
    const Model* mat): 
    Identifiable<Id>(id),
    Elem(lay, mat) 
{
    for (std::size_t i = 0; i < numberOfCoordinates(); i++) {
        v_[i] = coords[i];
    }
}

template<class T>
Quadrilateral4<T>::Quadrilateral4(const Id id,
    std::array<const coordinate::Coordinate<T, 3>*, 4> v,
    const Layer* lay,
    const Model* mat): 
    Identifiable<Id>(id),
    Elem(lay, mat),
    v_{v}
{}

template<class T>
bool Quadrilateral4<T>::isStructured(const Grid3& grid,
    const math::Real tol) const {
    if (!this->vertexInCell(grid, tol)) {
        return false;
    }
    if (!this->getBound().isSurface()) {
        return false;
    }
    if (!this->vertexInBound()) {
        return false;
    }
    return true;
}

template<class T>
const coordinate::Coordinate<T, 3>* Quadrilateral4<T>::getV(
    const std::size_t i) const {
    return v_[i];
}

template<class T>
const coordinate::Coordinate<T, 3>* Quadrilateral4<T>::getVertex(
    const std::size_t i) const {
    return v_[i];
}

template<class T>
const coordinate::Coordinate<T, 3>* Quadrilateral4<T>::getSideV(
    const std::size_t f,
    const std::size_t i) const {
    assert(f < this->numberOfFaces());
    assert(i < numberOfSideCoordinates());
    return v_[(f + i) % 4];
}

template<class T>
const coordinate::Coordinate<T, 3>* Quadrilateral4<T>::getSideVertex(
    const std::size_t f,
    const std::size_t i) const {
    assert(f < this->numberOfFaces());
    assert(i < this->numberOfSideVertices());
    return v_[(f + i) % 4];
}

template<class T>
void Quadrilateral4<T>::setV(const std::size_t i,
    const coordinate::Coordinate<T, 3>* coord) {
    v_[i] = coord;
}

template<class T>
std::unique_ptr<ElemI> Quadrilateral4<T>::toStructured(
    const CoordI3Group& cG,
    const Grid3& grid, const math::Real tol) const {
    return std::make_unique<Quadrilateral4<math::Int>>(this->getId(),
        this->vertexToStructured(cG, grid, tol).data(),
        this->getLayer(),
        this->getModel());
}

template<class T>
std::unique_ptr<ElemR> Quadrilateral4<T>::toUnstructured(
    const CoordR3Group& cG,
    const Grid3& grid) const {
    return std::make_unique<Quadrilateral4<math::Real>>(this->getId(),
        this->vertexToUnstructured(cG, grid).data(),
        this->getLayer(),
        this->getModel());
}

} 

typedef element::Quadrilateral4Base         Qua4;
typedef element::Quadrilateral4<math::Real> QuaR4;
typedef element::Quadrilateral4<math::Int > QuaI4;

} 
} 