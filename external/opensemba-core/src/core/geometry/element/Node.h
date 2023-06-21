

#pragma once

#include "Element.h"

namespace semba {
namespace geometry {
namespace element {

class NodeBase : public virtual Base {
public:
	static const std::size_t sizeOfCoordinates = 1;

    NodeBase() = default;
    virtual ~NodeBase() = default;

    inline std::size_t numberOfCoordinates() const { return 1; }
    inline std::size_t numberOfFaces   () const { return 1; }
    inline std::size_t numberOfVertices() const { return 1; }
    inline std::size_t numberOfSideVertices   (const std::size_t f = 0) const {
        return 1;
    }
    inline std::size_t numberOfSideCoordinates(const std::size_t f = 0) const {
        return 1;
    }
};

template<class T>
class Node : public virtual Element<T>,
             public virtual NodeBase {
public:
    Node() = default;
    Node(const Id id,
         const std::array<const coordinate::Coordinate<T, 3>*, 1> v,
         const Layer* lay = nullptr,
         const Model* mat = nullptr);
    Node(
        const Id id,
        const coordinate::Coordinate<T, 3>* v[1],
        const Layer* lay = nullptr,
        const Model* mat = nullptr);

    Node(const Node<T>& rhs) = default;
    Node& operator=(const Node<T>& rhs) = default;
    Node(Node<T>&&) = default;
    Node& operator=(Node&&) = default;
    virtual ~Node() = default;
    
    virtual std::unique_ptr<Base> clone() const override {
        return std::make_unique<Node>(*this);
    }

    bool isStructured(const Grid3&, const math::Real = Grid3::tolerance) const override;

    const coordinate::Coordinate<T,3>* getV(const std::size_t i) const override;
    const coordinate::Coordinate<T,3>* getSideV(
        const std::size_t f, const std::size_t i) const override;

    const coordinate::Coordinate<T,3>* getVertex(const std::size_t i) const override;
    const coordinate::Coordinate<T,3>* getSideVertex(
            const std::size_t f, const std::size_t i) const override;

    void setV(const std::size_t i, const coordinate::Coordinate<T,3>* coord) override;

    std::unique_ptr<ElemI> toStructured(const CoordI3Group&, const Grid3&,
                        const math::Real = Grid3::tolerance) const override;
    std::unique_ptr<ElemR> toUnstructured(const CoordR3Group&, const Grid3&) const override;

private:
    std::array<const coordinate::Coordinate<T,3>*, 1> v_;
};

template<class T>
Node<T>::Node(
    const Id id, std::array<const coordinate::Coordinate<T, 3>*,1> v,
    const Layer* lay, const Model* mat): Identifiable<Id>(id),
    Elem(lay, mat),
    v_{v}
{
}

template<class T>
Node<T>::Node(
    const Id id,
    const coordinate::Coordinate<T, 3>* v[1],
    const Layer* lay,
    const Model* mat)
{
    std::array<const coordinate::Coordinate<T, 3>*, 1> vArr;
    std::copy(v, v+1, vArr.begin());
    *this = Node<T>{id, vArr, lay, mat};
}

template<class T>
bool Node<T>::isStructured(const Grid3& grid, const math::Real tol) const 
{
    if (!this->vertexInCell(grid, tol)) {
        return false;
    }
    if (!this->getBound().isPoint()) {
        return false;
    }
    if (!this->vertexInBound()) {
        return false;
    }
    return true;
}

template<class T>
const coordinate::Coordinate<T, 3>* Node<T>::getV(const std::size_t i) const 
{
    return v_[i];
}

template<class T>
const coordinate::Coordinate<T, 3>* Node<T>::getSideV(
    const std::size_t f,
    const std::size_t i) const {
    assert(f == 0 && i == 0);
    return v_[i];
}

template<class T>
const coordinate::Coordinate<T, 3>* Node<T>::getVertex(
    const std::size_t i) const {
    assert(i == 0);
    return v_[i];
}

template<class T>
const coordinate::Coordinate<T, 3>* Node<T>::getSideVertex(
    const std::size_t f,
    const std::size_t i) const {
    assert(f == 0 && i == 0);
    return v_[i];
}

template<class T>
void Node<T>::setV(const std::size_t i,
    const coordinate::Coordinate<T, 3>* coord) {
    assert(i < numberOfCoordinates());
    v_[i] = coord;
}

template<class T>
std::unique_ptr<ElemI> Node<T>::toStructured(
    const CoordI3Group& cG,
    const Grid3& grid,
    const math::Real tol) const {
    return std::make_unique<Node<math::Int>>(this->getId(),
        this->vertexToStructured(cG, grid, tol).data(),
        this->getLayer(),
        this->getModel());
}

template<class T>
std::unique_ptr<ElemR> Node<T>::toUnstructured(
    const CoordR3Group& cG,
    const Grid3& grid) const {
    return std::make_unique<Node<math::Real>>(this->getId(),
        this->vertexToUnstructured(cG, grid).data(),
        this->getLayer(),
        this->getModel());
}

typedef Node<math::Real> NodR;
typedef Node<math::Int > NodI;

} 

typedef element::NodeBase         Nod;
typedef element::Node<math::Real> NodR;
typedef element::Node<math::Int > NodI;

} 
} 