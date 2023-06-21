#pragma once

#include "Line.h"

namespace semba {
namespace geometry {
namespace element {

class Line2Base : public virtual LineBase {
public:
    static const std::size_t sizeOfCoordinates = 2;
    virtual ~Line2Base() = default;

    inline std::size_t numberOfCoordinates() const { return sizeOfCoordinates; }
};

template<class T>
class Line2 : public virtual Line<T>,
              public virtual Line2Base {
public:
    Line2() = default;
    Line2(const Id id,
		std::array<const coordinate::Coordinate<T, 3>*,2> v,
		const Layer* lay = nullptr,
		const Model* mat = nullptr);
    Line2(const Id id,
        const coordinate::Coordinate<T, 3>* v[2],
        const Layer* lay = nullptr,
        const Model* mat = nullptr);
    Line2(std::array<const coordinate::Coordinate<T, 3>*, 2>);
    Line2(coordinate::Group<coordinate::Coordinate<T,3> >&,
          const Box<T,3>& box);

    Line2(const Line2<T>& rhs) = default;
    virtual ~Line2() = default;
    
    virtual std::unique_ptr<Base> clone() const override {
        return std::make_unique<Line2>(*this);
    }

    bool isStructured(const Grid3&, const math::Real = Grid3::tolerance) const override;

    const coordinate::Coordinate<T,3>* getV    (const std::size_t i) const override;
    const coordinate::Coordinate<T,3>* getSideV(const std::size_t f,
                                                const std::size_t i) const override;

    const coordinate::Coordinate<T,3>* getVertex    (const std::size_t i) const override;
    const coordinate::Coordinate<T,3>* getSideVertex(
            const std::size_t f,
            const std::size_t i) const override;

    void setV(const std::size_t i, const coordinate::Coordinate<T,3>* coord) override;

    std::unique_ptr<ElemI> toStructured(const CoordI3Group&,
        const Grid3&,
        const math::Real = Grid3::tolerance) const override;
    std::unique_ptr<ElemR> toUnstructured(const CoordR3Group&,
                          const Grid3&) const override;

    std::vector<std::unique_ptr<const Line2<T>>> splitByMiddle() const;

private:
    static const math::simplex::Line<1> lin;
    
    std::array<const coordinate::Coordinate<T,3>*,2> v_;

    void setCoordinates(std::array<const coordinate::Coordinate<T, 3>*, 2>);
    void setCoordinates(coordinate::Group<coordinate::Coordinate<T,3> >&,
                        const Box<T,3>& box);
};


template<class T>
const math::simplex::Line<1> Line2<T>::lin;

template<class T>
Line2<T>::Line2(const Id id,
    std::array<const coordinate::Coordinate<T, 3>*, 2> v,
    const Layer* lay, const Model* mat): 
    Identifiable<Id>(id),
    Elem(lay, mat),
    v_{v}
{
}

template<class T>
Line2<T>::Line2(
    const Id id,
    const coordinate::Coordinate<T, 3>* v[2],
    const Layer* lay,
    const Model* mat)
{
    std::array<const coordinate::Coordinate<T, 3>*, 2> vArr;
    std::copy(v, v+2, vArr.begin());
    *this = Line2<T>(id, vArr, lay, mat);
}

template<class T>
Line2<T>::Line2(std::array<const coordinate::Coordinate<T, 3>*, 2> v) {
    setCoordinates(v);
}

template<class T>
Line2<T>::Line2(coordinate::Group<coordinate::Coordinate<T, 3> >& cG,
    const Box<T, 3>& box) 
{
    setCoordinates(cG, box);
}

template<class T>
bool Line2<T>::isStructured(const Grid3& grid, const math::Real tol) const {
    if (!this->vertexInCell(grid, tol)) {
        return false;
    }
    if (!this->getBound().isLine()) {
        return false;
    }
    if (!this->vertexInBound()) {
        return false;
    }
    return true;
}

template<class T>
const coordinate::Coordinate<T, 3>* Line2<T>::getV(const std::size_t i) const {
    return v_[i];
}

template<class T>
const coordinate::Coordinate<T, 3>* Line2<T>::getSideV(
    const std::size_t f,
    const std::size_t i) const {
    return v_[i];
}

template<class T>
const coordinate::Coordinate<T, 3>* Line2<T>::getVertex(
    const std::size_t i) const {
    return v_[i];
}

template<class T>
const coordinate::Coordinate<T, 3>* Line2<T>::getSideVertex(
    const std::size_t f,
    const std::size_t i) const {
    return v_[i];
}

template<class T>
void Line2<T>::setV(const std::size_t i,
    const coordinate::Coordinate<T, 3>* coord) {

    assert(i < numberOfCoordinates());
    v_[i] = coord;
}

template<class T>
void Line2<T>::setCoordinates(std::array<const coordinate::Coordinate<T, 3>*, 2> v) {
    v_ = v;
}

template<class T>
std::unique_ptr<ElemI> Line2<T>::toStructured(
    const CoordI3Group& cG,
    const Grid3& grid, const math::Real tol) const 
{
    return std::make_unique<Line2<math::Int>>(this->getId(),
        this->vertexToStructured(cG, grid, tol).data(),
        this->getLayer(),
        this->getModel()
    );
}

template<class T>
std::unique_ptr<ElemR> Line2<T>::toUnstructured(
    const CoordR3Group& cG,
    const Grid3& grid) const {
    return std::make_unique<Line2<math::Real>>(
        this->getId(),
        this->vertexToUnstructured(cG, grid).data(),
        this->getLayer(),
        this->getModel()
    );
}

template<class T>
std::vector<std::unique_ptr<const Line2<T>>> Line2<T>::splitByMiddle() const {
    auto vertices = this->getVertices();
    const coordinate::Coordinate<T, 3>* left = vertices[0];
    const coordinate::Coordinate<T, 3>* right = vertices[1];

    const coordinate::Coordinate<T, 3>* middleCoordinate = new coordinate::Coordinate<T, 3>{ CoordId(), (*right - *left) / 2 + *left };

    std::vector<std::unique_ptr<const Line2<T>>> result{};

    result.push_back(
        std::make_unique<const Line2<T>>(
            ElemId(),
            std::array<const coordinate::Coordinate<T, 3>*, 2>{ left, middleCoordinate },
            this->getLayer(),
            this->getModel()
        )
    );
    result.push_back(
        std::make_unique<const Line2<T>>(
            ElemId(),
            std::array<const coordinate::Coordinate<T, 3>*, 2>{ middleCoordinate, right },
            this->getLayer(),
            this->getModel()
        )
    );

    return result;
}


} 

typedef element::Line2Base         Lin2;
typedef element::Line2<math::Real> LinR2;
typedef element::Line2<math::Int> LinI2;

} 
} 