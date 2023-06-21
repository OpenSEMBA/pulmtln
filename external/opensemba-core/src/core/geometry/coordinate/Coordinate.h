#pragma once

#include <memory>

#include "core/geometry/Grid.h"
#include "core/math/CartesianVector.h"
#include "core/util/Class.h"
#include "core/util/Identifiable.h"
#include "core/util/Identification.h"

namespace semba {
namespace geometry {
namespace coordinate {

class Base;
typedef util::Identification<Base> Id;

class Base : public virtual util::Class,
             public virtual util::Identifiable<Id> {
public:
    Base() = default;
    virtual ~Base() = default;
    
    virtual std::unique_ptr<Base> clone() const = 0;

    virtual bool operator==(const Base& rhs) const;
    virtual bool operator!=(const Base& rhs) const;
};

template <class T, std::size_t D>
class Coordinate : public virtual Base,
                   public virtual math::CartesianVector<T,D> {
public:
    Coordinate() = default;
    Coordinate(const Id id_, const math::CartesianVector<T,D>& pos);
    explicit Coordinate(const math::CartesianVector<T,D>& pos);
    Coordinate(const Coordinate& rhs);
    virtual ~Coordinate() = default;

    virtual std::unique_ptr<Base> clone() const override {
        return std::make_unique<Coordinate<T, D>>(*this);
    }

    Coordinate& operator=(const Coordinate& rhs);

    bool operator==(const Base& rhs) const override;
    bool operator!=(const Base& rhs) const override;

    virtual bool isStructured(const Grid<D>&,
                              const math::Real = Grid<D>::tolerance) const;

    math::CartesianVector<T,D>&       pos()       { return *this; }
    const math::CartesianVector<T,D>& pos() const { return *this; }

    virtual Coordinate<math::Int ,D>* toStructured  (const Grid<D>&) const;
    virtual Coordinate<math::Real,D>* toUnstructured(const Grid<D>&) const;

};


template<class T, std::size_t D>
Coordinate<T, D>::Coordinate(const Id id,
    const math::CartesianVector<T, D>& pos)
    : Identifiable<Id>(id),
    math::CartesianVector<T, D>(pos) {

}

template<class T, std::size_t D>
Coordinate<T, D>::Coordinate(const math::CartesianVector<T, D>& pos)
    : math::CartesianVector<T, D>(pos) {

}

template<class T, std::size_t D>
Coordinate<T, D>::Coordinate(const Coordinate& rhs)
    : Identifiable<Id>(rhs),
    math::CartesianVector<T, D>(rhs) {

}

template<class T, std::size_t D>
Coordinate<T, D>& Coordinate<T, D>::operator=(const Coordinate& rhs)
{
    setId(rhs.getId());
    math::CartesianVector<T, D>::operator=(rhs);
    return *this;
}

template<class T, std::size_t D>
bool Coordinate<T, D>::operator==(const Base& rhs) const {
    if (!Base::operator==(rhs)) {
        return false;
    }
    const Coordinate<T, D>* rhsPtr = rhs.castTo<Coordinate<T, D> >();
    bool res = true;
    res &= (this->getId() == rhsPtr->getId());
    res &= (this->pos() == rhsPtr->pos());
    return res;
}

template<class T, std::size_t D>
bool Coordinate<T, D>::operator!=(const Base& rhs) const {
    return Base::operator!=(rhs);
}

template<class T, std::size_t D>
bool Coordinate<T, D>::isStructured(const Grid<D>& grid,
    const math::Real tol) const {
    if (!grid.isCell(*this, tol)) {
        return false;
    }
    return true;
}

template<class T, std::size_t D>
Coordinate<math::Int, D>* Coordinate<T, D>::toStructured(
    const Grid<D>& grid) const {
    return new Coordinate<math::Int, D>(this->getId(), grid.getCell(*this));
}

template<class T, std::size_t D>
Coordinate<math::Real, D>* Coordinate<T, D>::toUnstructured(
    const Grid<D>& grid) const {
    return new Coordinate<math::Real, D>(this->getId(), grid.getPos(*this));
}


} 

typedef coordinate::Id                       CoordId;
typedef coordinate::Base                     Coord;
typedef coordinate::Coordinate<math::Real,3> CoordR3;
typedef coordinate::Coordinate<math::Int ,3> CoordI3;

} 
} 


