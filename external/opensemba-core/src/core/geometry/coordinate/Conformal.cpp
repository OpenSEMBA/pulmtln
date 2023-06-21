#include "Conformal.h"

namespace semba {
namespace geometry {
namespace coordinate {

Conformal::Conformal() {
    dir_ = math::Constants::CartesianAxis::x;
    length_ = 0.0;
}

Conformal::Conformal(const Id id,
                     const math::CVecI3& pos,
                     const math::Constants::CartesianAxis dir,
                     const math::Real length)
:   Identifiable<Id>(id),
    math::CVecI3(pos) {

    dir_    = dir;
    length_ = length;
}

Conformal::Conformal(const math::Constants::CartesianAxis dir,
                     const math::Real length) {

    dir_    = dir;
    length_ = length;
}

Conformal::Conformal(const Conformal& rhs)
:   Identifiable<Id>(rhs),
    math::CVecI3(rhs) {

    dir_    = rhs.dir_;
    length_ = rhs.length_;
}

Conformal::~Conformal() {

}

Conformal& Conformal::operator=(const Conformal& rhs) {
    if (this == &rhs)
        return *this;

    CoordI3::operator=(rhs);
    dir_    = rhs.dir_;
    length_ = rhs.length_;

    return *this;
}

bool Conformal::operator==(const Base& rhs) const {
    if (!Coordinate<math::Int,3>::operator==(rhs)) {
        return false;
    }
    const Conformal* rhsPtr = rhs.castTo<Conformal>();
    bool res = true;
    res &= (this->length_ == rhsPtr->length_);
    res &= (this->dir_ == rhsPtr->dir_);
    return res;
}

CoordR3* Conformal::toUnstructured(const Grid3& grid) const {
    math::CVecR3 pos = grid.getPos(*this);
    if (math::greater(getLength(), 0.0)) {
        math::Int dir = getDir();
        math::Real length = getLength();
        math::CVecI3 cellAux = *this;
        cellAux(dir)++;
        math::CVecR3 posAux = grid.getPos(cellAux);
        math::Real step = posAux(dir)-pos(dir);
        pos(dir) += step*length;
    }
    return new CoordR3(this->getId(), pos);
}

} 
} 
} 
